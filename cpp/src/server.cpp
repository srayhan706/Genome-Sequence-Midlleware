#include "server.hpp"
#include "dispatch.hpp"
#include "mpi_protocol.hpp"
#include "algorithms/kmer.hpp"
#include "genomic.grpc.pb.h"
#include "genomic.pb.h"
#include <grpcpp/grpcpp.h>
#include <grpcpp/health_check_service_interface.h>
#include <unistd.h>
#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <iostream>
#include <list>
#include <mutex>
#include <thread>
#include <unordered_map>

namespace genomic {

namespace {

using Clock = std::chrono::steady_clock;

// Input limits. Smith-Waterman keeps an (M+1)x(N+1) int matrix, so its limit
// is on the product of the lengths.
constexpr size_t kMaxSequenceLen = 20'000'000;
constexpr size_t kMaxAlignCells = 25'000'000;
constexpr int kMaxMessageBytes = 64 * 1024 * 1024;
constexpr size_t kCacheCapacityBytes = 64 * 1024 * 1024;

// Results cache keyed by the packed request (an unambiguous encoding of the
// operation's inputs). Least-recently-used entries are evicted once the total
// size of keys + values exceeds the byte budget.
class LruCache {
 public:
  explicit LruCache(size_t capacity_bytes) : capacity_(capacity_bytes) {}

  bool get(const std::string& key, std::string* value) {
    std::lock_guard<std::mutex> lock(mu_);
    auto it = index_.find(key);
    if (it == index_.end()) {
      ++misses_;
      return false;
    }
    items_.splice(items_.begin(), items_, it->second);
    *value = it->second->second;
    ++hits_;
    return true;
  }

  void put(const std::string& key, const std::string& value) {
    size_t size = key.size() + value.size();
    if (size > capacity_ / 4) return;  // don't let one entry flush the cache
    std::lock_guard<std::mutex> lock(mu_);
    if (index_.count(key)) return;
    items_.emplace_front(key, value);
    index_[key] = items_.begin();
    bytes_ += size;
    while (bytes_ > capacity_ && !items_.empty()) {
      auto& last = items_.back();
      bytes_ -= last.first.size() + last.second.size();
      index_.erase(last.first);
      items_.pop_back();
    }
  }

  void fill_stats(StatsResponse* rsp) {
    std::lock_guard<std::mutex> lock(mu_);
    rsp->set_cache_hits(hits_);
    rsp->set_cache_misses(misses_);
    rsp->set_cache_entries(static_cast<int64_t>(items_.size()));
    rsp->set_cache_bytes(static_cast<int64_t>(bytes_));
    rsp->set_cache_capacity_bytes(static_cast<int64_t>(capacity_));
  }

 private:
  using Item = std::pair<std::string, std::string>;
  std::mutex mu_;
  std::list<Item> items_;
  std::unordered_map<std::string, std::list<Item>::iterator> index_;
  size_t capacity_;
  size_t bytes_ = 0;
  int64_t hits_ = 0;
  int64_t misses_ = 0;
};

// Remove whitespace (pasted FASTA lines) and upper-case, so "acgt\nACGT"
// and "ACGTACGT" are the same query.
std::string clean_sequence(const std::string& s) {
  std::string out;
  out.reserve(s.size());
  for (unsigned char c : s) {
    if (!std::isspace(c)) out.push_back(static_cast<char>(std::toupper(c)));
  }
  return out;
}

grpc::Status invalid(const std::string& msg) {
  return grpc::Status(grpc::StatusCode::INVALID_ARGUMENT, msg);
}

void fill_info(const DispatchInfo& d, JobInfo* info) {
  info->set_worker_rank(d.worker_rank);
  info->set_compute_ms(d.compute_ms);
  info->set_total_ms(d.total_ms);
}

// Proto3 cannot tell "unset" from 0, so all-zero scoring means defaults.
grpc::Status make_align_params(const std::string& a, const std::string& b,
                               int match_s, int mismatch, int gap, AlignParams* p) {
  p->seq_a = clean_sequence(a);
  p->seq_b = clean_sequence(b);
  if (p->seq_a.empty() || p->seq_b.empty()) return invalid("sequences must be non-empty");
  if (p->seq_a.size() * p->seq_b.size() > kMaxAlignCells) {
    return invalid("alignment too large: len(a) * len(b) must be <= " +
                   std::to_string(kMaxAlignCells));
  }
  if (match_s == 0 && mismatch == 0 && gap == 0) {
    match_s = 2;
    mismatch = -1;
    gap = -1;
  }
  if (match_s <= 0) return invalid("match_score must be positive");
  p->match_s = match_s;
  p->mismatch = mismatch;
  p->gap = gap;
  return grpc::Status::OK;
}

class GenomicServiceImpl final : public ::genomic::GenomicService::Service {
 public:
  explicit GenomicServiceImpl(bool mpi_thread_multiple)
      : cache_(kCacheCapacityBytes),
        start_(Clock::now()),
        mpi_thread_multiple_(mpi_thread_multiple) {}

  grpc::Status SequenceAlignment(grpc::ServerContext* /*ctx*/,
                                 const AlignmentRequest* req,
                                 AlignmentResponse* rsp) override {
    RequestScope scope(this);
    AlignParams p;
    auto st = make_align_params(req->sequence_a(), req->sequence_b(), req->match_score(),
                                req->mismatch_penalty(), req->gap_penalty(), &p);
    if (!st.ok()) return st;
    return align(p, rsp);
  }

  grpc::Status PatternSearch(grpc::ServerContext* /*ctx*/,
                             const PatternSearchRequest* req,
                             PatternSearchResponse* rsp) override {
    RequestScope scope(this);
    PatternParams p;
    p.text = clean_sequence(req->text());
    p.pattern = clean_sequence(req->pattern());
    if (p.pattern.empty()) return invalid("pattern must be non-empty");
    if (p.text.size() > kMaxSequenceLen) return invalid("text too long");

    std::string result;
    bool hit = false;
    DispatchInfo d;
    auto st = run_job(TaskType::kPattern, pack_pattern_request(p), &result, &hit, &d);
    if (!st.ok()) return st;
    std::vector<int> pos;
    if (!parse_pattern_result(result, pos)) return internal("bad pattern result");
    rsp->mutable_positions()->Add(pos.begin(), pos.end());
    rsp->set_count(static_cast<int32_t>(pos.size()));
    rsp->set_cache_hit(hit);
    fill_info(d, rsp->mutable_info());
    return grpc::Status::OK;
  }

  grpc::Status KmerCount(grpc::ServerContext* /*ctx*/,
                         const KmerCountRequest* req,
                         KmerCountResponse* rsp) override {
    RequestScope scope(this);
    KmerParams p;
    p.sequence = clean_sequence(req->sequence());
    p.k = req->k();
    if (p.k < 1 || p.k > kMaxKmerK) {
      return invalid("k must be between 1 and " + std::to_string(kMaxKmerK));
    }
    if (p.sequence.size() > kMaxSequenceLen) return invalid("sequence too long");

    std::string result;
    bool hit = false;
    DispatchInfo d;
    auto st = run_job(TaskType::kKmer, pack_kmer_request(p), &result, &hit, &d);
    if (!st.ok()) return st;
    std::unordered_map<std::string, int64_t> counts;
    if (!parse_kmer_result(result, counts)) return internal("bad k-mer result");
    rsp->mutable_counts()->insert(counts.begin(), counts.end());
    rsp->set_cache_hit(hit);
    fill_info(d, rsp->mutable_info());
    return grpc::Status::OK;
  }

  // Either align every sequence against `reference`, or align consecutive
  // pairs (0-1, 2-3, ...). Jobs are spread over all workers concurrently.
  grpc::Status BatchAlignment(grpc::ServerContext* /*ctx*/,
                              const BatchAlignmentRequest* req,
                              BatchAlignmentResponse* rsp) override {
    RequestScope scope(this);
    auto t0 = Clock::now();
    const auto& seqs = req->sequences();

    std::vector<AlignParams> jobs;
    auto add = [&](const std::string& a, const std::string& b) {
      AlignParams p;
      auto st = make_align_params(a, b, req->match_score(), req->mismatch_penalty(),
                                  req->gap_penalty(), &p);
      if (st.ok()) jobs.push_back(std::move(p));
      return st;
    };
    if (!req->reference().empty()) {
      for (const auto& s : seqs) {
        if (auto st = add(req->reference(), s); !st.ok()) return st;
      }
    } else {
      if (seqs.size() % 2 != 0) {
        return invalid("without a reference, sequences must come in pairs");
      }
      for (int i = 0; i + 1 < seqs.size(); i += 2) {
        if (auto st = add(seqs[i], seqs[i + 1]); !st.ok()) return st;
      }
    }

    std::vector<AlignmentResponse> results(jobs.size());
    std::vector<grpc::Status> statuses(jobs.size());
    std::atomic<size_t> next{0};
    auto run = [&] {
      for (size_t i; (i = next.fetch_add(1)) < jobs.size();) {
        statuses[i] = align(jobs[i], &results[i]);
      }
    };
    // One thread per worker keeps every worker busy; dispatch_task hands each
    // job to whichever worker is idle.
    size_t parallel = std::max<size_t>(1, dispatch_worker_count());
    if (dispatch_worker_count() == 0) parallel = std::thread::hardware_concurrency();
    parallel = std::max<size_t>(1, std::min(parallel, jobs.size()));
    std::vector<std::thread> threads;
    for (size_t t = 1; t < parallel; ++t) threads.emplace_back(run);
    run();
    for (auto& t : threads) t.join();

    for (size_t i = 0; i < jobs.size(); ++i) {
      if (!statuses[i].ok()) return statuses[i];
      *rsp->add_results() = std::move(results[i]);
    }
    rsp->set_wall_ms(std::chrono::duration<double, std::milli>(Clock::now() - t0).count());
    return grpc::Status::OK;
  }

  grpc::Status HealthCheck(grpc::ServerContext* /*ctx*/,
                           const HealthRequest* /*req*/,
                           HealthResponse* rsp) override {
    rsp->set_worker_count(dispatch_worker_count());
    rsp->set_ok(true);
    return grpc::Status::OK;
  }

  grpc::Status GetStats(grpc::ServerContext* /*ctx*/,
                        const StatsRequest* /*req*/,
                        StatsResponse* rsp) override {
    for (const auto& w : dispatch_snapshot()) {
      auto* ws = rsp->add_workers();
      ws->set_rank(w.node.rank);
      ws->set_hostname(w.node.hostname);
      ws->set_omp_threads(w.node.omp_threads);
      ws->set_jobs(w.jobs);
      ws->set_failures(w.failures);
      ws->set_active(w.active);
      ws->set_current_task(w.current_task);
      ws->set_compute_ms(w.compute_ms);
      ws->set_bytes_sent(w.bytes_sent);
      ws->set_bytes_received(w.bytes_received);
    }
    rsp->set_requests(requests_.load());
    rsp->set_in_flight(in_flight_.load());
    cache_.fill_stats(rsp);
    rsp->set_uptime_s(std::chrono::duration<double>(Clock::now() - start_).count());
#ifdef _OPENMP
    rsp->set_openmp(true);
#endif
    rsp->set_mpi_thread_multiple(mpi_thread_multiple_);
    char host[256] = {};
    gethostname(host, sizeof(host) - 1);
    rsp->set_coordinator_host(host);
    return grpc::Status::OK;
  }

 private:
  // Counts a request for the lifetime of a handler.
  struct RequestScope {
    explicit RequestScope(GenomicServiceImpl* s) : s_(s) {
      ++s_->requests_;
      ++s_->in_flight_;
    }
    ~RequestScope() { --s_->in_flight_; }
    GenomicServiceImpl* s_;
  };

  static grpc::Status internal(const std::string& msg) {
    return grpc::Status(grpc::StatusCode::INTERNAL, msg);
  }

  // Serve from cache, or dispatch to a worker and cache the packed result.
  grpc::Status run_job(TaskType type, const std::string& payload, std::string* result,
                       bool* cache_hit, DispatchInfo* info) {
    auto t0 = Clock::now();
    std::string key = std::string(1, static_cast<char>(type)) + payload;
    if (cache_.get(key, result)) {
      *cache_hit = true;
      info->worker_rank = 0;
      info->compute_ms = 0;
      info->total_ms = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
      return grpc::Status::OK;
    }
    std::string error;
    if (!dispatch_task(type, payload, result, info, &error)) {
      return internal("worker failed: " + error);
    }
    *cache_hit = false;
    cache_.put(key, *result);
    info->total_ms = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
    return grpc::Status::OK;
  }

  grpc::Status align(const AlignParams& p, AlignmentResponse* rsp) {
    std::string result;
    bool hit = false;
    DispatchInfo d;
    auto st = run_job(TaskType::kAlign, pack_align_request(p), &result, &hit, &d);
    if (!st.ok()) return st;
    int score = 0;
    std::string a, b;
    if (!parse_align_result(result, score, a, b)) return internal("bad alignment result");
    rsp->set_aligned_a(a);
    rsp->set_aligned_b(b);
    rsp->set_score(score);
    rsp->set_cache_hit(hit);
    fill_info(d, rsp->mutable_info());
    return grpc::Status::OK;
  }

  LruCache cache_;
  Clock::time_point start_;
  bool mpi_thread_multiple_;
  std::atomic<int64_t> requests_{0};
  std::atomic<int32_t> in_flight_{0};
};

}  // namespace

int run_grpc_server(const std::string& bind_address, bool mpi_thread_multiple) {
  GenomicServiceImpl service(mpi_thread_multiple);
  grpc::EnableDefaultHealthCheckService(true);
  grpc::ServerBuilder builder;
  builder.AddListeningPort(bind_address, grpc::InsecureServerCredentials());
  builder.SetMaxReceiveMessageSize(kMaxMessageBytes);
  builder.SetMaxSendMessageSize(kMaxMessageBytes);
  builder.RegisterService(&service);
  std::unique_ptr<grpc::Server> server(builder.BuildAndStart());
  if (!server) {
    std::cerr << "Failed to start gRPC server on " << bind_address << std::endl;
    return 1;
  }
  std::cout << "Genomic middleware listening on " << bind_address
            << " (workers: " << dispatch_worker_count() << ")" << std::endl;
  server->Wait();
  return 0;
}

}  // namespace genomic
