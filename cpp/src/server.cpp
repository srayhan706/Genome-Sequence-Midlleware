#include "server.hpp"
#include "dispatch.hpp"
#include "mpi_protocol.hpp"
#include "genomic.grpc.pb.h"
#include "genomic.pb.h"
#include <grpcpp/grpcpp.h>
#include <grpcpp/health_check_service_interface.h>
#include <unordered_map>
#include <mutex>
#include <sstream>
#include <chrono>
#include <iostream>

namespace genomic {

namespace {

static constexpr size_t kMaxCacheEntries = 4096;
std::unordered_map<std::string, std::string> g_cache;
std::mutex g_cache_mutex;

std::string cache_key_align(const std::string& a, const std::string& b,
                            int match_s, int mismatch, int gap) {
  std::ostringstream os;
  os << "align:" << a << "|" << b << "|" << match_s << "|" << mismatch << "|" << gap;
  return os.str();
}

std::string cache_key_pattern(const std::string& text, const std::string& pattern) {
  return "pattern:" + text + "|" + pattern;
}

std::string cache_key_kmer(const std::string& seq, int k) {
  std::ostringstream os;
  os << "kmer:" << seq << "|" << k;
  return os.str();
}

bool cache_get(const std::string& key, std::string* value) {
  std::lock_guard<std::mutex> lock(g_cache_mutex);
  auto it = g_cache.find(key);
  if (it == g_cache.end()) return false;
  *value = it->second;
  return true;
}

void cache_put(const std::string& key, const std::string& value) {
  std::lock_guard<std::mutex> lock(g_cache_mutex);
  if (g_cache.size() >= kMaxCacheEntries) {
    g_cache.erase(g_cache.begin());
  }
  g_cache[key] = value;
}

class GenomicServiceImpl final : public ::genomic::GenomicService::Service {
 public:
  grpc::Status SequenceAlignment(grpc::ServerContext* /*ctx*/,
                                 const ::genomic::AlignmentRequest* req,
                                 ::genomic::AlignmentResponse* rsp) override {
    std::string key = cache_key_align(
        req->sequence_a(), req->sequence_b(),
        req->match_score(), req->mismatch_penalty(), req->gap_penalty());
    std::string cached;
    if (cache_get(key, &cached)) {
      int score;
      std::string a, b;
      if (parse_align_result(cached, score, a, b)) {
        rsp->set_aligned_a(a);
        rsp->set_aligned_b(b);
        rsp->set_score(score);
        rsp->set_cache_hit(true);
        return grpc::Status::OK;
      }
    }

    AlignParams p;
    p.seq_a = req->sequence_a();
    p.seq_b = req->sequence_b();
    p.match_s = req->match_score();
    p.mismatch = req->mismatch_penalty();
    p.gap = req->gap_penalty();
    int score = 0;
    std::string a, b;
    if (!dispatch_align(p, &score, &a, &b)) {
      return grpc::Status(grpc::StatusCode::INTERNAL, "dispatch failed");
    }
    rsp->set_aligned_a(a);
    rsp->set_aligned_b(b);
    rsp->set_score(score);
    rsp->set_cache_hit(false);
    cache_put(key, pack_align_result(score, a, b));
    return grpc::Status::OK;
  }

  grpc::Status PatternSearch(grpc::ServerContext* /*ctx*/,
                             const ::genomic::PatternSearchRequest* req,
                             ::genomic::PatternSearchResponse* rsp) override {
    std::string key = cache_key_pattern(req->text(), req->pattern());
    std::string cached;
    if (cache_get(key, &cached)) {
      std::vector<int> pos;
      if (parse_pattern_result(cached, pos)) {
        for (int i : pos) rsp->add_positions(i);
        rsp->set_count(static_cast<int32_t>(pos.size()));
        rsp->set_cache_hit(true);
        return grpc::Status::OK;
      }
    }

    PatternParams p;
    p.text = req->text();
    p.pattern = req->pattern();
    std::vector<int> pos;
    if (!dispatch_pattern(p, &pos)) {
      return grpc::Status(grpc::StatusCode::INTERNAL, "dispatch failed");
    }
    for (int i : pos) rsp->add_positions(i);
    rsp->set_count(static_cast<int32_t>(pos.size()));
    rsp->set_cache_hit(false);
    cache_put(key, pack_pattern_result(pos));
    return grpc::Status::OK;
  }

  grpc::Status KmerCount(grpc::ServerContext* /*ctx*/,
                         const ::genomic::KmerCountRequest* req,
                         ::genomic::KmerCountResponse* rsp) override {
    std::string key = cache_key_kmer(req->sequence(), req->k());
    std::string cached;
    if (cache_get(key, &cached)) {
      std::unordered_map<std::string, int64_t> counts;
      if (parse_kmer_result(cached, counts)) {
        for (const auto& kv : counts) (*rsp->mutable_counts())[kv.first] = kv.second;
        rsp->set_cache_hit(true);
        return grpc::Status::OK;
      }
    }

    KmerParams p;
    p.sequence = req->sequence();
    p.k = req->k();
    std::unordered_map<std::string, int64_t> counts;
    if (!dispatch_kmer(p, &counts)) {
      return grpc::Status(grpc::StatusCode::INTERNAL, "dispatch failed");
    }
    for (const auto& kv : counts) (*rsp->mutable_counts())[kv.first] = kv.second;
    rsp->set_cache_hit(false);
    cache_put(key, pack_kmer_result(counts));
    return grpc::Status::OK;
  }

  grpc::Status BatchAlignment(grpc::ServerContext* /*ctx*/,
                              const ::genomic::BatchAlignmentRequest* req,
                              ::genomic::BatchAlignmentResponse* rsp) override {
    int match_s = req->match_score();
    int mismatch = req->mismatch_penalty();
    int gap = req->gap_penalty();
    std::string ref = req->reference();
    const auto& seqs = req->sequences();

    if (!ref.empty()) {
      for (const std::string& s : seqs) {
        AlignParams p;
        p.seq_a = ref;
        p.seq_b = s;
        p.match_s = match_s;
        p.mismatch = mismatch;
        p.gap = gap;
        int score = 0;
        std::string a, b;
        if (!dispatch_align(p, &score, &a, &b)) {
          return grpc::Status(grpc::StatusCode::INTERNAL, "dispatch failed");
        }
        auto* res = rsp->add_results();
        res->set_aligned_a(a);
        res->set_aligned_b(b);
        res->set_score(score);
        res->set_cache_hit(false);
      }
    } else {
      for (size_t i = 0; i + 1 < seqs.size(); i += 2) {
        AlignParams p;
        p.seq_a = seqs[i];
        p.seq_b = seqs[i + 1];
        p.match_s = match_s;
        p.mismatch = mismatch;
        p.gap = gap;
        int score = 0;
        std::string a, b;
        if (!dispatch_align(p, &score, &a, &b)) {
          return grpc::Status(grpc::StatusCode::INTERNAL, "dispatch failed");
        }
        auto* res = rsp->add_results();
        res->set_aligned_a(a);
        res->set_aligned_b(b);
        res->set_score(score);
        res->set_cache_hit(false);
      }
    }
    return grpc::Status::OK;
  }

  grpc::Status HealthCheck(grpc::ServerContext* /*ctx*/,
                           const ::genomic::HealthRequest* /*req*/,
                           ::genomic::HealthResponse* rsp) override {
    rsp->set_worker_count(dispatch_worker_count());
    rsp->set_ok(true);
    return grpc::Status::OK;
  }
};

}  // namespace

int run_grpc_server(const std::string& bind_address) {
  GenomicServiceImpl service;
  grpc::EnableDefaultHealthCheckService(true);
  grpc::ServerBuilder builder;
  builder.AddListeningPort(bind_address, grpc::InsecureServerCredentials());
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
