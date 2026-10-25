#include "mpi_protocol.hpp"
#include <cstring>
#include <cstddef>

namespace genomic {

namespace {

template <typename T>
void append_pod(std::string& out, T v) {
  out.append(reinterpret_cast<const char*>(&v), sizeof(T));
}

template <typename T>
bool read_pod(const char*& p, const char* end, T* v) {
  if (end - p < static_cast<ptrdiff_t>(sizeof(T))) return false;
  std::memcpy(v, p, sizeof(T));
  p += sizeof(T);
  return true;
}

void append_len_str(std::string& out, const std::string& s) {
  append_pod<uint32_t>(out, static_cast<uint32_t>(s.size()));
  out.append(s);
}

bool read_len_str(const char*& p, const char* end, std::string& s) {
  uint32_t n = 0;
  if (!read_pod(p, end, &n)) return false;
  if (end - p < static_cast<ptrdiff_t>(n)) return false;
  s.assign(p, p + n);
  p += n;
  return true;
}

}  // namespace

const char* task_name(TaskType t) {
  switch (t) {
    case TaskType::kAlign: return "align";
    case TaskType::kPattern: return "pattern";
    case TaskType::kKmer: return "kmer";
    case TaskType::kShutdown: return "shutdown";
  }
  return "unknown";
}

std::string pack_align_request(const AlignParams& p) {
  std::string out;
  out.reserve(20 + p.seq_a.size() + p.seq_b.size());
  append_len_str(out, p.seq_a);
  append_len_str(out, p.seq_b);
  append_pod(out, p.match_s);
  append_pod(out, p.mismatch);
  append_pod(out, p.gap);
  return out;
}

std::string pack_pattern_request(const PatternParams& p) {
  std::string out;
  out.reserve(8 + p.text.size() + p.pattern.size());
  append_len_str(out, p.text);
  append_len_str(out, p.pattern);
  return out;
}

std::string pack_kmer_request(const KmerParams& p) {
  std::string out;
  out.reserve(8 + p.sequence.size());
  append_len_str(out, p.sequence);
  append_pod(out, p.k);
  return out;
}

bool parse_align_request(const std::string& payload, AlignParams& out) {
  const char* p = payload.data();
  const char* end = p + payload.size();
  return read_len_str(p, end, out.seq_a) && read_len_str(p, end, out.seq_b) &&
         read_pod(p, end, &out.match_s) && read_pod(p, end, &out.mismatch) &&
         read_pod(p, end, &out.gap);
}

bool parse_pattern_request(const std::string& payload, PatternParams& out) {
  const char* p = payload.data();
  const char* end = p + payload.size();
  return read_len_str(p, end, out.text) && read_len_str(p, end, out.pattern);
}

bool parse_kmer_request(const std::string& payload, KmerParams& out) {
  const char* p = payload.data();
  const char* end = p + payload.size();
  return read_len_str(p, end, out.sequence) && read_pod(p, end, &out.k);
}

std::string pack_align_result(int score, const std::string& a, const std::string& b) {
  std::string out;
  out.reserve(12 + a.size() + b.size());
  append_pod<int32_t>(out, score);
  append_len_str(out, a);
  append_len_str(out, b);
  return out;
}

std::string pack_pattern_result(const std::vector<int>& positions) {
  std::string out;
  out.reserve(4 + positions.size() * 4);
  append_pod<uint32_t>(out, static_cast<uint32_t>(positions.size()));
  for (int pos : positions) append_pod<int32_t>(out, pos);
  return out;
}

std::string pack_kmer_result(const std::unordered_map<std::string, int64_t>& counts) {
  std::string out;
  append_pod<uint32_t>(out, static_cast<uint32_t>(counts.size()));
  for (const auto& kv : counts) {
    append_len_str(out, kv.first);
    append_pod<int64_t>(out, kv.second);
  }
  return out;
}

bool parse_align_result(const std::string& payload, int& score,
                        std::string& aligned_a, std::string& aligned_b) {
  const char* p = payload.data();
  const char* end = p + payload.size();
  int32_t s = 0;
  if (!read_pod(p, end, &s)) return false;
  score = s;
  return read_len_str(p, end, aligned_a) && read_len_str(p, end, aligned_b);
}

bool parse_pattern_result(const std::string& payload, std::vector<int>& positions) {
  const char* p = payload.data();
  const char* end = p + payload.size();
  uint32_t n = 0;
  if (!read_pod(p, end, &n)) return false;
  if (end - p < static_cast<ptrdiff_t>(n) * 4) return false;
  positions.resize(n);
  for (uint32_t i = 0; i < n; ++i) {
    int32_t v = 0;
    read_pod(p, end, &v);
    positions[i] = v;
  }
  return true;
}

bool parse_kmer_result(const std::string& payload,
                       std::unordered_map<std::string, int64_t>& counts) {
  const char* p = payload.data();
  const char* end = p + payload.size();
  uint32_t n = 0;
  if (!read_pod(p, end, &n)) return false;
  counts.clear();
  counts.reserve(n);
  for (uint32_t i = 0; i < n; ++i) {
    std::string kmer;
    int64_t v = 0;
    if (!read_len_str(p, end, kmer) || !read_pod(p, end, &v)) return false;
    counts[kmer] = v;
  }
  return true;
}

std::string frame_task(TaskType type, const std::string& payload) {
  std::string out;
  out.reserve(kTaskHeaderSize + payload.size());
  append_pod<int32_t>(out, static_cast<int32_t>(type));
  out.append(payload);
  return out;
}

bool unframe_task(const std::string& msg, TaskType* type, std::string* payload) {
  const char* p = msg.data();
  const char* end = p + msg.size();
  int32_t t = 0;
  if (!read_pod(p, end, &t)) return false;
  *type = static_cast<TaskType>(t);
  payload->assign(p, end);
  return true;
}

std::string frame_result(int32_t status, int64_t compute_us, const std::string& body) {
  std::string out;
  out.reserve(kResultHeaderSize + body.size());
  append_pod(out, status);
  append_pod(out, compute_us);
  out.append(body);
  return out;
}

bool unframe_result(const std::string& msg, int32_t* status, int64_t* compute_us,
                    std::string* body) {
  const char* p = msg.data();
  const char* end = p + msg.size();
  if (!read_pod(p, end, status) || !read_pod(p, end, compute_us)) return false;
  body->assign(p, end);
  return true;
}

}  // namespace genomic
