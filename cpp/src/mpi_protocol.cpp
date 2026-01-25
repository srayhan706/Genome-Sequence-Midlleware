#include "mpi_protocol.hpp"
#include <cstring>
#include <sstream>
#include <algorithm>

namespace genomic {

namespace {
void append_len_str(std::string& out, const std::string& s) {
  uint32_t n = static_cast<uint32_t>(s.size());
  out.append(reinterpret_cast<const char*>(&n), 4);
  out.append(s);
}
void read_len_str(const char*& p, const char* end, std::string& s) {
  if (end - p < 4) return;
  uint32_t n;
  std::memcpy(&n, p, 4);
  p += 4;
  if (end - p < static_cast<ptrdiff_t>(n)) return;
  s.assign(p, p + n);
  p += n;
}
}  // namespace

std::string pack_align_request(const AlignParams& p) {
  std::string out;
  out.reserve(20 + p.seq_a.size() + p.seq_b.size());
  append_len_str(out, p.seq_a);
  append_len_str(out, p.seq_b);
  out.append(reinterpret_cast<const char*>(&p.match_s), 4);
  out.append(reinterpret_cast<const char*>(&p.mismatch), 4);
  out.append(reinterpret_cast<const char*>(&p.gap), 4);
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
  out.append(reinterpret_cast<const char*>(&p.k), 4);
  return out;
}

bool parse_align_request(const std::string& payload, AlignParams& out) {
  const char* p = payload.data();
  const char* end = p + payload.size();
  read_len_str(p, end, out.seq_a);
  read_len_str(p, end, out.seq_b);
  if (end - p < 12) return false;
  std::memcpy(&out.match_s, p, 4); p += 4;
  std::memcpy(&out.mismatch, p, 4); p += 4;
  std::memcpy(&out.gap, p, 4);
  return true;
}

bool parse_pattern_request(const std::string& payload, PatternParams& out) {
  const char* p = payload.data();
  const char* end = p + payload.size();
  read_len_str(p, end, out.text);
  read_len_str(p, end, out.pattern);
  return true;
}

bool parse_kmer_request(const std::string& payload, KmerParams& out) {
  const char* p = payload.data();
  const char* end = p + payload.size();
  read_len_str(p, end, out.sequence);
  if (end - p < 4) return false;
  std::memcpy(&out.k, p, 4);
  return true;
}

std::string pack_align_result(int score, const std::string& a, const std::string& b) {
  std::string out;
  out.reserve(12 + a.size() + b.size());
  out.append(reinterpret_cast<const char*>(&score), 4);
  append_len_str(out, a);
  append_len_str(out, b);
  return out;
}

std::string pack_pattern_result(const std::vector<int>& positions) {
  std::string out;
  uint32_t n = static_cast<uint32_t>(positions.size());
  out.append(reinterpret_cast<const char*>(&n), 4);
  out.reserve(4 + n * 4);
  for (int pos : positions) {
    out.append(reinterpret_cast<const char*>(&pos), 4);
  }
  return out;
}

std::string pack_kmer_result(const std::unordered_map<std::string, int64_t>& counts) {
  std::string out;
  uint32_t n = static_cast<uint32_t>(counts.size());
  out.append(reinterpret_cast<const char*>(&n), 4);
  for (const auto& kv : counts) {
    append_len_str(out, kv.first);
    out.append(reinterpret_cast<const char*>(&kv.second), 8);
  }
  return out;
}

bool parse_align_result(const std::string& payload, int& score,
                        std::string& aligned_a, std::string& aligned_b) {
  const char* p = payload.data();
  const char* end = p + payload.size();
  if (end - p < 4) return false;
  std::memcpy(&score, p, 4); p += 4;
  read_len_str(p, end, aligned_a);
  read_len_str(p, end, aligned_b);
  return true;
}

bool parse_pattern_result(const std::string& payload, std::vector<int>& positions) {
  const char* p = payload.data();
  const char* end = p + payload.size();
  if (end - p < 4) return false;
  uint32_t n;
  std::memcpy(&n, p, 4); p += 4;
  if (end - p < static_cast<ptrdiff_t>(n * 4)) return false;
  positions.resize(n);
  for (uint32_t i = 0; i < n; ++i) {
    std::memcpy(&positions[i], p, 4); p += 4;
  }
  return true;
}

bool parse_kmer_result(const std::string& payload,
                       std::unordered_map<std::string, int64_t>& counts) {
  const char* p = payload.data();
  const char* end = p + payload.size();
  if (end - p < 4) return false;
  uint32_t n;
  std::memcpy(&n, p, 4); p += 4;
  counts.clear();
  for (uint32_t i = 0; i < n; ++i) {
    std::string kmer;
    read_len_str(p, end, kmer);
    if (end - p < 8) return false;
    int64_t v;
    std::memcpy(&v, p, 8); p += 8;
    counts[kmer] = v;
  }
  return true;
}

}  // namespace genomic
