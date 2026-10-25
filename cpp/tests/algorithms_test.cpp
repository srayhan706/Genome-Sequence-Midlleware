// Unit tests for the algorithms and the MPI wire format.
// Each algorithm is checked against a brute-force reference, on small
// hand-picked cases and on inputs large enough to take the OpenMP paths.

#include "algorithms/kmer.hpp"
#include "algorithms/kmp.hpp"
#include "algorithms/smith_waterman.hpp"
#include "mpi_protocol.hpp"

#include <cstdio>
#include <map>
#include <random>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace genomic;

static int g_failures = 0;

#define CHECK(cond)                                                     \
  do {                                                                  \
    if (!(cond)) {                                                      \
      std::fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__,       \
                   __LINE__, #cond);                                    \
      ++g_failures;                                                     \
    }                                                                   \
  } while (0)

static std::string random_dna(std::mt19937& rng, size_t n, const char* alphabet = "ACGT") {
  std::string s(n, 'A');
  size_t k = std::char_traits<char>::length(alphabet);
  std::uniform_int_distribution<size_t> pick(0, k - 1);
  for (auto& c : s) c = alphabet[pick(rng)];
  return s;
}

static std::vector<int> naive_search(const std::string& t, const std::string& p) {
  std::vector<int> out;
  if (p.empty() || p.size() > t.size()) return out;
  for (size_t i = 0; i + p.size() <= t.size(); ++i)
    if (t.compare(i, p.size(), p) == 0) out.push_back(static_cast<int>(i));
  return out;
}

static std::map<std::string, int64_t> naive_kmers(const std::string& s, int k) {
  std::map<std::string, int64_t> out;
  for (size_t i = 0; i + k <= s.size(); ++i) {
    std::string sub = s.substr(i, k);
    bool ok = true;
    for (char& c : sub) {
      c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
      if (c != 'A' && c != 'C' && c != 'G' && c != 'T') ok = false;
    }
    if (ok) ++out[sub];
  }
  return out;
}

static std::map<std::string, int64_t> as_map(const std::unordered_map<std::string, int64_t>& m) {
  return {m.begin(), m.end()};
}

static void test_kmp() {
  // Regression: matches ending at the last character used to be missed.
  CHECK(kmp_search("AACGT", "ACGT") == std::vector<int>({1}));
  CHECK(kmp_search("AAAA", "AA") == std::vector<int>({0, 1, 2}));
  CHECK(kmp_search("ACGT", "ACGT") == std::vector<int>({0}));
  CHECK(kmp_search("ACG", "ACGT").empty());
  CHECK(kmp_search("ACGT", "").empty());

  std::mt19937 rng(1);
  for (int iter = 0; iter < 500; ++iter) {
    std::string t = random_dna(rng, 1 + rng() % 200, "AC");
    std::string p = random_dna(rng, 1 + rng() % 6, "AC");
    CHECK(kmp_search(t, p) == naive_search(t, p));
  }
  // Big enough for the OpenMP path; low-entropy alphabet gives matches that
  // straddle chunk boundaries.
  std::string big = random_dna(rng, (1 << 20) + 12345, "AC");
  for (const char* p : {"A", "ACCA", "AAAAAAAAAA", "CACACACAC"}) {
    CHECK(kmp_search(big, p) == naive_search(big, p));
  }
}

static void test_kmer() {
  auto m = kmer_count("ACGTACGT", 3);
  CHECK(m.size() == 4);
  CHECK(m["ACG"] == 2 && m["CGT"] == 2 && m["GTA"] == 1 && m["TAC"] == 1);
  CHECK(as_map(kmer_count("acgNacg", 3)) == (std::map<std::string, int64_t>{{"ACG", 2}}));
  CHECK(kmer_count("ACG", 4).empty());
  CHECK(kmer_count("ACG", 0).empty());
  CHECK(kmer_count(std::string(40, 'A'), 33).empty());
  CHECK(kmer_count(std::string(40, 'T'), 32)[std::string(32, 'T')] == 9);

  std::mt19937 rng(2);
  for (int iter = 0; iter < 200; ++iter) {
    std::string s = random_dna(rng, 1 + rng() % 300, "ACGTNacgt");
    int k = 1 + static_cast<int>(rng() % 8);
    CHECK(as_map(kmer_count(s, k)) == naive_kmers(s, k));
  }
  // Regression: the OpenMP path double-counted k-mers at chunk boundaries.
  std::string big(1 << 19, 'A');
  CHECK(kmer_count(big, 3)["AAA"] == static_cast<int64_t>(big.size()) - 2);
  std::string mixed = random_dna(rng, (1 << 19) + 77, "ACGTN");
  CHECK(as_map(kmer_count(mixed, 5)) == naive_kmers(mixed, 5));
}

static void test_smith_waterman() {
  int score = -1;
  auto [a, b] = smith_waterman("ACGTACGT", "ACGTACGT", 2, -1, -1, &score);
  CHECK(score == 16 && a == "ACGTACGT" && b == "ACGTACGT");

  auto [c, d] = smith_waterman("TTTTACGTCCCC", "GGACGTGG", 2, -1, -1, &score);
  CHECK(score == 8 && c == "ACGT" && d == "ACGT");

  smith_waterman("AAAA", "CCCC", 2, -1, -1, &score);
  CHECK(score == 0);

  auto [e, f] = smith_waterman("ACGTTACGT", "ACGTACGT", 2, -1, -1, &score);
  CHECK(score == 15 && e.size() == f.size());

#ifdef _OPENMP
  // Wavefront (parallel) and serial fills must give identical results.
  std::mt19937 rng(3);
  std::string x = random_dna(rng, 1500), y = random_dna(rng, 1300);
  int threads = omp_get_max_threads();
  omp_set_num_threads(1);
  int s1 = 0;
  auto serial = smith_waterman(x, y, 2, -1, -2, &s1);
  omp_set_num_threads(std::max(threads, 4));
  int s2 = 0;
  auto parallel = smith_waterman(x, y, 2, -1, -2, &s2);
  omp_set_num_threads(threads);
  CHECK(s1 == s2 && serial == parallel);
#endif
}

static void test_protocol() {
  AlignParams ap{"ACGT", "AGT", 3, -2, -1};
  AlignParams ap2;
  CHECK(parse_align_request(pack_align_request(ap), ap2));
  CHECK(ap2.seq_a == "ACGT" && ap2.seq_b == "AGT" && ap2.match_s == 3 &&
        ap2.mismatch == -2 && ap2.gap == -1);

  std::string packed = pack_kmer_request({"ACGTACGT", 4});
  KmerParams kp;
  CHECK(parse_kmer_request(packed, kp) && kp.sequence == "ACGTACGT" && kp.k == 4);
  CHECK(!parse_kmer_request(packed.substr(0, packed.size() - 1), kp));
  CHECK(!parse_kmer_request(packed.substr(0, 6), kp));

  PatternParams pp;
  CHECK(!parse_pattern_request("\xff\xff\xff\x7f", pp));

  std::vector<int> pos;
  CHECK(parse_pattern_result(pack_pattern_result({1, 5, 9}), pos) &&
        pos == std::vector<int>({1, 5, 9}));

  std::unordered_map<std::string, int64_t> counts{{"AC", 3}, {"GT", 1LL << 40}};
  std::unordered_map<std::string, int64_t> counts2;
  CHECK(parse_kmer_result(pack_kmer_result(counts), counts2) && counts2 == counts);

  int32_t status = -1;
  int64_t us = 0;
  std::string body;
  CHECK(unframe_result(frame_result(0, 1234, "xyz"), &status, &us, &body));
  CHECK(status == 0 && us == 1234 && body == "xyz");
  CHECK(!unframe_result("short", &status, &us, &body));

  TaskType t{};
  std::string payload;
  CHECK(unframe_task(frame_task(TaskType::kKmer, "abc"), &t, &payload));
  CHECK(t == TaskType::kKmer && payload == "abc");
}

int main() {
  test_kmp();
  test_kmer();
  test_smith_waterman();
  test_protocol();
#ifdef _OPENMP
  std::printf("OpenMP: on (%d threads)\n", omp_get_max_threads());
#else
  std::printf("OpenMP: off\n");
#endif
  if (g_failures) {
    std::printf("%d check(s) FAILED\n", g_failures);
    return 1;
  }
  std::printf("all tests passed\n");
  return 0;
}
