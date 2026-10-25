#include "kmer.hpp"
#include <algorithm>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace genomic {

namespace {

using CodeCounts = std::unordered_map<uint64_t, int64_t>;

int base_code(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

// Count k-mers starting in [first_start, last_start). Each k-mer is packed
// 2 bits per base into a uint64 and updated as the window slides; a non-ACGT
// character resets the window, so k-mers containing it are skipped.
void count_range(const std::string& s, int k, int first_start, int last_start,
                 CodeCounts& out) {
    const uint64_t mask = (k == 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
    const int end = std::min(static_cast<int>(s.size()), last_start + k - 1);
    uint64_t code = 0;
    int run = 0;
    for (int i = first_start; i < end; ++i) {
        int c = base_code(s[i]);
        if (c < 0) { run = 0; code = 0; continue; }
        code = ((code << 2) | static_cast<uint64_t>(c)) & mask;
        if (++run >= k) ++out[code];
    }
}

std::string decode(uint64_t code, int k) {
    static const char kBases[] = "ACGT";
    std::string s(static_cast<size_t>(k), 'A');
    for (int i = k - 1; i >= 0; --i, code >>= 2) s[i] = kBases[code & 3];
    return s;
}

}  // namespace

std::unordered_map<std::string, int64_t> kmer_count(const std::string& sequence, int k) {
    std::unordered_map<std::string, int64_t> total;
    const int n = static_cast<int>(sequence.size());
    if (k <= 0 || k > kMaxKmerK || k > n) return total;
    const int num_starts = n - k + 1;

    CodeCounts counts;
#ifdef _OPENMP
    const int nthreads = omp_get_max_threads();
    if (n >= (1 << 18) && nthreads > 1) {
        std::vector<CodeCounts> local(static_cast<size_t>(nthreads));
        const int chunk = (num_starts + nthreads - 1) / nthreads;
#pragma omp parallel for schedule(static) num_threads(nthreads)
        for (int t = 0; t < nthreads; ++t) {
            int first = t * chunk;
            int last = std::min(num_starts, first + chunk);
            if (first < last) count_range(sequence, k, first, last, local[t]);
        }
        for (const auto& m : local)
            for (const auto& p : m) counts[p.first] += p.second;
    } else
#endif
    {
        count_range(sequence, k, 0, num_starts, counts);
    }

    total.reserve(counts.size());
    for (const auto& p : counts) total.emplace(decode(p.first, k), p.second);
    return total;
}

} // namespace genomic
