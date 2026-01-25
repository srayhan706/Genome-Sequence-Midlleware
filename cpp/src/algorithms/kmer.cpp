#include "kmer.hpp"
#include <algorithm>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace genomic {

std::unordered_map<std::string, int64_t> kmer_count(const std::string& sequence, int k) {
    std::unordered_map<std::string, int64_t> total;
    int n = static_cast<int>(sequence.size());
    if (k <= 0 || k > n) return total;

#ifdef _OPENMP
    int nthreads = 1;
#pragma omp parallel
    { nthreads = omp_get_num_threads(); }
    if (n >= 20000 && nthreads > 1) {
        std::vector<std::unordered_map<std::string, int64_t>> local(
            static_cast<size_t>(nthreads));
        int chunk = (n - k + 1 + nthreads - 1) / nthreads;
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            int start = std::max(0, tid * chunk - (k - 1));
            int end = std::min(n - k + 1, (tid + 1) * chunk);
            if (start < 0) start = 0;
            for (int i = start; i < end; ++i) {
                std::string sub = sequence.substr(static_cast<size_t>(i),
                    static_cast<size_t>(k));
                bool ok = true;
                for (char c : sub) {
                    if (c != 'A' && c != 'C' && c != 'G' && c != 'T' &&
                        c != 'a' && c != 'c' && c != 'g' && c != 't') {
                        ok = false; break;
                    }
                }
                if (ok) ++local[tid][sub];
            }
        }
        for (const auto& m : local)
            for (const auto& p : m) total[p.first] += p.second;
        return total;
    }
#endif

    for (int i = 0; i <= n - k; ++i) {
        std::string sub = sequence.substr(static_cast<size_t>(i),
            static_cast<size_t>(k));
        bool ok = true;
        for (char c : sub) {
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T' &&
                c != 'a' && c != 'c' && c != 'g' && c != 't') {
                ok = false; break;
            }
        }
        if (ok) ++total[sub];
    }
    return total;
}

} // namespace genomic
