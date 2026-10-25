#include "kmp.hpp"
#include <vector>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace genomic {

static std::vector<int> compute_lps(const std::string& pat) {
    int n = static_cast<int>(pat.size());
    std::vector<int> lps(n, 0);
    for (int i = 1, len = 0; i < n; ) {
        if (pat[i] == pat[len]) { ++len; lps[i++] = len; }
        else if (len) len = lps[len - 1];
        else lps[i++] = 0;
    }
    return lps;
}

// Report matches whose start lies in [first_start, last_start). The scan has to
// read m - 1 characters past last_start to complete the final candidates.
static void kmp_range(const std::string& text, const std::string& pattern,
                      const std::vector<int>& lps, int first_start, int last_start,
                      std::vector<int>& out) {
    const int n = static_cast<int>(text.size());
    const int m = static_cast<int>(pattern.size());
    const int scan_end = std::min(n, last_start + m - 1);
    int j = 0;
    for (int i = first_start; i < scan_end; ++i) {
        while (j > 0 && text[i] != pattern[j]) j = lps[j - 1];
        if (text[i] == pattern[j]) ++j;
        if (j == m) {
            out.push_back(i - m + 1);
            j = lps[j - 1];
        }
    }
}

std::vector<int> kmp_search(const std::string& text, const std::string& pattern) {
    std::vector<int> out;
    const int n = static_cast<int>(text.size());
    const int m = static_cast<int>(pattern.size());
    if (m == 0 || n < m) return out;

    const std::vector<int> lps = compute_lps(pattern);
    const int num_starts = n - m + 1;

#ifdef _OPENMP
    const int min_parallel = 1 << 20;
    const int nthreads = omp_get_max_threads();
    if (n >= min_parallel && nthreads > 1) {
        // Each thread owns a contiguous block of start positions; results are
        // concatenated in block order, so the output stays sorted.
        std::vector<std::vector<int>> local(static_cast<size_t>(nthreads));
        const int chunk = (num_starts + nthreads - 1) / nthreads;
#pragma omp parallel for schedule(static) num_threads(nthreads)
        for (int t = 0; t < nthreads; ++t) {
            int first = t * chunk;
            int last = std::min(num_starts, first + chunk);
            if (first < last) kmp_range(text, pattern, lps, first, last, local[t]);
        }
        for (const auto& v : local) out.insert(out.end(), v.begin(), v.end());
        return out;
    }
#endif

    kmp_range(text, pattern, lps, 0, num_starts, out);
    return out;
}

} // namespace genomic
