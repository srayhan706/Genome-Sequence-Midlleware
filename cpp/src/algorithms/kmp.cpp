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

static void kmp_serial(const std::string& text, const std::string& pattern,
                       const std::vector<int>& lps, int i_start, int i_end,
                       std::vector<int>& out) {
    int n = static_cast<int>(text.size());
    int m = static_cast<int>(pattern.size());
    int i = i_start, j = 0;
    while (i < i_end) {
        if (pattern[j] == text[i]) { ++i; ++j; }
        if (j == m) {
            out.push_back(i - j);
            j = lps[j - 1];
        } else if (i < n && pattern[j] != text[i]) {
            if (j) j = lps[j - 1];
            else ++i;
        }
    }
}

std::vector<int> kmp_search(const std::string& text, const std::string& pattern) {
    std::vector<int> out;
    int n = static_cast<int>(text.size());
    int m = static_cast<int>(pattern.size());
    if (m == 0 || n < m) return out;

    std::vector<int> lps = compute_lps(pattern);

#ifdef _OPENMP
    const int min_parallel = 50000;
    if (n >= min_parallel && m <= 64) {
        int nthreads = 1;
#pragma omp parallel
        { nthreads = omp_get_num_threads(); }
        std::vector<std::vector<int>> local(static_cast<size_t>(nthreads));
        int search_len = n - m + 1;
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            int chunk = (search_len + nthreads - 1) / nthreads;
            int start = tid * chunk;
            int end = std::min(start + chunk, search_len);
            for (int i = start; i < end; ++i) {
                bool match = true;
                for (int j = 0; j < m; ++j)
                    if (text[i + j] != pattern[j]) { match = false; break; }
                if (match) local[tid].push_back(i);
            }
        }
        for (const auto& v : local) out.insert(out.end(), v.begin(), v.end());
        std::sort(out.begin(), out.end());
        return out;
    }
#endif

    kmp_serial(text, pattern, lps, 0, n - m + 1, out);
    return out;
}

} // namespace genomic
