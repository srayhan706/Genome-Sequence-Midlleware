#include "smith_waterman.hpp"
#include <algorithm>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace genomic {

std::pair<std::string, std::string> smith_waterman(
    const std::string& a,
    const std::string& b,
    int match_score,
    int mismatch_penalty,
    int gap_penalty,
    int* out_score
) {
    const int M = static_cast<int>(a.size());
    const int N = static_cast<int>(b.size());
    if (M == 0 || N == 0) {
        if (out_score) *out_score = 0;
        return {"", ""};
    }

    std::vector<int> H(static_cast<size_t>(M + 1) * (N + 1), 0);
    auto at = [&](int i, int j) { return static_cast<size_t>(i) * (N + 1) + j; };
    auto fill_cell = [&](int i, int j) {
        int s = (a[i - 1] == b[j - 1]) ? match_score : mismatch_penalty;
        H[at(i, j)] = std::max({
            H[at(i - 1, j - 1)] + s,
            H[at(i - 1, j)] + gap_penalty,
            H[at(i, j - 1)] + gap_penalty,
            0
        });
    };

    bool filled = false;
#ifdef _OPENMP
    // Cell (i,j) depends on its left, upper and upper-left neighbours, so the
    // matrix is cut into tiles and tiles on the same anti-diagonal are filled
    // in parallel (wavefront). Inside a tile cells go row by row, which keeps
    // memory access cache-friendly.
    constexpr int kTile = 256;
    const int BI = (M + kTile - 1) / kTile;
    const int BJ = (N + kTile - 1) / kTile;
    if (static_cast<long long>(M) * N >= (1LL << 20) && BI > 1 && BJ > 1 &&
        omp_get_max_threads() > 1) {
#pragma omp parallel
        for (int d = 0; d < BI + BJ - 1; ++d) {
            const int bi_lo = std::max(0, d - (BJ - 1));
            const int bi_hi = std::min(BI - 1, d);
#pragma omp for schedule(dynamic, 1)
            for (int bi = bi_lo; bi <= bi_hi; ++bi) {
                const int bj = d - bi;
                const int i_end = std::min(M, (bi + 1) * kTile);
                const int j_end = std::min(N, (bj + 1) * kTile);
                for (int i = bi * kTile + 1; i <= i_end; ++i)
                    for (int j = bj * kTile + 1; j <= j_end; ++j)
                        fill_cell(i, j);
            }
            // implicit barrier: the next anti-diagonal needs this one finished
        }
        filled = true;
    }
#endif
    if (!filled) {
        for (int i = 1; i <= M; ++i)
            for (int j = 1; j <= N; ++j)
                fill_cell(i, j);
    }

    // First maximum in row-major order, so serial and parallel runs agree.
    int best = 0, best_i = 0, best_j = 0;
    for (int i = 1; i <= M; ++i) {
        for (int j = 1; j <= N; ++j) {
            int h = H[at(i, j)];
            if (h > best) { best = h; best_i = i; best_j = j; }
        }
    }

    std::string ra, rb;
    int i = best_i, j = best_j;
    while (i > 0 && j > 0 && H[at(i, j)] > 0) {
        int s = (a[i - 1] == b[j - 1]) ? match_score : mismatch_penalty;
        if (H[at(i, j)] == H[at(i - 1, j - 1)] + s) {
            ra += a[i - 1]; rb += b[j - 1];
            --i; --j;
        } else if (H[at(i, j)] == H[at(i - 1, j)] + gap_penalty) {
            ra += a[i - 1]; rb += '-';
            --i;
        } else {
            ra += '-'; rb += b[j - 1];
            --j;
        }
    }
    std::reverse(ra.begin(), ra.end());
    std::reverse(rb.begin(), rb.end());

    if (out_score) *out_score = best;
    return {ra, rb};
}

} // namespace genomic
