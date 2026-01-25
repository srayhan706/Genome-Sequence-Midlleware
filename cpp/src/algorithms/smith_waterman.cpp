#include "smith_waterman.hpp"
#include <algorithm>
#include <vector>

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

    std::vector<int> H((M + 1) * (N + 1), 0);
    auto at = [&](int i, int j) { return i * (N + 1) + j; };

    int best = 0, best_i = 0, best_j = 0;
    for (int i = 1; i <= M; ++i) {
        for (int j = 1; j <= N; ++j) {
            int s = (a[i - 1] == b[j - 1]) ? match_score : mismatch_penalty;
            int h = std::max({
                H[at(i - 1, j - 1)] + s,
                H[at(i - 1, j)] + gap_penalty,
                H[at(i, j - 1)] + gap_penalty,
                0
            });
            H[at(i, j)] = h;
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
