// rips_explicit.hpp -- build the Vietoris-Rips complex explicitly (reference path, and flag filtrations with vertex births).
#pragma once
#include "io.hpp"
#include "complex.hpp"

namespace perseus {

// all simplices of dimension <= skel whose diameter is <= threshold (edges of infinite length absent)
inline FilteredComplex rips_explicit(const RipsInput& R, int skel) {
    int n = R.n;
    std::vector<std::vector<int>> nb(n);
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j)
        if (i != j && std::isfinite(R.dist(i, j)) && R.dist(i, j) <= R.threshold) nb[i].push_back(j);
    FilteredComplex C;
    std::map<std::vector<int>, int32_t> id;   // simplex -> cell id
    std::vector<std::vector<std::pair<std::vector<int>, double>>> layer(skel + 1);
    for (int i = 0; i < n; ++i) layer[0].push_back({{i}, R.vertex_birth.empty() ? 0.0 : R.vertex_birth[i]});
    for (int k = 1; k <= skel; ++k)
        for (auto& sd : layer[k - 1]) {
            int last = sd.first.back();
            for (int w : nb[last]) {
                if (w <= last) continue;
                double d = std::max(sd.second, R.vertex_birth.empty() ? 0.0 : R.vertex_birth[w]); bool ok = true;
                for (int v : sd.first) { double x = R.dist(v, w); if (!std::isfinite(x) || x > R.threshold) { ok = false; break; } d = std::max(d, x); }
                if (!ok) continue;
                std::vector<int> s = sd.first; s.push_back(w);
                layer[k].push_back({s, d});
            }
        }
    std::vector<std::pair<int32_t, int64_t>> bd; std::vector<int> f;
    for (int k = 0; k <= skel; ++k) for (auto& sd : layer[k]) {
        bd.clear();
        if (k > 0) for (int i = 0; i <= k; ++i) { f = sd.first; f.erase(f.begin() + i); bd.push_back({id.at(f), (i & 1) ? -1 : 1}); }
        id[sd.first] = C.add_cell(k, sd.second, bd);
    }
    return C;
}

} // namespace perseus
