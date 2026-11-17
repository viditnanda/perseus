// rips_sparse.hpp -- Vietoris-Rips / flag filtrations without n x n storage.
//
// The radius graph is built directly (a uniform grid for low-dimensional point clouds, a streaming scan
// otherwise), and the flag complex up to the requested skeleton is enumerated from sorted adjacency lists
// and handed to the explicit engine. Memory is linear in the number of simplices. This is the right path
// when the graph at the threshold is sparse; dense problems keep the implicit engine with edge collapse.
#pragma once
#include "io.hpp"
#include "complex.hpp"
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <array>

namespace perseus {

struct WEdge { int a, b; double w; };

inline double enclosing_radius_any(const RipsInput& R) {
    if (R.n <= 1) return 0.0;
    if (R.from_edges && R.lower.empty()) return std::numeric_limits<double>::infinity();
    double enc = std::numeric_limits<double>::infinity();
    for (int i = 0; i < R.n; ++i) {
        double mx = 0;
        for (int j = 0; j < R.n && mx < enc; ++j) if (i != j) mx = std::max(mx, R.dist(i, j));
        enc = std::min(enc, mx);
    }
    return enc;
}

// all edges of length <= thr (each once, a > b)
// the grid route needs every coordinate finite and |x| / thr <= 2^40 (see radius_edges)
inline bool grid_safe(const RipsInput& R, double thr) {
    const double lim = std::ldexp(thr, 40);
    for (double x : R.pts) if (!std::isfinite(x) || std::fabs(x) > lim) return false;
    return true;
}
inline std::vector<WEdge> radius_edges(const RipsInput& R, double thr) {
    std::vector<WEdge> out;
    const int n = R.n;
    if (R.from_edges && R.lower.empty()) {
        std::unordered_map<uint64_t, double> best;
        for (auto& e : R.edges) {
            int a = std::get<0>(e), b = std::get<1>(e); double w = std::get<2>(e);
            if (a < b) std::swap(a, b);
            if (!(w <= thr)) continue;
            uint64_t key = (uint64_t)a * (uint64_t)n + (uint64_t)b;
            auto it = best.find(key); if (it == best.end()) best.emplace(key, w); else it->second = std::min(it->second, w);
        }
        for (auto& kv : best) out.push_back({(int)(kv.first / n), (int)(kv.first % n), kv.second});
    } else if (R.lower.empty() && R.pdim > 0 && R.pdim <= 6 && std::isfinite(thr) && thr > 0 && grid_safe(R, thr)) {
        // Uniform grid: neighbours lie in the 3^d surrounding cells. Cells are floor(x / w) with w slightly
        // above thr and |x / w| <= 2^40 (checked by grid_safe), so the rounding of x / w (relative 2^-53, so
        // absolute at most 2^-13) cannot push two coordinates within thr of each other two cells apart; no
        // shift by a minimum, whose subtraction would round when coordinates differ greatly in size.
        const int d = R.pdim;
        const double w = thr * (1 + 1.0 / 1024);
        auto cell = [&](int i, int k) { return (int64_t)std::floor(R.pts[(size_t)i * d + k] / w); };
        auto hkey = [&](const int64_t* c) { uint64_t h = 1469598103934665603ULL; for (int k = 0; k < d; ++k) { h ^= (uint64_t)c[k] + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2); } return h; };
        std::unordered_map<uint64_t, std::vector<int>> grid; grid.reserve(n * 2);
        std::vector<int64_t> c(d), cc(d);
        for (int i = 0; i < n; ++i) { for (int k = 0; k < d; ++k) c[k] = cell(i, k); grid[hkey(c.data())].push_back(i); }
        int64_t n3 = 1; for (int k = 0; k < d; ++k) n3 *= 3;
        for (int i = 0; i < n; ++i) {
            for (int k = 0; k < d; ++k) c[k] = cell(i, k);
            for (int64_t s = 0; s < n3; ++s) {
                int64_t r = s; for (int k = 0; k < d; ++k) { cc[k] = c[k] + r % 3 - 1; r /= 3; }
                auto it = grid.find(hkey(cc.data())); if (it == grid.end()) continue;
                for (int j : it->second) {
                    if (j >= i) continue;
                    // the hash may collide: check the cell exactly
                    bool same = true; for (int k = 0; k < d && same; ++k) same = cell(j, k) == cc[k];
                    if (!same) continue;
                    double w = R.pdist(i, j);
                    if (w <= thr) out.push_back({i, j, w});
                }
            }
        }
    } else {
        for (int a = 1; a < n; ++a) for (int b = 0; b < a; ++b) { double w = R.dist(a, b); if (std::isfinite(w) && w <= thr) out.push_back({a, b, w}); }
    }
    return out;
}

// Whether the automatic engine choice should build the flag complex explicitly. The explicit route costs
// roughly 1.5 us per simplex, the implicit one roughly 30 ns per entry of the n x n distance table, so the
// explicit route is taken only for n >= 2048, average degree <= 64 and at most n^2 / 50 simplices up to
// dimension skel. Simplices are counted (not stored) and the count stops at the budget. On return E holds
// the radius graph when it was built.
inline bool prefer_sparse(const RipsInput& R, std::vector<WEdge>& E, int skel) {
    const int n = R.n;
    if (n < 2048) return false;
    E = radius_edges(R, R.threshold);
    if ((double)E.size() > 32.0 * n) return false;
    const double budget = (double)n * n / 50;
    std::vector<std::vector<int>> up(n);
    for (auto& e : E) up[std::min(e.a, e.b)].push_back(std::max(e.a, e.b));
    for (auto& v : up) std::sort(v.begin(), v.end());
    double count = n + (double)E.size();
    if (count > budget) return false;
    // cliques of size >= 3: extend each clique by common later neighbours
    std::vector<std::vector<int>> cand(skel + 2);
    auto rec = [&](auto&& self, int depth) -> bool {          // cand[depth]: common neighbours of the clique
        if (depth >= skel) return true;                       // clique of size depth + 1 = dimension depth
        for (int w : cand[depth]) {
            auto& nxt = cand[depth + 1]; nxt.clear();
            const auto& uw = up[w];
            std::set_intersection(cand[depth].begin(), cand[depth].end(), uw.begin(), uw.end(), std::back_inserter(nxt));
            count += (double)nxt.size();
            if (count > budget) return false;
            if (!nxt.empty() && !self(self, depth + 1)) return false;
        }
        return true;
    };
    for (int v = 0; v < n; ++v) {
        if (skel < 2) break;
        for (int w : up[v]) {                                 // edge v-w: triangles through common later neighbours
            auto& c2 = cand[2]; c2.clear();
            std::set_intersection(up[v].begin(), up[v].end(), up[w].begin(), up[w].end(), std::back_inserter(c2));
            count += (double)c2.size();
            if (count > budget) return false;
            if (!c2.empty() && !rec(rec, 2)) return false;
        }
    }
    return true;
}

// the flag complex of a weighted graph up to dimension skel (simplex value = max of its edge values and
// vertex births). Cells are numbered by dimension; faces are found through a hash of the vertex tuple.
inline FilteredComplex flag_complex(int n, const std::vector<WEdge>& E, const std::vector<double>& vbirth, int skel) {
    std::vector<std::vector<std::pair<int, double>>> up(n);          // neighbours with larger index
    for (auto& e : E) { int a = std::min(e.a, e.b), b = std::max(e.a, e.b); up[a].push_back({b, e.w}); }
    for (auto& v : up) std::sort(v.begin(), v.end());
    typedef std::array<int32_t, 8> Key;
    struct KH { size_t operator()(const Key& k) const { uint64_t h = 1469598103934665603ULL; for (int x : k) { h ^= (uint64_t)(uint32_t)x + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2); } return h; } };
    if (skel > 7) throw std::runtime_error("flag complex: dimension above 7 not supported");
    // enumerate simplices by dimension: lists of (vertices, value)
    std::vector<std::vector<std::pair<Key, double>>> S(skel + 1);
    for (int v = 0; v < n; ++v) { Key k; k.fill(-1); k[0] = v; S[0].push_back({k, vbirth.empty() ? 0.0 : vbirth[v]}); }
    // cand: (vertex, max edge value to the current clique)
    std::vector<std::vector<std::pair<int, double>>> cand(skel + 2);
    Key cur; cur.fill(-1);
    auto rec = [&](auto&& self, int depth, double val) -> void {
        // cur[0..depth] is a clique of dimension depth with value val; cand[depth] its common upper neighbours
        if (depth == skel) return;
        for (auto& cw : cand[depth]) {
            int u = cw.first; double nv = std::max(val, cw.second);
            cur[depth + 1] = u;
            S[depth + 1].push_back({cur, nv});
            if (depth + 1 < skel) {
                auto& nc = cand[depth + 1]; nc.clear();
                const auto& nu = up[u];
                size_t i = 0, j = 0;
                // candidates after u: intersect cand[depth] (elements > u) with up[u]
                auto& cd = cand[depth];
                while (i < cd.size() && cd[i].first <= u) ++i;
                while (i < cd.size() && j < nu.size()) {
                    if (cd[i].first < nu[j].first) ++i; else if (nu[j].first < cd[i].first) ++j;
                    else { nc.push_back({cd[i].first, std::max(cd[i].second, nu[j].second)}); ++i; ++j; }
                }
                self(self, depth + 1, nv);
            }
            cur[depth + 1] = -1;
        }
    };
    if (skel >= 1) for (int v = 0; v < n; ++v) {
        cur.fill(-1); cur[0] = v;
        cand[0].assign(up[v].begin(), up[v].end());
        rec(rec, 0, vbirth.empty() ? 0.0 : vbirth[v]);
    }
    FilteredComplex C;
    std::unordered_map<Key, int32_t, KH> id; id.reserve(S.size() ? S[0].size() * 4 : 0);
    std::vector<std::pair<int32_t, int64_t>> bd;
    for (int k = 0; k <= skel; ++k) {
        for (auto& sv : S[k]) {
            bd.clear();
            if (k > 0) for (int i = 0; i <= k; ++i) {
                Key f; f.fill(-1); int m = 0; for (int j = 0; j <= k; ++j) if (j != i) f[m++] = sv.first[j];
                bd.push_back({id.at(f), (i & 1) ? -1 : 1});
            }
            id[sv.first] = C.add_cell(k, sv.second, bd);
        }
        std::vector<std::pair<Key, double>>().swap(S[k]);
    }
    return C;
}

} // namespace perseus
