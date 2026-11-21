// pc_maps.hpp -- maps induced on mod-p homology by homomorphisms of finite solvable groups, and the persistent
// homology of a chain of such maps (Ellis-King persistence of a normal or subgroup series).
//
// A homomorphism phi: H -> G of pc groups lifts to a chain map between Anick resolutions through G's contracting
// homotopy, Phi_0 = e, Phi_n(c) = h(Phi_{n-1}(d c)), the group coefficients pushed along phi (sylow.hpp,
// chain_map_eps). On F_p (x) P it carries the cycles f0(s) of H's split to cycles of G, and G's cocycles g0^*(t) give
// their coordinates. So phi_*: H_n(H; F_p) -> H_n(G; F_p) is a matrix in the bases dual to the cohomology bases of
// pc_group.hpp, and phi^*: H^n(G; F_p) -> H^n(H; F_p) is its transpose.
//
// For a chain G_1 -> G_2 -> ... -> G_m, r(i, j) = rank of the composite H_n(G_i) -> H_n(G_j) (r(i, i) = dim) is
// the persistence matrix of Ellis and King (Persistent homology of groups, J. Group Theory 14 (2011); arXiv
// 1006.2237), and the barcode follows: the number of bars born at i and alive through j is
//   r(i, j) - r(i - 1, j) - r(i, j + 1) + r(i - 1, j + 1),
// with r = 0 outside 1..m.
#pragma once
#include "perseus/sylow.hpp"

namespace perseus {
namespace pcg {

struct PcSeriesResult {
    int p = 0, maxdim = 0;
    std::vector<std::vector<size_t>> dims;                                  // dims[g][n] = dim H_n(G_g; F_p)
    std::vector<std::vector<std::vector<std::vector<uint32_t>>>> maps;     // maps[i][n]: dim H_n(G_{i+1}) x dim H_n(G_i)
    std::vector<std::vector<std::vector<size_t>>> rank;                     // rank[n][i][j], i <= j
    struct Bar { int birth, death; size_t mult; };                          // alive at birth..death (indices of groups)
    std::vector<std::vector<Bar>> barcode;                                  // barcode[n]
    double seconds = 0;
};

inline size_t rank_mod_p(std::vector<std::vector<i64>> A, i64 p) {
    size_t r = 0, rows = A.size(), cols = rows ? A[0].size() : 0;
    for (size_t c = 0; c < cols && r < rows; ++c) {
        size_t i = r; while (i < rows && A[i][c] % p == 0) ++i;
        if (i == rows) continue;
        std::swap(A[i], A[r]);
        i64 inv = inv_mod(A[r][c] % p, p);
        for (size_t i2 = r + 1; i2 < rows; ++i2) if (A[i2][c] % p) {
            i64 f = A[i2][c] % p * inv % p;
            for (size_t j = c; j < cols; ++j) A[i2][j] = ((A[i2][j] - f * A[r][j]) % p + p) % p;
        }
        ++r;
    }
    return r;
}

// phi_*: H_n(H; F_p) -> H_n(G; F_p) for n = 0..maxdim, in the bases of the two RingCtx; [n][t][s]
inline std::vector<std::vector<std::vector<uint32_t>>> induced_homology_map(RingCtx& XH, RingCtx& XG, const std::vector<int>& phi, int maxdim) {
    const i64 p = XG.p;
    if ((int)phi.size() != XH.G.n) throw std::invalid_argument("pc_maps: the map must give an image for every element");
    for (int x : phi) if (x < 0 || x >= XG.G.n) throw std::invalid_argument("pc_maps: an image is not an element of the target");
    for (int a = 0; a < XH.G.n; ++a) for (int b = 0; b < XH.G.n; ++b)
        if (phi[XH.G.mul(a, b)] != XG.G.mul(phi[a], phi[b])) throw std::invalid_argument("pc_maps: the map is not a homomorphism");
    auto F = chain_map_eps(XH.A, XG.A, phi, maxdim);
    std::vector<std::vector<std::vector<uint32_t>>> M(maxdim + 1);
    for (int n = 0; n <= maxdim; ++n) {
        M[n].assign(XG.dim(n), std::vector<uint32_t>(XH.dim(n), 0));
        for (size_t s = 0; s < XH.dim(n); ++s) {
            std::vector<i64> w(XG.A.C.size(n), 0);
            for (auto& kv : XH.reps[n][s]) for (auto& d : F[n][kv.first]) w[d.first] = (w[d.first] + kv.second * d.second) % p;
            for (size_t t = 0; t < XG.dim(n); ++t) { i64 x = 0; for (auto& kv : XG.basis[n][t]) x = (x + w[kv.first] * kv.second) % p; M[n][t][s] = (uint32_t)x; }
        }
    }
    return M;
}

inline PcSeriesResult pc_series(const std::vector<PcGroup>& groups_in, const std::vector<std::vector<int>>& maps, int maxdim, int p, bool choose = true, const std::vector<PcgsList>& alternatives = {}) {
    auto t0 = std::chrono::steady_clock::now();
    if (p < 2) throw std::invalid_argument("pc_maps: needs a prime");
    for (i64 d = 2; d * d <= p; ++d) if (p % d == 0) throw std::invalid_argument("pc_maps: needs a prime");
    if (maxdim < 0) throw std::invalid_argument("pc_maps: maxdim must be >= 0");
    size_t m = groups_in.size();
    std::vector<PcGroup> groups;                       // the cheapest pcgs of each group (element ids, so the maps, unchanged)
    for (size_t g = 0; g < groups_in.size(); ++g) groups.push_back(choose ? choose_pcgs(groups_in[g], maxdim, p, nullptr, g < alternatives.size() ? alternatives[g] : PcgsList{}) : groups_in[g]);
    if (m == 0 || maps.size() + 1 != m) throw std::invalid_argument("pc_maps: m groups need m - 1 maps");
    PcSeriesResult R; R.p = p; R.maxdim = maxdim;
    std::vector<std::unique_ptr<RingCtx>> X;
    for (auto& G : groups) { X.emplace_back(new RingCtx(G, maxdim, p)); std::vector<size_t> d; for (int n = 0; n <= maxdim; ++n) d.push_back(X.back()->dim(n)); R.dims.push_back(d); }
    for (size_t i = 0; i + 1 < m; ++i) R.maps.push_back(induced_homology_map(*X[i], *X[i + 1], maps[i], maxdim));
    R.rank.resize(maxdim + 1); R.barcode.resize(maxdim + 1);
    for (int n = 0; n <= maxdim; ++n) {
        auto& r = R.rank[n]; r.assign(m, std::vector<size_t>(m, 0));
        for (size_t i = 0; i < m; ++i) {
            r[i][i] = R.dims[i][n];
            std::vector<std::vector<i64>> C(R.dims[i][n], std::vector<i64>(R.dims[i][n], 0));   // composite i -> j, as rows of columns
            for (size_t a = 0; a < C.size(); ++a) C[a][a] = 1;                 // C[s] = image of basis vector s, a column
            for (size_t j = i + 1; j < m; ++j) {
                const auto& Mj = R.maps[j - 1][n];                               // H_n(G_{j-1}) -> H_n(G_j)
                std::vector<std::vector<i64>> C2(C.size(), std::vector<i64>(R.dims[j][n], 0));
                for (size_t s = 0; s < C.size(); ++s) for (size_t u = 0; u < C[s].size(); ++u) if (C[s][u]) for (size_t t = 0; t < R.dims[j][n]; ++t) C2[s][t] = (C2[s][t] + C[s][u] * Mj[t][u]) % p;
                C.swap(C2);
                r[i][j] = rank_mod_p(C, p);
            }
        }
        auto rr = [&](long i, long j) -> long long { return (i < 0 || j >= (long)m) ? 0 : (long long)r[i][j]; };
        for (long i = 0; i < (long)m; ++i) for (long j = i; j < (long)m; ++j) {
            long long k = rr(i, j) - rr(i - 1, j) - rr(i, j + 1) + rr(i - 1, j + 1);
            if (k < 0) throw std::logic_error("pc_maps: negative bar multiplicity");
            if (k) R.barcode[n].push_back({(int)i, (int)j, (size_t)k});
        }
    }
    R.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return R;
}

// ---- parsing the pcseries format -----------------------------------------------------------------------------------
struct PcSeriesData { std::vector<PcGroup> groups; std::vector<std::vector<int>> maps; std::vector<PcgsList> alternatives; };
inline PcSeriesData parse_pcseries(const std::string& text) {
    std::vector<std::string> tok; size_t pos = 0;
    while (pos <= text.size()) {
        size_t e = text.find('\n', pos); if (e == std::string::npos) e = text.size();
        std::string line = text.substr(pos, e - pos); pos = e + 1;
        auto h = line.find('#'); if (h != std::string::npos) line = line.substr(0, h);
        size_t i = 0; while (i < line.size()) { while (i < line.size() && isspace((unsigned char)line[i])) ++i; size_t j = i; while (j < line.size() && !isspace((unsigned char)line[j])) ++j; if (j > i) tok.push_back(line.substr(i, j - i)); i = j; }
    }
    size_t at = 0;
    auto next = [&]() -> std::string { if (at >= tok.size()) throw std::invalid_argument("pcseries: unexpected end of input"); return tok[at++]; };
    auto num = [&]() -> long long { std::string t = next(); try { size_t q; long long v = std::stoll(t, &q); if (q != t.size()) throw 0; return v; } catch (...) { throw std::invalid_argument("pcseries: expected a number, got '" + t + "'"); } };
    auto word = [&](const char* w) { std::string t = next(); if (t != w) throw std::invalid_argument(std::string("pcseries: expected '") + w + "', got '" + t + "'"); };
    PcSeriesData D; word("groups"); long long m = num();
    if (m < 1) throw std::invalid_argument("pcseries: needs at least one group");
    for (long long g = 0; g < m; ++g) {
        int n = (int)num(), k = (int)num();
        if (n < 1 || k < 0) throw std::invalid_argument("pcseries: bad sizes");
        std::vector<int> pcgs(k), rel(k); for (auto& x : pcgs) x = (int)num(); for (auto& x : rel) x = (int)num();
        int e = (int)num();
        std::vector<std::vector<int>> T(n, std::vector<int>(n)); for (auto& r : T) for (auto& x : r) x = (int)num();
        D.groups.emplace_back(T, pcgs, rel, e);
        PcgsList alts;
        while (at < tok.size() && tok[at] == "alternative") {      // alternative k / pcgs / relative orders
            ++at; int ka = (int)num(); std::vector<int> pa(ka), ra(ka); for (auto& x : pa) x = (int)num(); for (auto& x : ra) x = (int)num();
            alts.push_back({pa, ra});
        }
        D.alternatives.push_back(alts);
    }
    word("maps");
    for (long long g = 0; g + 1 < m; ++g) { std::vector<int> im(D.groups[g].n); for (auto& x : im) x = (int)num(); D.maps.push_back(im); }
    if (at != tok.size()) throw std::invalid_argument("pcseries: trailing input");
    return D;
}

}  // namespace pcg
}  // namespace perseus
