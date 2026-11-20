// sylow.hpp -- homology and mod-p cohomology of any finite group G from its Sylow subgroups (stable elements).
//
// Input (tools/export_sylow.g): for each prime p, a Sylow p-subgroup P (a pc group) and, for every double coset P g P
// with g outside P, the subgroup Q = P meet g P g^-1 with the two maps Q -> P: the inclusion iota and the conjugation
// kappa(y) = g^-1 y g. G itself is never needed, so G need not be solvable.
//
// Cartan-Eilenberg: for n >= 1 and p-local trivial coefficients M, H^n(G; M) is the set of x in H^n(P; M) with
// iota^* x = kappa^* x in H^n(Q; M) for every g. Over R = Z/p^v, |P| = p^v, Hom_R(-, R) is exact and dualises
// homology and cohomology, so H_n(G; R) = coker( sum_g H_n(Q; R) -> H_n(P; R), y -> iota_* y - kappa_* y ). Since
// p^v kills A_n = H_n(G; Z)_(p) for n >= 1, universal coefficients give H_n(G; R) = A_n (+) A_{n-1} as R-modules
// (A_0 = Z contributes nothing to H_1), so the cyclic factors of H_n(G; R), less those of A_{n-1}, are those of A_n.
//
// The maps: a homomorphism phi: Q -> P lifts to a chain map Phi between the Anick resolutions through P's contracting
// homotopy, Phi_0 = e, Phi_n(c) = h(Phi_{n-1}(d c)) with group coefficients pushed along phi. Homology over R is read on
// the split of R (x) P by unit pivots (pc_group.hpp, fp_split with Q = p^k): the survivors carry a residual
// differential with no unit entries, the column transforms f0 carry its cycles into P, and the functionals g0 carry
// cycles of P back. Lengths are sums of Smith valuations of small matrices over R.
//
// Mod-p cohomology ring: the stable subspace of H^*(P; F_p) is a subring. Its minimal generators are found degree by
// degree as in pc_group.hpp: a stable generator gamma is lifted once (as a cocycle on P's Anick chains) and multiplied
// with every stable class; the decomposables of degree n are the span of these products.
#pragma once
#include "perseus/pc_group.hpp"
#include <set>
#include <memory>

namespace perseus {
namespace pcg {

struct SylowCoset { std::vector<int> elements, pcgs, images; };          // ids of P's elements
struct SylowPrime { int p = 0; std::vector<std::vector<int>> mult; std::vector<int> pcgs, rel; int identity = 0; std::vector<SylowCoset> cosets; };
struct SylowData { long long order = 0; std::vector<SylowPrime> primes; };

struct StableRing {
    int p = 0;
    std::vector<size_t> dims, indecomposables;
    std::vector<std::vector<std::vector<uint32_t>>> basis;     // basis[n][i]: a stable class, in the basis of H^n(P; F_p)
    std::vector<std::pair<int, std::vector<uint32_t>>> generators;   // (degree, coordinates in basis[degree])
    // (generator index, b, j) -> coordinates of generator . basis[b][j] in basis[a + b]
    std::map<std::tuple<int, int, int>, std::vector<uint32_t>> products;
    std::vector<size_t> sylow_dims;                            // dim H^n(P; F_p)
};
struct SylowResult {
    long long order = 0; std::vector<int> primes; bool complete = false;
    bool has_homology = false; GroupHomology H; std::vector<std::vector<std::vector<u64>>> ppart;   // ppart[n][i]: p_i-part, prime powers
    bool has_ring = false; StableRing R;
    double seconds = 0;
};

// ---- parsing ---------------------------------------------------------------------------------------------------
inline SylowData parse_sylow(const std::string& text) {
    std::vector<std::string> tok; std::string line; size_t pos = 0;
    while (pos <= text.size()) {
        size_t e = text.find('\n', pos); if (e == std::string::npos) e = text.size();
        line = text.substr(pos, e - pos); pos = e + 1;
        auto h = line.find('#'); if (h != std::string::npos) line = line.substr(0, h);
        size_t i = 0; while (i < line.size()) { while (i < line.size() && isspace((unsigned char)line[i])) ++i; size_t j = i; while (j < line.size() && !isspace((unsigned char)line[j])) ++j; if (j > i) tok.push_back(line.substr(i, j - i)); i = j; }
    }
    size_t at = 0;
    auto next = [&]() -> std::string { if (at >= tok.size()) throw std::invalid_argument("sylow: unexpected end of input"); return tok[at++]; };
    auto num = [&]() -> long long { std::string t = next(); try { size_t q; long long v = std::stoll(t, &q); if (q != t.size()) throw 0; return v; } catch (...) { throw std::invalid_argument("sylow: expected a number, got '" + t + "'"); } };
    auto word = [&](const char* w) { std::string t = next(); if (t != w) throw std::invalid_argument(std::string("sylow: expected '") + w + "', got '" + t + "'"); };
    SylowData D; word("order"); D.order = num();
    while (at < tok.size()) {
        word("prime"); SylowPrime S; S.p = (int)num();
        int n = (int)num(), k = (int)num();
        if (n < 1 || k < 0) throw std::invalid_argument("sylow: bad sizes");
        for (int i = 0; i < k; ++i) S.pcgs.push_back((int)num());
        for (int i = 0; i < k; ++i) S.rel.push_back((int)num());
        S.identity = (int)num();
        S.mult.assign(n, std::vector<int>(n));
        for (auto& r : S.mult) for (auto& x : r) x = (int)num();
        word("cosets"); int m = (int)num();
        for (int c = 0; c < m; ++c) {
            SylowCoset C; int q = (int)num(), kq = (int)num();
            for (int i = 0; i < q; ++i) C.elements.push_back((int)num());
            for (int i = 0; i < kq; ++i) C.pcgs.push_back((int)num());
            for (int i = 0; i < q; ++i) C.images.push_back((int)num());
            S.cosets.push_back(C);
        }
        D.primes.push_back(S);
    }
    return D;
}

// ---- a subgroup of P as a pc group, and the two maps into P ---------------------------------------------------
struct SubgroupMaps { PcGroup Q; std::vector<int> iota, kappa; };
inline SubgroupMaps make_subgroup(const PcGroup& P, const SylowCoset& C) {
    int q = (int)C.elements.size();
    std::vector<int> pos(P.n, -1);
    for (int i = 0; i < q; ++i) {
        if (C.elements[i] < 0 || C.elements[i] >= P.n || pos[C.elements[i]] >= 0) throw std::invalid_argument("sylow: bad subgroup element list");
        pos[C.elements[i]] = i;
    }
    if ((int)C.images.size() != q) throw std::invalid_argument("sylow: the images do not match the subgroup");
    std::vector<std::vector<int>> T(q, std::vector<int>(q));
    for (int i = 0; i < q; ++i) for (int j = 0; j < q; ++j) {
        int x = pos[P.mul(C.elements[i], C.elements[j])];
        if (x < 0) throw std::invalid_argument("sylow: the subgroup is not closed under multiplication");
        T[i][j] = x;
    }
    std::vector<int> pcgs;
    for (int g : C.pcgs) { if (g < 0 || g >= P.n || pos[g] < 0) throw std::invalid_argument("sylow: a pcgs element of Q is not in Q"); pcgs.push_back(pos[g]); }
    if (pos[P.e] < 0) throw std::invalid_argument("sylow: the subgroup has no identity");
    SubgroupMaps S{PcGroup(T, pcgs, {}, pos[P.e]), C.elements, C.images};
    std::vector<char> seen(P.n, 0);
    for (int i = 0; i < q; ++i) {
        int y = S.kappa[i]; if (y < 0 || y >= P.n || seen[y]++) throw std::invalid_argument("sylow: the conjugation map is not injective");
        for (int j = 0; j < q; ++j) if (S.kappa[T[i][j]] != P.mul(S.kappa[i], S.kappa[j])) throw std::invalid_argument("sylow: the conjugation map is not a homomorphism");
    }
    return S;
}

// ---- chain maps -------------------------------------------------------------------------------------------------
// The augmented images eps(Phi_n(c)) on P's chains, n = 0..D, for the chain map over phi: Q -> P.
typedef std::vector<std::pair<int32_t, i64>> CV;
inline std::vector<std::vector<CV>> chain_map_eps(const Anick& AQ, Anick& AP, const std::vector<int>& phi, int D) {
    const u64 NQ = (u64)AQ.G.n, NP = (u64)AP.G.n; const i64 M = AP.M;
    std::vector<Vec> prev(1, Vec{{(u64)AP.G.e, 1}});            // Phi_0 on the single 0-chain
    std::vector<std::vector<CV>> out(D + 1);
    auto aug = [&](const std::vector<Vec>& Ph, std::vector<CV>& o) {
        o.resize(Ph.size());
        for (size_t c = 0; c < Ph.size(); ++c) {
            std::map<int32_t, i64> agg; for (auto& kv : Ph[c]) { i64& s = agg[(int32_t)(kv.first / NP)]; s = (s + kv.second) % M; }
            for (auto& kv : agg) if (kv.second) o[c].push_back(kv);
        }
    };
    aug(prev, out[0]);
    for (int n = 1; n <= D; ++n) {
        std::vector<Vec> cur(AQ.C.size(n));
        for (size_t c = 0; c < AQ.C.size(n); ++c) {
            std::vector<std::pair<u64, i64>> acc;
            for (auto& kv : AQ.d[n][c]) {
                size_t cp = kv.first / NQ; int y = phi[kv.first % NQ];
                for (auto& kw : prev[cp]) acc.push_back({(kw.first / NP) * NP + (u64)AP.G.mul((int)(kw.first % NP), y), mulmod(kv.second, kw.second, M)});
            }
            Vec X; normalize(acc, M, X);
            AP.hvec(n - 1, X, cur[c]);
        }
        prev.swap(cur);
        aug(prev, out[n]);
    }
    return out;
}

// ---- linear algebra over Z/p^k ------------------------------------------------------------------------------------
inline int pval(i64 x, i64 p, int k) { if (!x) return k; int v = 0; while (x % p == 0) { x /= p; ++v; } return v; }
// length of the submodule of (Z/p^k)^nrows spanned by the columns
inline long long span_length(std::vector<CV> cols, size_t nrows, i64 p, int k) {
    for (auto& c : cols) std::sort(c.begin(), c.end());
    cols.erase(std::remove_if(cols.begin(), cols.end(), [](const CV& c) { return c.empty(); }), cols.end());
    if (cols.empty() || nrows == 0) return 0;
    long long l = 0; for (int v : local_elementary_valuations(cols, nrows, p, k)) l += k - v;
    return l;
}
// generators of the kernel of (Z/p^k)^c -> (Z/p^k)^r given by the columns: Smith form with the column transform V;
// the kernel is spanned by p^(k - e_i) V e_i on the pivots and V e_j off them
inline std::vector<std::vector<i64>> local_kernel(const std::vector<CV>& cols, size_t nrows, i64 p, int k) {
    i64 Q = 1; for (int i = 0; i < k; ++i) Q *= p;
    size_t c = cols.size();
    std::vector<std::vector<i64>> A(nrows, std::vector<i64>(c, 0)), V(c, std::vector<i64>(c, 0));
    for (size_t j = 0; j < c; ++j) { V[j][j] = 1; for (auto& kv : cols[j]) A[kv.first][j] = ((kv.second % Q) + Q) % Q; }
    std::vector<int> ev; size_t t = 0;
    while (t < std::min(nrows, c)) {
        int best = k; size_t bi = 0, bj = 0;
        for (size_t i = t; i < nrows && best > 0; ++i) for (size_t j = t; j < c; ++j) { int v = pval(A[i][j], p, k); if (v < best) { best = v; bi = i; bj = j; if (!v) break; } }
        if (best == k) break;
        std::swap(A[t], A[bi]);
        for (auto& r : A) std::swap(r[t], r[bj]);
        for (auto& r : V) std::swap(r[t], r[bj]);
        i64 pe = 1; for (int i = 0; i < best; ++i) pe *= p;
        i64 uinv = inv_mod((A[t][t] / pe) % Q, Q);
        for (size_t j = t + 1; j < c; ++j) if (A[t][j]) {
            i64 f = (A[t][j] / pe) % Q * uinv % Q;
            for (size_t i = 0; i < nrows; ++i) if (A[i][t]) A[i][j] = ((A[i][j] - f * A[i][t]) % Q + Q) % Q;
            for (size_t i = 0; i < c; ++i) if (V[i][t]) V[i][j] = ((V[i][j] - f * V[i][t]) % Q + Q) % Q;
        }
        for (size_t i = t + 1; i < nrows; ++i) A[i][t] = 0;   // row operations: they do not change the kernel
        ev.push_back(best); ++t;
    }
    std::vector<std::vector<i64>> out;
    for (size_t i = 0; i < c; ++i) {
        i64 m = 1;
        if (i < t) { if (ev[i] == 0) continue; for (int r = 0; r < k - ev[i]; ++r) m *= p; }
        std::vector<i64> g(c); bool nz = false;
        for (size_t r = 0; r < c; ++r) { g[r] = V[r][i] * m % Q; if (g[r]) nz = true; }
        if (nz) out.push_back(g);
    }
    return out;
}
// the residual differential d'_n over Z/p^k on the survivors: columns S_n, rows S_{n-1}
inline std::vector<CV> residual(const FpSplit& F, int n) {
    std::vector<CV> cols;
    for (int32_t s : F.S[n]) {
        CV col;
        for (auto& kv : F.E[n].resid[s]) { int32_t r = F.sindex[n - 1][kv.first]; if (r < 0) throw std::logic_error("sylow: residual entry off the survivors"); col.push_back({r, (i64)kv.second}); }
        cols.push_back(col);
    }
    return cols;
}

// ---- the p-part of H_n(G; Z), n = 1..D ---------------------------------------------------------------------------
inline std::vector<std::vector<u64>> stable_homology(const SylowPrime& S, int D, bool stats = false) {
    PcGroup P(S.mult, S.pcgs, S.rel, S.identity);
    const i64 p = S.p; int v = 0; { long long x = P.n; while (x % p == 0) { x /= p; ++v; } if (x != 1) throw std::invalid_argument("sylow: P is not a p-group"); }
    std::vector<std::vector<u64>> A(D + 1);
    if (v == 0) return A;
    i64 M = 1; for (int i = 0; i < v; ++i) M *= p;
    auto t0 = std::chrono::steady_clock::now();
    Anick AP(P, D + 2, M);
    struct Cos {                                         // owns Q: the Anick resolution keeps a reference to its group
        SubgroupMaps S; std::unique_ptr<Anick> AQ; std::vector<std::vector<CV>> delta;
        explicit Cos(SubgroupMaps s) : S(std::move(s)) {}
    };
    std::vector<std::unique_ptr<Cos>> cos;
    for (auto& C : S.cosets) {
        if (C.elements.size() <= 1) continue;
        std::unique_ptr<Cos> c(new Cos(make_subgroup(P, C)));
        c->AQ.reset(new Anick(c->S.Q, D + 1, M));
        auto fi = chain_map_eps(*c->AQ, AP, c->S.iota, D), fk = chain_map_eps(*c->AQ, AP, c->S.kappa, D);
        c->delta.resize(D + 1);
        for (int n = 0; n <= D; ++n) {
            c->delta[n].resize(fi[n].size());
            for (size_t q = 0; q < fi[n].size(); ++q) {
                std::map<int32_t, i64> acc;
                for (auto& kv : fi[n][q]) acc[kv.first] = (acc[kv.first] + kv.second) % M;
                for (auto& kv : fk[n][q]) acc[kv.first] = (acc[kv.first] + M - kv.second) % M;
                for (auto& kv : acc) if (kv.second) c->delta[n][q].push_back(kv);
            }
        }
        cos.push_back(std::move(c));
    }
    if (stats) fprintf(stderr, "  [sylow p=%lld: |P| = %d, %zu double cosets, resolutions and chain maps %.3f s]\n", (long long)p, P.n, cos.size(), std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    // One ring R = Z/p^v. Every exponent of A_n is at most v (|P| kills the p-part), so H_n(G; R) = A_n (+) A_{n-1}
    // as R-modules (universal coefficients, A_0 = Z contributing nothing to H_1). With H_n(G; R) = K / S,
    // K = ker d'_n on P's survivors and S = im d'_{n+1} + the images of (iota - kappa)_*, the number of cyclic factors
    // of order > p^j is l(p^j M) - l(p^(j+1) M), where l(p^j M) = l(p^j K + S) - l(S).
    const int k = v; const i64 Qk = M;
    FpSplit FP = fp_split(AP, D + 1, p, Qk, false);           // P: cocycles and residual only, no transforms
    std::vector<FpSplit> FQ; for (auto& c : cos) FQ.push_back(fp_split(*c->AQ, D, p, Qk));
    if (stats) fprintf(stderr, "  [sylow p=%lld: splits mod p^%d, %.3f s]\n", (long long)p, v, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    for (int n = 1; n <= D; ++n) {
        size_t Nn = FP.S[n].size();
        auto cocP = fp_cocycles(FP, n, AP.C.size(n));
        std::vector<CV> Sg = residual(FP, n + 1);
        for (size_t ci = 0; ci < cos.size(); ++ci) {
            const FpSplit& F = FQ[ci]; const Cos& c = *cos[ci];
            for (auto& y : local_kernel(residual(F, n), F.S[n - 1].size(), p, k)) {
                std::map<int32_t, i64> z;                     // the cycle f0(y) on Q's chains
                for (size_t s2 = 0; s2 < y.size(); ++s2) if (y[s2]) for (auto& kv : F.E[n].trans[F.S[n][s2]]) z[kv.first] = (z[kv.first] + y[s2] * kv.second) % Qk;
                std::vector<i64> w(AP.C.size(n), 0);         // (iota - kappa)_* on P's chains
                for (auto& kv : z) if (kv.second) for (auto& d : c.delta[n][kv.first]) w[d.first] = (w[d.first] + kv.second * d.second) % Qk;
                CV col;                                       // g0: coordinates on P's survivors
                for (size_t t = 0; t < Nn; ++t) { i64 s3 = 0; for (auto& kv : cocP[t]) s3 = (s3 + w[kv.first] * kv.second) % Qk; if (s3) col.push_back({(int32_t)t, s3}); }
                if (!col.empty()) Sg.push_back(col);
            }
        }
        std::vector<CV> Kg;                                   // generators of K = ker d'_n
        for (auto& y : local_kernel(residual(FP, n), FP.S[n - 1].size(), p, k)) { CV col; for (size_t t = 0; t < y.size(); ++t) if (y[t]) col.push_back({(int32_t)t, y[t]}); Kg.push_back(col); }
        long long lS = span_length(Sg, Nn, p, k);
        std::vector<long long> lpM(k + 1, 0);                 // l(p^j M), j = 0..k
        for (int j = 0; j < k; ++j) {
            i64 pj = 1; for (int i = 0; i < j; ++i) pj *= p;
            std::vector<CV> cols(Sg);
            for (auto& g : Kg) { CV c2; for (auto& kv : g) { i64 x = kv.second * pj % Qk; if (x) c2.push_back({kv.first, x}); } if (!c2.empty()) cols.push_back(c2); }
            lpM[j] = span_length(cols, Nn, p, k) - lS;
            if (lpM[j] < 0) throw std::logic_error("sylow: impossible length");
        }
        std::multiset<int> Mexp;                              // exponents of the cyclic factors of H_n(G; R)
        for (int j = 0; j < k; ++j) {
            long long gt_j = lpM[j] - lpM[j + 1];                // factors of exponent > j
            long long gt_j1 = j + 1 < k ? lpM[j + 1] - lpM[j + 2] : 0;
            long long exact = gt_j - gt_j1;                       // factors of exponent exactly j + 1
            if (exact < 0) throw std::logic_error("sylow: inconsistent lengths");
            for (long long i = 0; i < exact; ++i) Mexp.insert(j + 1);
        }
        for (u64 x : A[n - 1]) {                               // remove A_{n-1}
            int e = 0; while (x > 1) { x /= (u64)p; ++e; }
            auto it = Mexp.find(e); if (it == Mexp.end()) throw std::logic_error("sylow: H_n(G; Z/p^v) does not contain A_{n-1}");
            Mexp.erase(it);
        }
        for (int e : Mexp) { u64 q = 1; for (int i = 0; i < e; ++i) q *= (u64)p; A[n].push_back(q); }
        std::sort(A[n].begin(), A[n].end());
    }
    if (stats) fprintf(stderr, "  [sylow p=%lld: homology %.3f s]\n", (long long)p, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    return A;
}

// ---- the mod-p cohomology ring: stable elements ------------------------------------------------------------------
// dense row reduction over F_p; returns the pivot columns, rows reduced in place (reduced row echelon form)
inline std::vector<size_t> rref(std::vector<std::vector<i64>>& R, size_t ncols, i64 p) {
    std::vector<size_t> piv; size_t r = 0;
    for (size_t c = 0; c < ncols && r < R.size(); ++c) {
        size_t i = r; while (i < R.size() && R[i][c] % p == 0) ++i;
        if (i == R.size()) continue;
        std::swap(R[r], R[i]);
        i64 inv = inv_mod(R[r][c] % p, p); for (auto& x : R[r]) x = x * inv % p;
        for (size_t i2 = 0; i2 < R.size(); ++i2) if (i2 != r && R[i2][c] % p) { i64 f = R[i2][c] % p; for (size_t j = 0; j < ncols; ++j) R[i2][j] = ((R[i2][j] - f * R[r][j]) % p + p) % p; }
        piv.push_back(c); ++r;
    }
    R.resize(r);
    return piv;
}
inline StableRing stable_ring(const SylowPrime& S, int D, bool stats = false) {
    PcGroup P(S.mult, S.pcgs, S.rel, S.identity);
    const i64 p = S.p;
    auto t0 = std::chrono::steady_clock::now();
    RingCtx X(P, D, p);
    StableRing R; R.p = (int)p;
    // restriction rows: for each coset and each cycle f0(s) of Q, the class (iota - kappa)_* f0(s) in H_n(P; F_p);
    // x in H^n(P) is stable iff it vanishes on all of them
    std::vector<std::vector<std::vector<i64>>> rows(D + 1);
    for (auto& C : S.cosets) {
        if (C.elements.size() <= 1) continue;
        SubgroupMaps sm = make_subgroup(P, C);
        Anick AQ(sm.Q, D + 1, p);
        FpSplit F = fp_split(AQ, D, p);
        auto fi = chain_map_eps(AQ, X.A, sm.iota, D), fk = chain_map_eps(AQ, X.A, sm.kappa, D);
        for (int n = 1; n <= D; ++n) for (int32_t s : F.S[n]) {
            std::vector<i64> w(X.A.C.size(n), 0);
            for (auto& kv : F.E[n].trans[s]) {
                for (auto& d : fi[n][kv.first]) w[d.first] = (w[d.first] + kv.second * d.second) % p;
                for (auto& d : fk[n][kv.first]) w[d.first] = (w[d.first] + (p - kv.second) * d.second) % p;
            }
            std::vector<i64> row(X.dim(n), 0);
            for (size_t i = 0; i < X.dim(n); ++i) { i64 t = 0; for (auto& kv : X.basis[n][i]) t = (t + w[kv.first] * kv.second) % p; row[i] = t; }
            rows[n].push_back(row);
        }
    }
    // the stable subspace: the null space of the rows, one basis vector per free column (so a stable vector's
    // coordinates in this basis are its entries on the free columns)
    std::vector<std::vector<size_t>> freecols(D + 1);
    R.basis.resize(D + 1);
    for (int n = 0; n <= D; ++n) {
        size_t m = X.dim(n); R.sylow_dims.push_back(m);
        auto piv = rref(rows[n], m, p);
        std::vector<char> isp(m, 0); for (auto c : piv) isp[c] = 1;
        for (size_t f = 0; f < m; ++f) if (!isp[f]) {
            std::vector<uint32_t> x(m, 0); x[f] = 1;
            for (size_t r = 0; r < piv.size(); ++r) x[piv[r]] = (uint32_t)((p - rows[n][r][f] % p) % p);
            R.basis[n].push_back(x); freecols[n].push_back(f);
        }
        R.dims.push_back(R.basis[n].size());
    }
    if (stats) fprintf(stderr, "  [sylow ring p=%lld: stable subspaces %.3f s]\n", (long long)p, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    auto stable_coords = [&](int n, const std::vector<uint32_t>& x) {
        std::vector<uint32_t> c; for (size_t f : freecols[n]) c.push_back(x[f]);
        for (auto& r : rows[n]) { i64 t = 0; for (size_t j = 0; j < x.size(); ++j) t = (t + r[j] * x[j]) % p; if (t) throw std::logic_error("sylow: a product of stable classes is not stable"); }
        return c;
    };
    std::vector<Echelon> decomp; for (int n = 0; n <= D; ++n) decomp.emplace_back(p);   // in stable coordinates
    R.indecomposables.assign(D + 1, 0);
    for (int a = 1; a <= D; ++a) {
        Echelon E = decomp[a];
        std::vector<size_t> gens;
        for (size_t i = 0; i < R.dims[a]; ++i) { SV e1{{(int32_t)i, 1u}}; if (E.insert(e1) >= 0) gens.push_back(i); }
        R.indecomposables[a] = gens.size();
        for (size_t i : gens) {
            int gi = (int)R.generators.size();
            std::vector<uint32_t> ec(R.dims[a], 0); ec[i] = 1; R.generators.push_back({a, ec});
            if (a >= D) continue;
            X.lift(a, X.cocycle(a, R.basis[a][i]), 1, [&](int b, const std::vector<std::vector<uint32_t>>& prods) {
                for (size_t j = 0; j < R.basis[b].size(); ++j) {          // gamma . y = sum_t y_t (gamma . beta_t)
                    std::vector<uint32_t> pr(X.dim(a + b), 0);
                    for (size_t t = 0; t < prods.size(); ++t) if (R.basis[b][j][t]) for (size_t u = 0; u < pr.size(); ++u) pr[u] = (uint32_t)((pr[u] + (i64)R.basis[b][j][t] * prods[t][u]) % p);
                    auto c = stable_coords(a + b, pr);
                    R.products[std::make_tuple(gi, b, (int)j)] = c;
                    SV v; for (size_t t = 0; t < c.size(); ++t) if (c[t]) v.push_back({(int32_t)t, c[t]});
                    decomp[a + b].insert(v);
                }
            });
        }
    }
    if (stats) fprintf(stderr, "  [sylow ring p=%lld: products %.3f s]\n", (long long)p, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    return R;
}

inline SylowResult sylow_compute(const SylowData& Dt, int maxdim, bool homology, int ring_prime) {
    auto t0 = std::chrono::steady_clock::now();
    bool stats = std::getenv("PERSEUS_PC_STATS") != nullptr;
    if (maxdim < 0) throw std::invalid_argument("sylow: maxdim must be >= 0");
    SylowResult res; res.order = Dt.order;
    long long prod = 1; std::set<int> seen;
    for (auto& S : Dt.primes) {
        if (S.p < 2 || seen.count(S.p)) throw std::invalid_argument("sylow: a prime is missing or repeated");
        seen.insert(S.p); res.primes.push_back(S.p);
        long long n = (long long)S.mult.size(); prod *= n;
        if (Dt.order % n || (Dt.order / n) % S.p == 0) throw std::invalid_argument("sylow: P is not a Sylow subgroup of a group of the given order");
    }
    res.complete = prod == Dt.order;
    if (homology) {
        res.has_homology = true;
        res.ppart.assign(maxdim + 1, {});
        std::vector<std::vector<std::vector<u64>>> parts;
        for (auto& S : Dt.primes) parts.push_back(stable_homology(S, maxdim, stats));
        res.H.free.assign(maxdim + 1, 0); res.H.free[0] = 1; res.H.torsion.assign(maxdim + 1, {});
        for (int n = 1; n <= maxdim; ++n) {
            size_t len = 0;
            for (auto& pp : parts) { res.ppart[n].push_back(pp[n]); len = std::max(len, pp[n].size()); }
            std::vector<u64> inv(len, 1);                         // invariant factors, each dividing the next
            for (auto& pp : parts) { std::vector<u64> l(pp[n]); std::sort(l.rbegin(), l.rend()); for (size_t i = 0; i < l.size(); ++i) inv[i] *= l[i]; }
            std::reverse(inv.begin(), inv.end());
            res.H.torsion[n] = inv;
        }
    }
    if (ring_prime) {
        const SylowPrime* S = nullptr; for (auto& x : Dt.primes) if (x.p == ring_prime) S = &x;
        if (!S) throw std::invalid_argument("sylow: no Sylow subgroup for the ring prime in the input");
        res.R = stable_ring(*S, maxdim, stats); res.has_ring = true;
    }
    res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return res;
}

}  // namespace pcg
}  // namespace perseus
