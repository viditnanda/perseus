// Tests for Laurent (Z[t, t^-1]) coefficients: polynomial arithmetic, Alexander modules of knots (Fox calculus),
// the integral content (bad primes), and the certified barcode over Q(t) with its bad irreducibles and bad
// primes against independent oracles (evaluation at integers t = t0 through the explicit engine over Q, and the
// rank function of the boundary matrix over F_p(t) and Q(t)).
//   g++ -std=c++17 -O2 -pthread -Iinclude tests/t_laurent.cpp -o t_laurent && ./t_laurent [trials]
#include "perseus/laurent.hpp"
#include "perseus/zlin.hpp"
#include "perseus/chain.hpp"
#include "perseus/explicit_engine.hpp"
#include "perseus/certify.hpp"
#include <cstdio>
#include <random>

using namespace perseus;
using namespace perseus::laurent;
static int fails = 0, checks = 0;
static void check(bool ok, const std::string& what) { ++checks; if (!ok) { ++fails; printf("FAIL %s\n", what.c_str()); } }
static Poly P(const std::string& s) { LPoly l = parse_lpoly(s); return pshift(l.c, l.v); }   // nonnegative powers only
static std::string show(const LHomology& H) {
    std::string s;
    for (size_t k = 0; k < H.free.size(); ++k) { s += " H" + std::to_string(k) + "=L^" + std::to_string(H.free[k]); for (auto& f : H.torsion[k]) s += "+L/(" + to_string(f) + ")"; }
    return s;
}
static std::string showb(const LBarcode& B) {
    std::string s;
    for (size_t k = 0; k < B.size(); ++k) { s += " d" + std::to_string(k) + ":"; for (auto& b : B[k]) { char t[64]; snprintf(t, sizeof t, "[%g,%g)", b.first, b.second); s += t; } }
    return s;
}

// the chain complex C_2 -> C_1 -> C_0 of the infinite cyclic cover of a presentation complex: generators
// x_1..x_n all mapping to t, relators as words (+i = x_i, -i = x_i^-1); d_2 by Fox calculus, d_1 x_i = (t - 1) *
static LaurentComplex fox(int n, const std::vector<std::vector<int>>& rel) {
    LaurentComplex X;
    int32_t v = X.add_cell(0, 0, {});
    std::vector<int32_t> e;
    for (int i = 0; i < n; ++i) e.push_back(X.add_cell(1, 0, {{v, ladd(lmonomial(1, 1), lmonomial(-1, 0))}}));
    for (auto& r : rel) {
        std::vector<LPoly> d(n);
        int pref = 0;    // abelianised prefix: t^pref
        for (int x : r) {
            int i = std::abs(x) - 1;
            if (x > 0) { d[i] = ladd(d[i], lmonomial(1, pref)); pref += 1; }
            else { pref -= 1; d[i] = ladd(d[i], lmonomial(-1, pref)); }
        }
        std::vector<std::pair<int32_t, LPoly>> b;
        for (int i = 0; i < n; ++i) if (!d[i].zero()) b.push_back({e[i], d[i]});
        X.add_cell(2, 0, b);
    }
    return X;
}

// evaluate at an integer t0 != 0: columns are scaled by t0^m to clear negative powers (rank-preserving)
static FilteredComplex evaluate(const LaurentComplex& X, int64_t t0) {
    FilteredComplex C;
    for (size_t c = 0; c < X.size(); ++c) {
        int m = 0; for (auto& f : X.bd[c]) m = std::min(m, f.second.v);
        std::vector<std::pair<int32_t, int64_t>> b;
        for (auto& f : X.bd[c]) {
            Poly p = pshift(f.second.c, f.second.v - m);
            BigInt val, pw(1); for (auto& x : p) { val = val + x * pw; pw = pw * BigInt(t0); }
            if (val.bits() > 62) throw std::overflow_error("evaluate: too large");
            int64_t v = (int64_t)val.to_u64(); if (val.neg) v = -v;
            if (v) b.push_back({f.first, v});
        }
        C.add_cell(X.dim[c], X.birth[c], b);
    }
    return C;
}
static LBarcode to_lbar(const Barcode& b, int maxdim) {
    LBarcode r(maxdim + 1);
    for (size_t k = 0; k < b.size() && (int)k <= maxdim; ++k) for (auto& x : b[k]) if (x.death > x.birth) r[k].push_back({x.birth, x.death});
    for (auto& v : r) std::sort(v.begin(), v.end());
    return r;
}
static bool peval(const Poly& g, int64_t t0) { BigInt v, pw(1); for (auto& x : g) { v = v + x * pw; pw = pw * BigInt(t0); } return v.is_zero(); }

// barcode over K(t) (K = Q or F_p) from the rank function of the boundary matrix in filtration order
static LBarcode rank_oracle(const LaurentComplex& X, int maxdim, const Dom& D) {
    size_t n = X.size();
    std::vector<int32_t> ord(n); for (size_t i = 0; i < n; ++i) ord[i] = (int32_t)i;
    std::stable_sort(ord.begin(), ord.end(), [&](int32_t a, int32_t b) { return X.birth[a] != X.birth[b] ? X.birth[a] < X.birth[b] : X.dim[a] < X.dim[b]; });
    std::vector<int32_t> pos(n); for (size_t i = 0; i < n; ++i) pos[ord[i]] = (int32_t)i;
    PMat M(n, std::vector<Poly>(n));
    for (size_t j = 0; j < n; ++j) {
        int32_t c = ord[j]; int m = 0; for (auto& f : X.bd[c]) m = std::min(m, f.second.v);
        for (auto& f : X.bd[c]) { Poly p = pshift(f.second.c, f.second.v - m); D.trim(p); M[pos[f.first]][j] = p; }
    }
    auto r = [&](size_t i, long j) -> size_t {   // rank of rows >= i, columns <= j
        if (j < 0 || i >= n) return 0;
        PMat B; for (size_t a = i; a < n; ++a) B.push_back(std::vector<Poly>(M[a].begin(), M[a].begin() + j + 1));
        return smith(B, D).rank;
    };
    std::vector<char> paired(n, 0), death(n, 0);
    LBarcode out(maxdim + 1);
    std::vector<std::vector<size_t>> R(n + 1, std::vector<size_t>(n + 1, 0));
    for (size_t i = 0; i <= n; ++i) for (long j = -1; j < (long)n; ++j) R[i][j + 1] = r(i, j);
    for (size_t j = 0; j < n; ++j) for (size_t i = 0; i < n; ++i) {
        long v = (long)R[i][j + 1] - (long)R[i + 1][j + 1] - (long)R[i][j] + (long)R[i + 1][j];
        if (v == 1) {
            paired[i] = 1; death[j] = 1;
            int k = X.dim[ord[i]];
            if (k <= maxdim && X.birth[ord[j]] > X.birth[ord[i]]) out[k].push_back({X.birth[ord[i]], X.birth[ord[j]]});
        }
    }
    for (size_t i = 0; i < n; ++i) if (!paired[i] && !death[i] && X.dim[ord[i]] <= maxdim) out[X.dim[ord[i]]].push_back({X.birth[ord[i]], std::numeric_limits<double>::infinity()});
    for (auto& v : out) std::sort(v.begin(), v.end());
    return out;
}

// rank over F_p of an integer matrix evaluated at t0 (dims of H(C (x) F_p, t = t0))
static size_t rank_mod_p(std::vector<std::vector<int64_t>> A, uint64_t p) {
    size_t R = A.size(), C = R ? A[0].size() : 0, rk = 0;
    for (auto& row : A) for (auto& x : row) { x %= (int64_t)p; if (x < 0) x += p; }
    Dom D(p);
    for (size_t c = 0; c < C && rk < R; ++c) {
        size_t piv = R; for (size_t r = rk; r < R; ++r) if (A[r][c]) { piv = r; break; }
        if (piv == R) continue;
        std::swap(A[rk], A[piv]);
        uint64_t iv = D.inv((uint64_t)A[rk][c]);
        for (size_t r = rk + 1; r < R; ++r) if (A[r][c]) {
            uint64_t f = (unsigned __int128)A[r][c] * iv % p;
            for (size_t k = c; k < C; ++k) A[r][k] = (int64_t)(((unsigned __int128)A[r][k] + (unsigned __int128)(p - f) * (uint64_t)A[rk][k]) % p);
        }
        ++rk;
    }
    return rk;
}

// a random based complex over Z[t, t^-1] with known pieces, disguised by unimodular base changes
static LaurentComplex random_complex(std::mt19937_64& g, int top, bool filtered) {
    // pieces per degree k: free classes, elementary torsion L --f--> L, and L^2 --[2, t+1]--> L
    std::vector<std::vector<std::vector<LPoly>>> d(top + 2);   // d[k][row][col]: k-cells -> (k-1)-cells
    std::vector<int> nk(top + 2, 0);
    auto randpoly = [&](int dg) { LPoly r; for (int i = 0; i <= dg; ++i) r = ladd(r, lmonomial((int64_t)(g() % 5) - 2, i)); if (r.zero()) r = lmonomial(1, 0); return r; };
    std::vector<std::vector<std::pair<int, int>>> cells(top + 2);   // per degree: (piece id, slot)
    struct Ent { int k, r, c; LPoly v; };
    std::vector<Ent> ents;
    for (int k = 0; k <= top; ++k) {
        int fr = (int)(g() % 2); for (int i = 0; i < fr; ++i) ++nk[k];
        if (k < top) {
            int tor = (int)(g() % 3);
            for (int i = 0; i < tor; ++i) {
                LPoly f;
                switch (g() % 5) {
                    case 0: f = ladd(lmonomial(1, 1), lmonomial(-1, 0)); break;              // t - 1
                    case 1: f = ladd(ladd(lmonomial(1, 2), lmonomial(-1, 1)), lmonomial(1, 0)); break;   // t^2 - t + 1
                    case 2: f = ladd(lmonomial(1, 1), lmonomial(-2, 0)); break;              // t - 2
                    case 3: f = lmonomial(3, 0); break;                                     // 3 (Z-torsion)
                    default: f = randpoly(1 + (int)(g() % 2));
                }
                ents.push_back({k + 1, nk[k], nk[k + 1], f}); ++nk[k]; ++nk[k + 1];
            }
            if (g() % 4 == 0) {
                ents.push_back({k + 1, nk[k], nk[k + 1], lmonomial(2, 0)});
                ents.push_back({k + 1, nk[k], nk[k + 1] + 1, ladd(lmonomial(1, 1), lmonomial(1, 0))});
                ++nk[k]; nk[k + 1] += 2;
            }
        }
    }
    for (int k = 1; k <= top; ++k) d[k].assign(nk[k - 1], std::vector<LPoly>(nk[k]));
    for (auto& e : ents) d[e.k][e.r][e.c] = ladd(d[e.k][e.r][e.c], e.v);
    // base changes: in degree k, column operations col_a += m col_b on d_k and the inverse row operation
    // row_b -= m row_a on d_{k+1}
    for (int k = 0; k <= top; ++k) {
        if (nk[k] < 2) continue;
        for (int it = 0; it < 3; ++it) {
            int a = (int)(g() % nk[k]), b = (int)(g() % nk[k]); if (a == b) continue;
            LPoly m = lmonomial((g() % 2) ? 1 : -1, (int)(g() % 3) - 1);
            if (k >= 1) for (int r = 0; r < nk[k - 1]; ++r) d[k][r][a] = ladd(d[k][r][a], lmul(m, d[k][r][b]));
            if (k + 1 <= top) for (int c = 0; c < nk[k + 1]; ++c) d[k + 1][b][c] = ladd(d[k + 1][b][c], lmul(lmonomial(-1, 0), lmul(m, d[k + 1][a][c])));
        }
    }
    LaurentComplex X; std::vector<std::vector<int32_t>> id(top + 1);
    for (int k = 0; k <= top; ++k) {
        for (int c = 0; c < nk[k]; ++c) {
            std::vector<std::pair<int32_t, LPoly>> b; double bb = 0;
            if (k >= 1) for (int r = 0; r < nk[k - 1]; ++r) if (!d[k][r][c].zero()) { b.push_back({id[k - 1][r], d[k][r][c]}); bb = std::max(bb, X.birth[id[k - 1][r]]); }
            if (filtered) bb += (double)(g() % 3);
            id[k].push_back(X.add_cell(k, bb, b));
        }
    }
    return X;
}

int main(int argc, char** argv) {
    int trials = argc > 1 ? atoi(argv[1]) : 60;
    Dom Z;
    // ---- arithmetic ----
    {
        LPoly a = parse_lpoly("-t^-1 + 2t^3 - 5"), b = parse_lpoly("0t^2 + 3");
        check(a.str() == "2t^3 - 5 - t^-1", "parse/print " + a.str());
        check(b.str() == "3", "parse 0 coefficient " + b.str());
        Poly x = pmul(P("t-1"), P("t+2"), Z), y = pmul(P("t-1"), P("t-3"), Z);
        check(to_string(pgcd(x, y, Z)) == "t - 1", "gcd " + to_string(pgcd(x, y, Z)));
        check(to_string(pgcd(pscale(BigInt(6), x, Z), pscale(BigInt(4), y, Z), Z)) == "2t - 2", "gcd with contents");
        Poly s = pmul(pmul(P("t-1"), P("t-1"), Z), P("2t+1"), Z);
        check(to_string(sqfree(s)) == "2t^2 - t - 1", "squarefree part " + to_string(sqfree(s)));
        auto B = coprime_basis({pmul(P("t-1"), P("t+1"), Z), pmul(P("t-1"), P("t+2"), Z)});
        check(B.size() == 3, "gcd-free basis of (t-1)(t+1), (t-1)(t+2) has 3 elements");
        Poly f = P("3t^2+t+5"), q = P("2t-7");
        Poly r4 = prem_fixed(f, q, 4, Z);   // 2^4 f mod (2t - 7) = 16 f(7/2)
        check(r4.size() == 1 && r4[0] == BigInt(16 * 3 * 49 / 4 + 16 * 7 / 2 + 80), "prem_fixed value " + to_string(r4));
        check(proved_irreducible(P("t^2-t+1")) && !proved_irreducible(P("t^2-1")) && proved_irreducible(P("2t-7")), "irreducibility of quadratics");
    }
    // ---- knots: Alexander polynomials over Q[t, t^-1] ----
    {
        struct K { const char* name; int n; std::vector<std::vector<int>> rel; const char* delta; };
        std::vector<K> knots = {
            {"trefoil", 2, {{1, 2, 1, -2, -1, -2}}, "t^2 - t + 1"},                                       // xyx = yxy
            {"figure-eight", 2, {{-2, 1, 2, -1, 2, -1, -2, 1, 2, -1, -2, -1, 2, 1, -2, 1}}, ""},          // filled below
            {"T(2,5)", 2, {{1, 2, 1, 2, 1, -2, -1, -2, -1, -2}}, "t^4 - t^3 + t^2 - t + 1"},              // (xy)^2 x = y (xy)^2
            {"unknot", 1, {}, ""},
        };
        // figure-eight: <x, y | y x^-1 y x y^-1 = x^-1 y x y^-1 x> written as w x^-1 = 1, w = y x^-1 y x y^-1 ... use the
        // standard relator x y^-1 x^-1 y x^-1 y x y^-1 x^-1 ... replaced by the Wirtinger form a^-1 b a b^-1 a b a^-1 b^-1 ... keep the
        // classical: [x^-1, y] x [y^-1, x]-type relator r = y^-1 x y x^-1 y x y^-1 x^-1 ... we use w = x y^-1 x^-1 y, r = w x w^-1 y^-1
        knots[1].rel = {{1, -2, -1, 2, 1, -2, 1, 2, -1, -2}};
        knots[1].delta = "t^2 - 3t + 1";
        for (auto& k : knots) {
            LaurentComplex X = fox(k.n, k.rel);
            check(X.validate().empty(), std::string(k.name) + ": d d = 0 (" + X.validate() + ")");
            LaurentResult L = laurent_persistence(X, 2);
            const LHomology& H = L.homology;
            bool ok = H.torsion[0].size() == 1 && to_string(H.torsion[0][0]) == "t - 1" && H.free[0] == 0 && H.free[1] == 0 && H.free[2] == 0;
            if (std::string(k.delta).empty()) ok = ok && H.torsion[1].empty();
            else ok = ok && H.torsion[1].size() == 1 && to_string(H.torsion[1][0]) == k.delta;
            check(ok, std::string(k.name) + ": Alexander module" + show(H));
            check(L.homology_bad_primes.empty() && L.bad_primes.empty(), std::string(k.name) + ": no bad primes");
        }
    }
    // ---- integral content: L/(2, t + 1) and Z-torsion ----
    {
        LaurentComplex X; int32_t v = X.add_cell(0, 0, {});
        X.add_cell(1, 0, {{v, lmonomial(2, 0)}}); X.add_cell(1, 0, {{v, ladd(lmonomial(1, 1), lmonomial(1, 0))}});
        LaurentResult L = laurent_persistence(X, 1);
        check(L.homology.free == std::vector<size_t>({0, 1}) && L.homology.torsion[0].empty(), "L/(2, t+1): zero over Q[t, t^-1]" + show(L.homology));
        check(L.homology_bad_primes.size() == 1 && L.homology_bad_primes.count(2) && L.homology_bad_primes[2].torsion[0].size() == 1 && to_string(L.homology_bad_primes[2].torsion[0][0]) == "t + 1",
              "L/(2, t+1): bad prime 2 with F_2[t, t^-1]/(t + 1)");
        check(L.bad_primes.empty(), "L/(2, t+1): the barcode over F_2(t) does not see it");
    }
    // ---- a filtered complex with bad irreducibles t - 2, t^2 - t + 1 and the bad prime 2 ----
    {
        for (int twice = 0; twice < 2; ++twice) {
            LaurentComplex X; int32_t e = X.add_cell(1, 0, {});
            LPoly f = lmul(parse_lpoly("t-2"), parse_lpoly("t^2-t+1")); if (twice) f = lmul(f, lmonomial(2, -3));
            X.add_cell(2, 1, {{e, f}});
            LaurentResult L = laurent_persistence(X, 2);
            check(showb(L.bars) == " d0: d1:[0,1) d2:", "Q(t) barcode" + showb(L.bars));
            std::string bi; for (auto& b : L.bad_irreducibles) bi += "(" + to_string(b.g) + (b.irreducible ? ",irr" : "") + ")" + showb(b.bars);
            check(bi == "(t - 2,irr) d0: d1:[0,inf) d2:[1,inf)(t^2 - t + 1,irr) d0: d1:[0,inf) d2:[1,inf)", "bad irreducibles " + bi);
            check(L.bad_primes.size() == (size_t)twice && (!twice || L.bad_primes.count(2)), "bad prime 2 iff the content is 2");
            check(L.homology.torsion[1].size() == 1 && to_string(L.homology.torsion[1][0]) == "t^3 - 3t^2 + 3t - 2", "H_1 over Q[t^+-1]" + show(L.homology));
            check(L.homology_bad_primes.size() == (size_t)twice, "homology bad prime 2 iff the content is 2");
        }
        // a composite candidate that must split: pivot values (t - 1)(t + 1) in one place, t - 1 elsewhere
        LaurentComplex X; int32_t a = X.add_cell(1, 0, {}), b = X.add_cell(1, 0, {});
        X.add_cell(2, 1, {{a, parse_lpoly("t^2-1")}}); X.add_cell(2, 2, {{b, parse_lpoly("t-1")}, {a, parse_lpoly("t+1")}});
        LaurentResult L = laurent_persistence(X, 2);
        std::string bi; for (auto& x : L.bad_irreducibles) bi += "(" + to_string(x.g) + ")";
        check(bi == "(t + 1)(t - 1)", "bad irreducibles t - 1, t + 1: " + bi);
    }
    // ---- the infinite cyclic cover of the torus along one circle (a cylinder), triangulated: heavy reduction ----
    {
        std::mt19937_64 g(5);
        for (int trial = 0; trial < 6; ++trial) {
            int N = 3 + trial % 3;                           // N x N grid torus
            auto vid = [&](int i, int j) { return ((i % N + N) % N) * N + ((j % N + N) % N); };
            std::map<std::vector<int>, std::vector<int>> lifts;  // simplex (sorted ids) -> lifted x-coordinates
            for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j)
                for (auto tri : {std::vector<std::pair<int, int>>{{i, j}, {i + 1, j}, {i + 1, j + 1}}, std::vector<std::pair<int, int>>{{i, j}, {i, j + 1}, {i + 1, j + 1}}}) {
                    std::vector<std::pair<int, int>> vs;               // (id, lifted x)
                    for (auto& p : tri) vs.push_back({vid(p.first, p.second), p.first});
                    std::sort(vs.begin(), vs.end());
                    for (int m = 1; m < 8; ++m) { std::vector<int> ids, xs; for (int b = 0; b < 3; ++b) if (m >> b & 1) { ids.push_back(vs[b].first); xs.push_back(vs[b].second); } lifts[ids] = xs; }
                }
            // phi(a, b) = (x_a + (x_b - x_a) - x_b(mod N lift)) / N: the number of times the edge wraps
            auto phi = [&](int xa, int xb, int ia, int ib) { int ca = ia / N, cb = ib / N; return ((xb - xa) - (cb - ca)) / N; };
            std::vector<std::vector<int>> L; for (auto& kv : lifts) L.push_back(kv.first);
            std::stable_sort(L.begin(), L.end(), [](const std::vector<int>& a, const std::vector<int>& b) { return a.size() < b.size(); });
            std::map<std::vector<int>, int32_t> id; LaurentComplex X;
            bool filtered = trial >= 3;
            for (auto& sx : L) {
                std::vector<std::pair<int32_t, LPoly>> b; double birth = filtered ? (double)(g() % 3) : 0;
                for (size_t i = 0; i < sx.size() && sx.size() > 1; ++i) {
                    auto f = sx; f.erase(f.begin() + i);
                    int e = i == 0 ? phi(lifts[sx][0], lifts[sx][1], sx[0], sx[1]) : 0;
                    b.push_back({id[f], lmonomial((i & 1) ? -1 : 1, e)}); birth = std::max(birth, X.birth[id[f]]);
                }
                id[sx] = X.add_cell((int)sx.size() - 1, birth, b);
            }
            std::string err = X.validate();
            check(err.empty(), "twisted torus valid: " + err);
            LaurentResult Lr = laurent_persistence(X, 2), U = laurent_persistence(X, 2, true, true, false);
            check(show(Lr.homology) == " H0=L^0+L/(t - 1) H1=L^0+L/(t - 1) H2=L^0", "cover of the torus is a cylinder:" + show(Lr.homology));
            check(Lr.after_unfiltered < Lr.cells / 4, "twisted torus reduces (" + std::to_string(Lr.cells) + " -> " + std::to_string(Lr.after_unfiltered) + ")");
            bool same = U.bars == Lr.bars && same_homology(U.homology, Lr.homology) && U.bad_irreducibles.size() == Lr.bad_irreducibles.size() && U.bad_primes == Lr.bad_primes;
            for (size_t i = 0; same && i < U.bad_irreducibles.size(); ++i) same = U.bad_irreducibles[i].g == Lr.bad_irreducibles[i].g && U.bad_irreducibles[i].bars == Lr.bad_irreducibles[i].bars;
            check(same, "twisted torus: reduction preserves the barcode, bad irreducibles and homology");
            for (int64_t t0 : {-2, -1, 1, 2, 3}) {
                FilteredComplex C = evaluate(X, t0); Prepared PC(C); ExplicitBackend be(PC, 2);
                LBarcode got = to_lbar(certify(be).first.bars, 2), want = Lr.bars;
                for (auto& b : Lr.bad_irreducibles) if (peval(b.g, t0)) want = b.bars;
                check(got == want, "twisted torus: barcode at t = " + std::to_string(t0) + showb(got) + " vs" + showb(want));
            }
        }
    }
    // ---- random complexes against the oracles ----
    {
        std::mt19937_64 g(23);
        int nbadirr = 0, nbadp = 0, nhbad = 0, nred = 0;
        for (int t = 0; t < trials; ++t) {
            int top = 2 + (int)(g() % 2);
            LaurentComplex X = random_complex(g, top, t % 4 != 0);
            std::string err = X.validate(); check(err.empty(), "random complex valid: " + err);
            if (!err.empty()) continue;
            int md = top;
            LaurentResult L = laurent_persistence(X, md);
            {   // the unit-pair reduction changes nothing
                LaurentResult U = laurent_persistence(X, md, true, true, false);
                bool same = U.bars == L.bars && U.bad_primes == L.bad_primes && same_homology(U.homology, L.homology) && U.bad_irreducibles.size() == L.bad_irreducibles.size();
                for (size_t i = 0; same && i < U.bad_irreducibles.size(); ++i) same = U.bad_irreducibles[i].g == L.bad_irreducibles[i].g && U.bad_irreducibles[i].bars == L.bad_irreducibles[i].bars;
                same = same && U.homology_bad_primes.size() == L.homology_bad_primes.size();
                for (auto& kv : U.homology_bad_primes) same = same && L.homology_bad_primes.count(kv.first) && same_homology(kv.second, L.homology_bad_primes.at(kv.first));
                check(same, "unit-pair reduction preserves everything, trial " + std::to_string(t) + " (" + std::to_string(L.cells) + " -> " + std::to_string(L.after_filtered) + ", " + std::to_string(L.after_unfiltered) + ")");
                nred += L.after_unfiltered < L.cells;
            }
            nbadirr += !L.bad_irreducibles.empty(); nbadp += !L.bad_primes.empty(); nhbad += !L.homology_bad_primes.empty();
            // Q(t) and F_p(t) barcodes against the rank function
            if (X.size() <= 22) {
                check(rank_oracle(X, md, Dom()) == L.bars, "Q(t) barcode = rank oracle, trial " + std::to_string(t) + showb(L.bars) + " vs" + showb(rank_oracle(X, md, Dom())));
                for (uint64_t p : {2ULL, 3ULL, 5ULL}) {
                    LBarcode want = L.bad_primes.count(p) ? L.bad_primes.at(p) : L.bars;
                    check(rank_oracle(X, md, Dom(p)) == want, "F_" + std::to_string(p) + "(t) barcode = rank oracle (bad-prime list complete), trial " + std::to_string(t));
                }
            }
            // evaluation at t = t0: the barcode over Q is that of the bad irreducible vanishing at t0, if any
            for (int64_t t0 : {-3, -2, -1, 1, 2, 3, 5}) {
                FilteredComplex C = evaluate(X, t0);
                Prepared PC(C); ExplicitBackend be(PC, std::max(0, std::min(md, PC.D)));
                auto q = certify(be);
                LBarcode got = to_lbar(q.first.bars, md);
                LBarcode want = L.bars; int hits = 0;
                for (auto& b : L.bad_irreducibles) if (peval(b.g, t0)) { want = b.bars; ++hits; }
                check(hits <= 1 && got == want, "barcode at t = " + std::to_string(t0) + " vs bad irreducibles, trial " + std::to_string(t) + ":" + showb(got) + " vs" + showb(want));
                // homology over Q at t0: free ranks plus roots of the invariant factors (universal coefficients)
                IntegralHomology Hq = integral_homology(C, md);
                for (int k = 0; k <= md; ++k) {
                    size_t pred = L.homology.free[k];
                    for (auto& f : L.homology.torsion[k]) pred += peval(f, t0);
                    if (k > 0) for (auto& f : L.homology.torsion[k - 1]) pred += peval(f, t0);
                    check(Hq.free[k] == pred, "rank of H_" + std::to_string(k) + " at t = " + std::to_string(t0) + ", trial " + std::to_string(t));
                }
            }
            // F_p[t, t^-1] at t0 in F_p^*: dimensions predicted by the reduced (or the bad-prime) answer
            for (uint64_t p : {2ULL, 3ULL, 5ULL, 7ULL}) {
                LHomology Hp = L.homology_bad_primes.count(p) ? L.homology_bad_primes.at(p) : reduce_homology(L.homology, p);
                for (uint64_t t0 = 1; t0 < p; ++t0) {
                    FilteredComplex C = evaluate(X, (int64_t)t0);
                    std::vector<std::vector<int32_t>> cells(md + 3);
                    for (size_t c = 0; c < C.size(); ++c) if (C.dim[c] <= md + 1) cells[C.dim[c]].push_back((int32_t)c);
                    std::vector<size_t> rk(md + 3, 0);
                    for (int k = 1; k <= md + 1; ++k) if (!cells[k].empty() && !cells[k - 1].empty()) rk[k] = rank_mod_p(perseus::detail::boundary_block(C, cells[k - 1], cells[k]), p);
                    auto root = [&](const Poly& f) { BigInt v, pw(1); for (auto& x : f) { v = v + x * pw; pw = pw * BigInt((int64_t)t0); } return v.mod_small(p) == 0; };
                    for (int k = 0; k <= md; ++k) {
                        size_t dimk = cells[k].size() - rk[k] - rk[k + 1];
                        size_t pred = Hp.free[k];
                        for (auto& f : Hp.torsion[k]) pred += root(f);
                        if (k > 0) for (auto& f : Hp.torsion[k - 1]) pred += root(f);
                        check(dimk == pred, "F_" + std::to_string(p) + " dimension of H_" + std::to_string(k) + " at t = " + std::to_string(t0) + ", trial " + std::to_string(t) + show(Hp));
                    }
                }
            }
        }
        printf("(random complexes: %d with bad irreducibles, %d with bad primes over F_p(t), %d with homology bad primes, %d reduced)\n", nbadirr, nbadp, nhbad, nred);
    }
    printf("%d checks, %d failures\n", checks, fails);
    return fails ? 1 : 0;
}
