// t_pcgroup.cpp -- tests for pc_group.hpp: homology and mod-p cohomology rings of finite solvable groups.
//   1. closed forms: cyclic groups, elementary abelian groups, generalised quaternion and dihedral groups, S_3;
//   2. the Kunneth formula on direct products (mixed primes and mixed relative orders included);
//   3. the universal coefficient theorem: the mod-p ring engine against the integral engine, for every group
//      and every prime dividing its order;
//   4. ground truth from HAP and the Green-King database (tests/fixtures/pcgroups/expected.txt);
//   5. invalid inputs are rejected;
//   6. the Sylow route (sylow.hpp): against HAP and the literature for non-solvable groups, against the direct
//      engine for solvable ones, universal coefficients between its integral and mod-p results, invalid inputs.
#include "perseus/pc_group.hpp"
#include "perseus/sylow.hpp"
#include <cstdio>
#include <fstream>
#include <sstream>
#include <functional>
#include <array>
#include <map>
#include <set>
using namespace perseus::pcg;

static int fails = 0, checks = 0;
#define CHECK(c, msg) do { ++checks; if (!(c)) { ++fails; printf("FAIL: %s  [%s:%d]\n", std::string(msg).c_str(), __FILE__, __LINE__); } } while (0)

struct Grp { std::vector<std::vector<int>> T; std::vector<int> pcgs, rel; int e = 0; std::string name; };

static Grp make(int n, std::function<int(int, int)> mul, int e, std::vector<int> pcgs, std::vector<int> rel, std::string name) {
    Grp G; G.T.assign(n, std::vector<int>(n)); for (int a = 0; a < n; ++a) for (int b = 0; b < n; ++b) G.T[a][b] = mul(a, b);
    G.e = e; G.pcgs = pcgs; G.rel = rel; G.name = name; return G;
}
static std::vector<int> prime_factors(int n) { std::vector<int> f; for (int p = 2; p * p <= n; ++p) while (n % p == 0) { f.push_back(p); n /= p; } if (n > 1) f.push_back(n); return f; }
static Grp cyclic(int n) {                     // pcgs g, g^p1, g^(p1 p2), ...
    std::vector<int> pcgs, rel; int s = 1; for (int p : prime_factors(n)) { pcgs.push_back(s % n); rel.push_back(p); s *= p; }
    return make(n, [n](int a, int b) { return (a + b) % n; }, 0, pcgs, rel, "C" + std::to_string(n));
}
static Grp elem_abelian(int p, int k) {
    int n = 1; for (int i = 0; i < k; ++i) n *= p;
    std::vector<int> pcgs, rel; int s = n / p; for (int i = 0; i < k; ++i) { pcgs.push_back(s); rel.push_back(p); s /= p; }
    return make(n, [p, k](int a, int b) { int r = 0, m = 1; for (int i = 0; i < k; ++i) { r += ((a / m % p + b / m % p) % p) * m; m *= p; } return r; }, 0, pcgs, rel, "(Z/" + std::to_string(p) + ")^" + std::to_string(k));
}
// x^a y^b, a < N, b < 2, with y x y^-1 = x^-1 and y^2 = x^(N/2) (quaternion) or 1 (dihedral); code a + N b
static Grp metacyclic(int N, bool quaternion, std::string name) {
    auto mul = [N, quaternion](int u, int v) { int a = u % N, b = u / N, c = v % N, d = v / N; int a2 = a + (b ? N - c : c); int b2 = b ^ d; if (quaternion && b && d) a2 += N / 2; return (a2 % N) + N * b2; };
    std::vector<int> pcgs{N}, rel{2}; int s = 1; for (int p : prime_factors(N)) { pcgs.push_back(s); rel.push_back(p); s *= p; }
    return make(2 * N, mul, 0, pcgs, rel, name);
}
static Grp product(const Grp& A, const Grp& B) {
    int na = (int)A.T.size(), nb = (int)B.T.size();
    std::vector<int> pcgs, rel;
    for (size_t i = 0; i < A.pcgs.size(); ++i) { pcgs.push_back(A.pcgs[i] * nb + B.e); rel.push_back(A.rel[i]); }
    for (size_t i = 0; i < B.pcgs.size(); ++i) { pcgs.push_back(A.e * nb + B.pcgs[i]); rel.push_back(B.rel[i]); }
    return make(na * nb, [&A, &B, nb](int u, int v) { return A.T[u / nb][v / nb] * nb + B.T[u % nb][v % nb]; }, A.e * nb + B.e, pcgs, rel, A.name + " x " + B.name);
}
static Grp read_fixture(const std::string& path) {
    std::ifstream in(path); std::string line, all;
    while (std::getline(in, line)) { size_t h = line.find('#'); if (h != std::string::npos) line.resize(h); all += line + " "; }
    std::istringstream is(all); int n, k; is >> n >> k;
    Grp G; G.pcgs.resize(k); G.rel.resize(k);
    for (auto& x : G.pcgs) is >> x;
    for (auto& x : G.rel) is >> x;
    is >> G.e;
    G.T.assign(n, std::vector<int>(n)); for (auto& r : G.T) for (auto& x : r) is >> x;
    G.name = path; return G;
}

// cyclic summands as a sorted list of prime powers
static std::vector<u64> prime_powers(const std::vector<u64>& inv) {
    std::vector<u64> out;
    for (u64 x : inv) for (u64 p = 2; x > 1; ++p) if (x % p == 0) { u64 q = 1; while (x % p == 0) { x /= p; q *= p; } out.push_back(q); }
    std::sort(out.begin(), out.end()); return out;
}
static u64 gcd(u64 a, u64 b) { while (b) { u64 t = a % b; a = b; b = t; } return a; }

static GroupHomology homology(const Grp& G, int D) { PcGroup P(G.T, G.pcgs, G.rel, G.e); return group_homology(P, D); }

static void check_uct(const Grp& G, const GroupHomology& H, int D) {
    for (int p : prime_factors((int)G.T.size())) {
        static std::vector<int> done; (void)done;
        PcGroup P(G.T, G.pcgs, G.rel, G.e);
        CohomologyRing R = cohomology_ring(P, D, p);
        for (int n = 0; n <= D; ++n) {
            auto t = [&](int m) { size_t c = 0; if (m < 0) return c; for (u64 x : H.torsion[m]) if (x % (u64)p == 0) ++c; return c; };
            size_t expect = H.free[n] + t(n) + t(n - 1);
            CHECK(R.dims[n] == expect, G.name + ": UCT mod " + std::to_string(p) + " in degree " + std::to_string(n) + ": ring " + std::to_string(R.dims[n]) + " vs " + std::to_string(expect));
        }
    }
}


// The ring products against the Alexander-Whitney cup product, on cocycles: pull the basis cocycles back to the
// normalised bar resolution through the comparison map f_m = h o f_{m-1} o d, and compare the AW cup of the pull-backs
// with the pull-back of our product modulo bar coboundaries; also against (-1)^(ab) times it. Returns
// {agree only with the cup, agree only with the twisted product, agree with both (sign invisible), neither}.
static std::array<int, 4> cup_versus_aw(const Grp& g, int D) {
    PcGroup P(g.T, g.pcgs, g.rel, g.e);
    int p = P.rel[0], n = P.n; const u64 N = (u64)n;
    CohomologyRing R = cohomology_ring(P, D, p);
    Anick A(P, D + 1, p);
    std::vector<int> ne; for (int x = 0; x < n; ++x) if (x != P.e) ne.push_back(x);
    std::map<std::vector<int>, Vec> fmemo;
    std::function<Vec(const std::vector<int>&)> f = [&](const std::vector<int>& b) -> Vec {
        auto it = fmemo.find(b); if (it != fmemo.end()) return it->second;
        size_t m = b.size(); Vec out;
        if (m == 0) out = Vec{{(u64)P.e, 1}};
        else {
            std::vector<std::pair<u64, i64>> acc;
            auto add = [&](const Vec& v, i64 s, int right) { for (auto& kv : v) acc.push_back({(kv.first / N) * N + (u64)P.mul((int)(kv.first % N), right), ((s * kv.second) % p + p) % p}); };
            add(f(std::vector<int>(b.begin(), b.end() - 1)), 1, b[m - 1]);         // d[g1..gm] = [g1..g_{m-1}] g_m + ...
            for (size_t i = 0; i + 1 < m; ++i) { int y = P.mul(b[i], b[i + 1]); if (y == P.e) continue; std::vector<int> c(b); c[i] = y; c.erase(c.begin() + i + 1); add(f(c), ((m - (i + 1)) % 2) ? -1 : 1, P.e); }
            add(f(std::vector<int>(b.begin() + 1, b.end())), (m % 2) ? -1 : 1, P.e);
            Vec X; normalize(acc, p, X); A.hvec((int)m - 1, X, out);
        }
        return fmemo[b] = out;
    };
    auto pull = [&](int deg, const std::vector<std::pair<int32_t, uint32_t>>& coc, const std::vector<int>& b) {
        std::vector<i64> c(A.C.size(deg), 0); for (auto& kv : coc) c[kv.first] = kv.second;
        i64 s = 0; for (auto& kv : f(b)) s += kv.second * c[kv.first / N]; return s % p;
    };
    auto cells = [&](int m) { std::vector<std::vector<int>> out(1); for (int k = 0; k < m; ++k) { std::vector<std::vector<int>> o2; for (auto& c : out) for (int x : ne) { auto d = c; d.push_back(x); o2.push_back(d); } out.swap(o2); } return out; };
    std::array<int, 4> res{0, 0, 0, 0};
    for (int b = 1; 1 + b <= D; ++b) {
        int a = 1, m = a + b; auto cm = cells(m);
        std::map<std::vector<int>, std::map<int32_t, i64>> rows;               // coboundaries of indicators of (m-1)-cells
        for (size_t j = 0; j < cm.size(); ++j) {
            const auto& c0 = cm[j];
            auto term = [&](const std::vector<int>& face, int s) { i64& v = rows[face][(int32_t)j]; v = ((v + s) % p + p) % p; };
            term(std::vector<int>(c0.begin(), c0.end() - 1), 1);
            for (size_t i = 0; i + 1 < (size_t)m; ++i) { int y = P.mul(c0[i], c0[i + 1]); if (y == P.e) continue; std::vector<int> c(c0); c[i] = y; c.erase(c.begin() + i + 1); term(c, ((m - (i + 1)) % 2) ? -1 : 1); }
            term(std::vector<int>(c0.begin() + 1, c0.end()), (m % 2) ? -1 : 1);
        }
        Echelon E(p); for (auto& kv : rows) { SV v; for (auto& t : kv.second) if (t.second) v.push_back({t.first, (uint32_t)t.second}); E.insert(v); }
        for (size_t i = 0; i < R.basis[a].size(); ++i) for (size_t j = 0; j < R.basis[b].size(); ++j) {
            auto it = R.products.find(std::make_tuple(a, (int)i, b, (int)j)); if (it == R.products.end()) continue;
            std::vector<i64> ours(A.C.size(m), 0);
            for (size_t t = 0; t < it->second.size(); ++t) if (it->second[t]) for (auto& kv : R.basis[m][t]) ours[kv.first] = (ours[kv.first] + (i64)it->second[t] * kv.second) % p;
            std::vector<std::pair<int32_t, uint32_t>> oc; for (size_t c = 0; c < ours.size(); ++c) if (ours[c]) oc.push_back({(int32_t)c, (uint32_t)ours[c]});
            SV d1, d2;
            for (size_t c = 0; c < cm.size(); ++c) {
                const auto& bc = cm[c];
                i64 cup = pull(a, R.basis[a][i], std::vector<int>(bc.begin(), bc.begin() + a)) * pull(b, R.basis[b][j], std::vector<int>(bc.begin() + a, bc.end())) % p;
                i64 pr = pull(m, oc, bc), tw = (a * b) % 2 ? (p - pr) % p : pr;
                i64 x = ((cup - pr) % p + p) % p, y = ((cup - tw) % p + p) % p;
                if (x) d1.push_back({(int32_t)c, (uint32_t)x});
                if (y) d2.push_back({(int32_t)c, (uint32_t)y});
            }
            E.reduce(d1); E.reduce(d2);
            if (d1.empty() && d2.empty()) ++res[2]; else if (d1.empty()) ++res[0]; else if (d2.empty()) ++res[1]; else ++res[3];
        }
    }
    return res;
}

int main() {
    // 1. closed forms
    for (int n : {2, 3, 4, 5, 6, 8, 9, 12, 30, 32}) {
        Grp G = cyclic(n); int D = 7; GroupHomology H = homology(G, D);
        CHECK(H.free[0] == 1, G.name + " H_0");
        for (int d = 1; d <= D; ++d) {
            bool ok = H.free[d] == 0 && (d % 2 ? (H.torsion[d] == std::vector<u64>{(u64)n}) : H.torsion[d].empty());
            CHECK(ok, G.name + ": H_" + std::to_string(d));
        }
        check_uct(G, H, 6);
        for (int p : prime_factors(n)) {
            PcGroup P(G.T, G.pcgs, G.rel, G.e); CohomologyRing R = cohomology_ring(P, 6, p);
            bool poly = (p == 2 && n % 4 != 0);    // F_2[x] when 4 does not divide n, else Lambda(x) (x) F_p[y]
            CHECK(R.indecomposables[1] == 1 && R.indecomposables[2] == (poly ? 0u : 1u) && R.indecomposables[3] == 0, G.name + ": ring generators mod " + std::to_string(p));
        }
    }
    for (auto pk : std::vector<std::pair<int, int>>{{2, 1}, {2, 3}, {2, 4}, {3, 2}, {3, 3}, {5, 2}}) {
        Grp G = elem_abelian(pk.first, pk.second); int D = 6, p = pk.first, k = pk.second;
        PcGroup P(G.T, G.pcgs, G.rel, G.e); CohomologyRing R = cohomology_ring(P, D, p);
        for (int n = 0; n <= D; ++n) { u64 c = 1; for (int i = 1; i <= k - 1; ++i) c = c * (n + i) / i; CHECK(R.dims[n] == c, G.name + ": dim H^" + std::to_string(n)); }
        CHECK(R.indecomposables[1] == (size_t)k && R.indecomposables[2] == (p == 2 ? 0u : (size_t)k) && R.indecomposables[3] == 0 && R.indecomposables[4] == 0, G.name + ": ring generators");
        check_uct(G, homology(G, D), D);
    }
    for (int m = 3; m <= 6; ++m) {                // generalised quaternion Q_{2^m}: period 4
        int N = 1 << (m - 1); Grp G = metacyclic(N, true, "Q" + std::to_string(2 * N)); int D = 8; GroupHomology H = homology(G, D);
        for (int d = 1; d <= D; ++d) {
            std::vector<u64> want = d % 4 == 1 ? std::vector<u64>{2, 2} : d % 4 == 3 ? std::vector<u64>{(u64)(2 * N)} : std::vector<u64>{};
            CHECK(H.free[d] == 0 && H.torsion[d] == want, G.name + ": H_" + std::to_string(d));
        }
        check_uct(G, H, D);
    }
    for (int m = 3; m <= 5; ++m) {                // dihedral D_{2N}, N = 2^(m-1) (Handel)
        int N = 1 << (m - 1); Grp G = metacyclic(N, false, "D" + std::to_string(2 * N)); int D = 8; GroupHomology H = homology(G, D);
        for (int d = 1; d <= D; ++d) {
            int k = (d - 1) / 4, r = (d - 1) % 4; size_t twos = r == 0 ? 2 * k + 2 : r == 1 ? 2 * k + 1 : 2 * k + 2;
            std::vector<u64> want(twos, 2); if (r == 2) want.push_back(N); if (N == 2 && r == 2) { want.pop_back(); }
            std::vector<u64> pp = prime_powers(H.torsion[d]), wp = prime_powers(want);
            CHECK(H.free[d] == 0 && pp == wp, G.name + ": H_" + std::to_string(d));
        }
        check_uct(G, H, D);
    }
    {                                             // S_3 = D_6: mixed relative orders (2, 3); period 4
        Grp G = metacyclic(3, false, "S3"); int D = 8; GroupHomology H = homology(G, D);
        for (int d = 1; d <= D; ++d) {
            std::vector<u64> want = d % 4 == 1 ? std::vector<u64>{2} : d % 4 == 3 ? std::vector<u64>{6} : std::vector<u64>{};
            CHECK(H.free[d] == 0 && H.torsion[d] == want, "S3: H_" + std::to_string(d));
        }
        check_uct(G, H, 6);
    }
    // 2. Kunneth on direct products
    std::vector<std::pair<Grp, Grp>> pairs = { {cyclic(2), cyclic(4)}, {metacyclic(4, true, "Q8"), cyclic(3)}, {metacyclic(3, false, "S3"), cyclic(2)},
                                               {metacyclic(4, false, "D8"), cyclic(2)}, {cyclic(6), cyclic(4)}, {metacyclic(3, false, "S3"), metacyclic(3, false, "S3")} };
    for (auto& ab : pairs) {
        int D = 6; Grp P = product(ab.first, ab.second);
        GroupHomology A = homology(ab.first, D), B = homology(ab.second, D), C = homology(P, D);
        for (int n = 1; n <= D; ++n) {
            std::vector<u64> want;
            auto summands = [](const GroupHomology& H, int i) { std::vector<u64> s(H.free[i], 0); for (u64 x : H.torsion[i]) s.push_back(x); return s; };   // 0 = Z
            for (int i = 0; i <= n; ++i) for (u64 x : summands(A, i)) for (u64 y : summands(B, n - i)) {
                u64 t = x == 0 ? y : y == 0 ? x : gcd(x, y); if (t != 1) want.push_back(t); }
            for (int i = 0; i <= n - 1; ++i) for (u64 x : summands(A, i)) for (u64 y : summands(B, n - 1 - i)) {
                if (x && y) { u64 t = gcd(x, y); if (t != 1) want.push_back(t); } }
            size_t wfree = std::count(want.begin(), want.end(), (u64)0);
            std::vector<u64> wt; for (u64 x : want) if (x) wt.push_back(x);
            CHECK(C.free[n] == wfree && prime_powers(C.torsion[n]) == prime_powers(wt), P.name + ": Kunneth in degree " + std::to_string(n));
        }
        check_uct(P, C, 5);
    }
    // 3b. the products are cup products: against Alexander-Whitney on the bar resolution, on cocycles
    for (auto ng : std::vector<std::pair<std::string, Grp>>{{"(Z/3)^2", elem_abelian(3, 2)}, {"(Z/5)^2", elem_abelian(5, 2)}, {"(Z/2)^3", elem_abelian(2, 3)},
                                                           {"Q8", metacyclic(4, true, "Q8")}, {"SmallGroup(27,3)", read_fixture("tests/fixtures/pcgroups/g27_3.txt")}}) {
        auto r = cup_versus_aw(ng.second, 3);
        CHECK(r[3] == 0, ng.first + ": a product agrees with neither the cup product nor its twist");
        CHECK(r[1] == 0, ng.first + ": a product agrees only with (-1)^(ab) times the cup product");
        if (ng.first == "(Z/3)^2" || ng.first == "(Z/5)^2") CHECK(r[0] > 0, ng.first + ": no pair where the sign is visible was tested");
    }
    // 3c. the F_p splitting behind the cohomology bases, and the two product modes
    for (auto ng : std::vector<std::tuple<std::string, Grp, int, int>>{{"Q8", metacyclic(4, true, "Q8"), 2, 6}, {"SmallGroup(27,3)", read_fixture("tests/fixtures/pcgroups/g27_3.txt"), 3, 6},
                                                                       {"SmallGroup(64,138)", read_fixture("tests/fixtures/pcgroups/g64_138.txt"), 2, 6}, {"S4 mod 2", read_fixture("tests/fixtures/pcgroups/g24_12.txt"), 2, 6},
                                                                       {"S4 mod 3", read_fixture("tests/fixtures/pcgroups/g24_12.txt"), 3, 6}, {"(Z/3)^3", elem_abelian(3, 3), 3, 5}}) {
        const std::string& nm = std::get<0>(ng); const Grp& g = std::get<1>(ng); int p = std::get<2>(ng), D = std::get<3>(ng);
        PcGroup P(g.T, g.pcgs, g.rel, g.e);
        Anick A(P, D + 1, p); FpSplit F = fp_split(A, D, p);
        bool cyc = true, coc = true, dual = true;
        for (int n = 0; n <= D; ++n) {
            auto Ep = A.eps(n + 1);                                    // eps(d_{n+1}) as columns over chains_n
            auto Ed = n ? A.eps(n) : decltype(Ep)();
            auto cocs = fp_cocycles(F, n, A.C.size(n));
            std::vector<SV> reps; for (int32_t c : F.S[n]) reps.push_back(n == 0 ? SV{{c, 1u}} : F.E[n].trans[c]);
            for (auto& r : reps) if (n) {                               // eps(d_n) f0(s) = 0
                std::map<int32_t, i64> acc; for (auto& kv : r) for (auto& e : Ed[kv.first]) acc[e.first] = (acc[e.first] + e.second * kv.second) % p;
                for (auto& kv : acc) if (kv.second) cyc = false;
            }
            for (size_t i = 0; i < cocs.size(); ++i) {
                std::vector<i64> phi(A.C.size(n), 0); for (auto& kv : cocs[i]) phi[kv.first] = kv.second;
                for (auto& col : Ep) { i64 t = 0; for (auto& e : col) t += phi[e.first] * e.second; if (t % p) coc = false; }   // phi o eps(d_{n+1}) = 0
                for (size_t j = 0; j < reps.size(); ++j) { i64 t = 0; for (auto& kv : reps[j]) t += phi[kv.first] * kv.second; if ((t % p) != (i == j ? 1 : 0)) dual = false; }
            }
        }
        CHECK(cyc, nm + ": the representatives f0(s) are cycles");
        CHECK(coc, nm + ": the basis cochains g0^*(s) are cocycles");
        CHECK(dual, nm + ": the cocycles and the representatives are dual bases");
        CohomologyRing Rg = cohomology_ring(P, D, p), Ra = cohomology_ring(P, D, p, nullptr, true);
        CHECK(Rg.dims == Ra.dims && Rg.indecomposables == Ra.indecomposables && Rg.generators == Ra.generators, nm + ": generators-only and all-products modes find the same generators");
        std::set<std::pair<int, int>> gens(Rg.generators.begin(), Rg.generators.end());
        bool keys = true, same = true; size_t common = 0;
        for (auto& kv : Rg.products) {
            if (!gens.count({std::get<0>(kv.first), std::get<1>(kv.first)})) keys = false;
            auto it = Ra.products.find(kv.first); if (it != Ra.products.end()) { ++common; if (it->second != kv.second) same = false; }
        }
        for (auto& gn : gens) for (int b = 1; gn.first + b <= D; ++b) for (size_t j = 0; j < Rg.dims[b]; ++j) if (!Rg.products.count(std::make_tuple(gn.first, gn.second, b, (int)j))) keys = false;
        CHECK(keys, nm + ": generators-only mode lists exactly the products of each generator with every class");
        CHECK(same && common > 0, nm + ": the two modes agree on their common products");
    }
    // 4. HAP and Green-King ground truth
    std::string dir = "tests/fixtures/pcgroups/";
    std::ifstream ex(dir + "expected.txt"); std::string line;
    while (std::getline(ex, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ls(line); std::string g, kind; ls >> g >> kind;
        Grp G = read_fixture(dir + g + ".txt");
        if (kind == "H") {
            std::string spec; ls >> spec; std::vector<std::vector<u64>> want; std::string field; std::istringstream fs(spec + ";");
            while (std::getline(fs, field, ';')) { std::vector<u64> v; std::istringstream xs(field); std::string t; while (std::getline(xs, t, ',')) if (!t.empty()) v.push_back(std::stoull(t)); want.push_back(v); }
            want.pop_back();                      // the sentinel field added above
            int D = (int)want.size(); GroupHomology H = homology(G, D);
            for (int d = 1; d <= D; ++d) CHECK(H.free[d] == 0 && H.torsion[d] == want[d - 1], g + ": H_" + std::to_string(d) + " against HAP");
            check_uct(G, H, std::min(D, 5));
        } else if (kind == "R") {
            int p; std::string w; std::vector<size_t> dims, gens; ls >> p;
            while (ls >> w) { std::string list; ls >> list; std::istringstream xs(list); std::string t; while (std::getline(xs, t, ',')) (w == "dims" ? dims : gens).push_back(std::stoul(t)); }
            int D = 10;                           // the whole Green-King range in the fixtures
            PcGroup P(G.T, G.pcgs, G.rel, G.e); CohomologyRing R = cohomology_ring(P, D, p);
            for (int n = 0; n <= D; ++n) CHECK(R.dims[n] == dims[n], g + ": dim H^" + std::to_string(n) + "(F_" + std::to_string(p) + ") against Green-King");
            for (int n = 1; n <= D; ++n) CHECK(R.indecomposables[n] == gens[n - 1], g + ": generators in degree " + std::to_string(n) + " against Green-King");
        }
    }
    // 5. invalid inputs
    auto rejects = [](std::function<void()> f) { try { f(); } catch (const std::invalid_argument&) { return true; } return false; };
    { Grp G = cyclic(4); CHECK(rejects([&] { PcGroup(G.T, {2, 1}, {2, 2}, 0); }), "C4 with the pcgs in the wrong order is rejected"); }
    { Grp G = cyclic(6); CHECK(rejects([&] { PcGroup(G.T, G.pcgs, {2, 2}, 0); }), "relative orders that do not multiply to |G| are rejected"); }
    { Grp G = metacyclic(4, false, "D8"); std::swap(G.T[1][2], G.T[1][3]); CHECK(rejects([&] { PcGroup(G.T, G.pcgs, G.rel, G.e); }), "a table that is not a group is rejected"); }
    { Grp G = cyclic(8); auto T = G.T; T[3][5] = 1; CHECK(rejects([&] { PcGroup(T, G.pcgs, G.rel, 0); }), "a corrupted table is rejected"); }
    { Grp G = metacyclic(4, false, "D8"); CHECK(rejects([&] { PcGroup(G.T, {1, 2, 4}, {4, 2, 2}, 0); }), "a sequence that is not a pcgs is rejected"); }
    {   // a loop that is not a group: (Z/2)^3 with an intercalate swapped. It has an identity, is a Latin square and passes
        // every pcgs check, so only the associativity test (Light's, on the pcgs) can reject it
        Grp G = elem_abelian(2, 3); auto& T = G.T;
        std::swap(T[1][4], T[1][7]); std::swap(T[2][4], T[2][7]);
        bool latin = true; for (int a = 0; a < 8; ++a) { std::vector<int> r(8, 0), c(8, 0); for (int b = 0; b < 8; ++b) { if (r[T[a][b]]++) latin = false; if (c[T[b][a]]++) latin = false; } }
        CHECK(latin, "the intercalate swap keeps a Latin square");
        std::string msg; try { PcGroup(T, G.pcgs, G.rel, 0); } catch (const std::invalid_argument& ex) { msg = ex.what(); }
        CHECK(msg.find("associative") != std::string::npos, "a non-associative loop is rejected by the associativity test: " + msg);
    }
    // 6. the Sylow route
    {
        auto load = [](const std::string& path) { std::ifstream in(path); std::stringstream ss; ss << in.rdbuf(); return parse_sylow(ss.str()); };
        std::string sd = "tests/fixtures/sylow/";
        auto pp = [](const std::vector<u64>& inv) { return prime_powers(inv); };
        std::ifstream ex2(sd + "expected.txt"); std::string line2;
        while (std::getline(ex2, line2)) {
            if (line2.empty() || line2[0] == '#') continue;
            std::istringstream ls(line2); std::string g, kind; ls >> g >> kind;
            SylowData Dt = load(sd + g + ".txt");
            if (kind == "H") {
                std::string spec; ls >> spec; std::vector<std::vector<u64>> want; std::string field; std::istringstream fs(spec);
                while (std::getline(fs, field, ';')) { std::vector<u64> v; std::istringstream xs(field); std::string t; while (std::getline(xs, t, ',')) if (t != "0") v.push_back(std::stoull(t)); std::sort(v.begin(), v.end()); want.push_back(v); }
                int D = (int)want.size(); SylowResult R = sylow_compute(Dt, D, true, 0);
                CHECK(R.complete, g + ": every prime is present");
                for (int d = 1; d <= D; ++d) CHECK(R.H.free[d] == 0 && pp(R.H.torsion[d]) == want[d - 1], g + ": H_" + std::to_string(d) + " through Sylow subgroups against HAP");
            } else if (kind == "R") {
                int p; std::string w; std::vector<size_t> dims, gens; ls >> p;
                while (ls >> w) { std::string list; ls >> list; std::istringstream xs(list); std::string t; while (std::getline(xs, t, ',')) (w == "dims" ? dims : gens).push_back(std::stoul(t)); }
                SylowResult R = sylow_compute(Dt, 9, false, p);
                for (int n = 0; n <= 9; ++n) CHECK(R.R.dims[n] == dims[n], g + ": dim H^" + std::to_string(n) + "(F_" + std::to_string(p) + ") against the literature");
                for (int n = 1; n <= 9; ++n) CHECK(R.R.indecomposables[n] == gens[n - 1], g + ": mod-" + std::to_string(p) + " generators in degree " + std::to_string(n) + " against the literature");
            }
        }
        // universal coefficients between the integral and the mod-p results of the Sylow route, every prime
        for (std::string g : {"a5", "psl27", "a6", "m11", "a7", "s4", "gl23", "c4c4s3"}) {
            SylowData Dt = load(sd + g + ".txt"); int D = 6;
            SylowResult H = sylow_compute(Dt, D, true, 0);
            for (auto& S : Dt.primes) {
                SylowResult R = sylow_compute(Dt, D, false, S.p);
                auto t = [&](int m) { size_t c = 0; if (m < 1) return c; for (u64 x : H.H.torsion[m]) if (x % (u64)S.p == 0) ++c; return c; };
                for (int n = 1; n <= D; ++n) CHECK(R.R.dims[n] == t(n) + t(n - 1), g + ": UCT mod " + std::to_string(S.p) + " in degree " + std::to_string(n) + " (Sylow route)");
            }
        }
        // against the direct engine on solvable groups
        for (auto gp : std::vector<std::pair<std::string, std::string>>{{"s4", "g24_12"}, {"gl23", "g48_29"}, {"c4c4s3", "g96_64"}}) {
            SylowData Dt = load(sd + gp.first + ".txt"); Grp G = read_fixture("tests/fixtures/pcgroups/" + gp.second + ".txt");
            PcGroup P(G.T, G.pcgs, G.rel, G.e); int D = 6;
            SylowResult S = sylow_compute(Dt, D, true, 0); GroupHomology H = group_homology(P, D);
            bool same = true; for (int d = 1; d <= D; ++d) if (S.H.torsion[d] != H.torsion[d]) same = false;
            CHECK(same, gp.first + ": integral homology through Sylow subgroups agrees with the direct engine");
            for (auto& Sp : Dt.primes) {
                SylowResult R = sylow_compute(Dt, D, false, Sp.p); CohomologyRing C = cohomology_ring(P, D, Sp.p);
                CHECK(R.R.dims == C.dims && R.R.indecomposables == C.indecomposables, gp.first + ": mod-" + std::to_string(Sp.p) + " ring (dims, generators) through Sylow subgroups agrees with the direct engine");
            }
        }
        // a prime left out: only the other p-parts, flagged incomplete
        { SylowData Dt = load(sd + "a5.txt"); Dt.primes.pop_back(); SylowResult R = sylow_compute(Dt, 3, true, 0);
          CHECK(!R.complete && R.H.torsion[3] == std::vector<u64>{6}, "A5 without its Sylow 5-subgroup: H_3 is the 2- and 3-parts, Z/6, flagged incomplete"); }
        // invalid inputs
        auto rejects = [](std::function<void()> f) { try { f(); } catch (const std::invalid_argument&) { return true; } return false; };
        { SylowData Dt = load(sd + "a5.txt"); auto& C = Dt.primes[0].cosets[0]; std::swap(C.images[0], C.images[1]);
          CHECK(rejects([&] { sylow_compute(Dt, 3, true, 0); }), "a conjugation map that is not a homomorphism is rejected"); }
        { SylowData Dt = load(sd + "a5.txt"); Dt.primes[0].cosets[0].elements.pop_back();
          CHECK(rejects([&] { sylow_compute(Dt, 3, true, 0); }), "a subgroup list that does not match its images is rejected"); }
        { SylowData Dt = load(sd + "a5.txt"); Dt.order = 120;
          CHECK(rejects([&] { sylow_compute(Dt, 3, true, 0); }), "a subgroup that is not Sylow for the stated order is rejected"); }
        CHECK(rejects([&] { parse_sylow("order 60\nprime 2\n4 2\n"); }), "a truncated file is rejected");
    }
    printf("t_pcgroup: %d checks, %d failures\n", checks, fails);
    return fails ? 1 : 0;
}
