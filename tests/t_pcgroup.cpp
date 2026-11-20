// t_pcgroup.cpp -- tests for pc_group.hpp: homology and mod-p cohomology rings of finite solvable groups.
//   1. closed forms: cyclic groups, elementary abelian groups, generalised quaternion and dihedral groups, S_3;
//   2. the Kunneth formula on direct products (mixed primes and mixed relative orders included);
//   3. the universal coefficient theorem: the mod-p ring engine against the integral engine, for every group
//      and every prime dividing its order;
//   4. ground truth from HAP and the Green-King database (tests/fixtures/pcgroups/expected.txt);
//   5. invalid inputs are rejected.
#include "perseus/pc_group.hpp"
#include <cstdio>
#include <fstream>
#include <sstream>
#include <functional>
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
            int D = G.T.size() >= 128 ? 7 : 9;
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
    printf("t_pcgroup: %d checks, %d failures\n", checks, fails);
    return fails ? 1 : 0;
}
