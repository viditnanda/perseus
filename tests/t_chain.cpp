// Tests for chain complexes over Z (validation, algebraic Morse reduction, integral outputs) and finite
// simplicial sets (identities, normalised chains, nerves, BG, quotients, suspensions).
//   g++ -std=c++17 -O2 -pthread -Iinclude tests/t_chain.cpp -o t_chain && ./t_chain
#include "perseus/chain.hpp"
#include "perseus/simplicial_set.hpp"
#include "perseus/explicit_engine.hpp"
#include "perseus/certify.hpp"
#include <cstdio>
#include <random>
#include <string>

using namespace perseus;
static int fails = 0, checks = 0;
static void check(bool ok, const std::string& what) { ++checks; if (!ok) { ++fails; printf("FAIL %s\n", what.c_str()); } }

static std::string show(const IntegralHomology& H) {
    std::string s;
    for (size_t k = 0; k < H.free.size(); ++k) {
        s += " H" + std::to_string(k) + "=Z^" + std::to_string(H.free[k]);
        for (auto& t : H.torsion[k]) s += "+Z/" + t.str();
    }
    return s;
}
// expected: per degree (free rank, torsion list)
static bool same(const IntegralHomology& H, const std::vector<std::pair<size_t, std::vector<int>>>& want) {
    if (H.free.size() < want.size()) return false;
    for (size_t k = 0; k < want.size(); ++k) {
        if (H.free[k] != want[k].first || H.torsion[k].size() != want[k].second.size()) return false;
        for (size_t i = 0; i < want[k].second.size(); ++i) if (!(H.torsion[k][i] == BigInt(want[k].second[i]))) return false;
    }
    return true;
}
static IntegralHomology homology_of(const SimplicialSet& S, int maxdim) {
    std::string e = S.validate(); if (!e.empty()) { printf("  invalid: %s\n", e.c_str()); return IntegralHomology(); }
    FilteredComplex C = S.chains();
    e = validate_chain_complex(C); if (!e.empty()) { printf("  invalid chains: %s\n", e.c_str()); return IntegralHomology(); }
    ReducedComplex R = algebraic_morse(C);
    return integral_homology(R.C, maxdim);
}
static void expect(const std::string& name, const IntegralHomology& H, const std::vector<std::pair<size_t, std::vector<int>>>& want) {
    check(same(H, want), name + ":" + show(H));
}

// RP^2, 6-vertex triangulation (hemi-icosahedron)
static const int RP2[10][3] = {{1,2,4},{1,2,6},{1,3,4},{1,3,5},{1,5,6},{2,3,5},{2,3,6},{2,4,5},{3,4,6},{4,5,6}};

int main(int argc, char** argv) {
    int trials = argc > 1 ? atoi(argv[1]) : 200;
    // ---- Smith forms ----
    {
        auto d = zlin::invariant_factors({{2, 4}, {6, 8}});
        check(d.size() == 2 && d[0] == BigInt(2) && d[1] == BigInt(4), "SNF [[2,4],[6,8]] = diag(2,4)");
        auto d2 = zlin::invariant_factors({{0, 0}, {0, 0}});
        check(d2.empty(), "SNF of zero");
        // overflow path: entries near 2^62
        int64_t big = (int64_t)1 << 61;
        auto d3 = zlin::invariant_factors({{big, big + 2}, {big + 6, big + 4}});
        // det = big(big+4) - (big+2)(big+6) = -4 big - 12 => factors gcd(entries)=2 and |det|/2
        BigInt det = BigInt(big) * BigInt(big + 4) - BigInt(big + 2) * BigInt(big + 6);
        check(d3.size() == 2 && d3[0] == BigInt(2) && d3[1] == det.abs() / BigInt(2), "SNF with int64 overflow (BigInt fallback)");
        auto K = zlin::kernel_basis({{1, 2, 3}, {2, 4, 6}}, 3);
        check(K.size() == 2, "kernel of a rank-1 3-column matrix has rank 2");
    }
    // ---- RP^2: integral homology, Morse reduction, two-stage image ----
    {
        SimplicialBuilder B; for (auto& t : RP2) B.add({t[0], t[1], t[2]}, 0.0);
        FilteredComplex C = B.build();
        check(validate_chain_complex(C).empty(), "RP^2 chains valid");
        expect("RP^2 (simplicial)", integral_homology(C, 2), {{1, {}}, {0, {2}}, {0, {}}});
        ReducedComplex R = algebraic_morse(C);
        expect("RP^2 after Morse reduction", integral_homology(R.C, 2), {{1, {}}, {0, {2}}, {0, {}}});
        check(R.C.size() <= 3, "RP^2 reduces to at most 3 cells (got " + std::to_string(R.C.size()) + ")");
        // two-stage: the triangle 1-2-3 (an empty triangle: a projective line) at s = 0, all of RP^2 at t = 1
        SimplicialBuilder B2; B2.add({1, 2}, 0); B2.add({2, 3}, 0); B2.add({1, 3}, 0);
        for (auto& t : RP2) B2.add({t[0], t[1], t[2]}, 1.0);
        FilteredComplex C2 = B2.build();
        IntegralHomology im = integral_image(C2, 0, 1, 2);
        expect("image H(RP^1) -> H(RP^2)", im, {{1, {}}, {0, {2}}, {0, {}}});
        ReducedComplex R2 = algebraic_morse(C2);
        expect("image after Morse reduction", integral_image(R2.C, 0, 1, 2), {{1, {}}, {0, {2}}, {0, {}}});
        // the boundary of a Moebius band maps to twice the core: the circle 1-2-3 in a cone over RP^2 dies
        SimplicialBuilder B3; B3.add({1, 2}, 0); B3.add({2, 3}, 0); B3.add({1, 3}, 0);
        for (auto& t : RP2) B3.add({t[0], t[1], t[2], 7}, 1.0);
        expect("image into the cone is trivial", integral_image(B3.build(), 0, 1, 2), {{1, {}}, {0, {}}, {0, {}}});
        // H_1(S^1) -> H_1(RP^2) for the edge-path 1-2-4-1 (a face boundary): image 0
        SimplicialBuilder B4; B4.add({1, 2}, 0); B4.add({2, 4}, 0); B4.add({1, 4}, 0);
        for (auto& t : RP2) B4.add({t[0], t[1], t[2]}, 1.0);
        expect("image of a contractible loop", integral_image(B4.build(), 0, 1, 2), {{1, {}}, {0, {}}, {0, {}}});
    }
    // ---- simplicial sets ----
    for (int n = 1; n <= 6; ++n) {
        std::vector<std::pair<size_t, std::vector<int>>> w(n + 1, {0, {}}); w[0].first = 1; w[n].first += 1;
        expect("minimal S^" + std::to_string(n), homology_of(minimal_sphere(n), n), w);
    }
    {   // torus and Klein bottle: one vertex, three edges, two triangles (faces d_0, d_1, d_2)
        for (int klein = 0; klein < 2; ++klein) {
            SimplicialSet S; int32_t v = S.add(0, 0, {});
            SSimp V; V.nd = v;
            int32_t a = S.add(1, 0, {V, V}), b = S.add(1, 0, {V, V}), c = S.add(1, 0, {V, V});
            SSimp A{a, {}}, Bb{b, {}}, Cc{c, {}};
            S.add(2, 0, {Bb, Cc, A});                       // d U = b - c + a
            if (!klein) S.add(2, 0, {A, Cc, Bb});           // d L = a - c + b
            else S.add(2, 0, {A, Bb, Cc});                  // d L = a - b + c
            if (!klein) expect("torus", homology_of(S, 2), {{1, {}}, {2, {}}, {1, {}}});
            else expect("Klein bottle", homology_of(S, 2), {{1, {}}, {1, {2}}, {0, {}}});
        }
        // RP^2: one vertex, one edge a, one triangle with faces (a, s_0 *, a)
        SimplicialSet S; int32_t v = S.add(0, 0, {}); SSimp V{v, {}};
        int32_t a = S.add(1, 0, {V, V}); SSimp A{a, {}};
        S.add(2, 0, {A, make_degenerate(v, {0}), A});
        expect("minimal RP^2", homology_of(S, 2), {{1, {}}, {0, {2}}, {0, {}}});
        // a broken simplicial set must be rejected
        SimplicialSet T; int32_t u = T.add(0, 0, {}), u2 = T.add(0, 0, {});
        int32_t e1 = T.add(1, 0, {SSimp{u, {}}, SSimp{u2, {}}});
        T.add(2, 0, {SSimp{e1, {}}, SSimp{e1, {}}, SSimp{e1, {}}});     // d_0 d_1 != d_0 d_0
        check(!T.validate().empty(), "an inconsistent simplicial set is rejected");
    }
    {   // BG for Z/2, Z/3, S_3
        auto cyclic = [](int n) { std::vector<std::vector<int>> m(n, std::vector<int>(n)); for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) m[i][j] = (i + j) % n; return m; };
        int K = 7;
        std::vector<std::pair<size_t, std::vector<int>>> w2(K, {0, {}}), w3(K, {0, {}});
        w2[0].first = w3[0].first = 1;
        for (int k = 1; k < K; k += 2) { w2[k].second = {2}; w3[k].second = {3}; }
        expect("BZ/2", homology_of(bar_construction(cyclic(2), 0, K), K - 1), w2);
        expect("BZ/3", homology_of(bar_construction(cyclic(3), 0, K), K - 1), w3);
        // S_3 as permutations of {0,1,2}
        std::vector<std::vector<int>> perms = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
        std::vector<std::vector<int>> mult(6, std::vector<int>(6));
        for (int i = 0; i < 6; ++i) for (int j = 0; j < 6; ++j) {
            std::vector<int> p(3); for (int x = 0; x < 3; ++x) p[x] = perms[i][perms[j][x]];
            for (int r = 0; r < 6; ++r) if (perms[r] == p) mult[i][j] = r;
        }
        // H_n(S_3; Z): Z, Z/2, 0, Z/6, 0, Z/2 (n = 0..5)
        expect("BS_3", homology_of(bar_construction(mult, 0, 6), 5), {{1, {}}, {0, {2}}, {0, {}}, {0, {6}}, {0, {}}, {0, {2}}});
    }
    {   // a category with a non-identity endomorphism: the monoid {1, e}, e e = e (e is a zero): contractible
        Category C; C.objects = 1; C.src = {0, 0}; C.tgt = {0, 0}; C.identity = {0};
        C.compose = [](int g, int f) { return (g == 1 || f == 1) ? 1 : 0; };
        expect("B{1,e} (e idempotent)", homology_of(nerve(C, 6), 5), {{1, {}}, {0, {}}, {0, {}}, {0, {}}, {0, {}}, {0, {}}});
        // two objects, a non-identity endomorphism f of x with f f = 1 (a Z/2 inside) and an arrow x -> y
        // absorbing it (g f = g): the nerve is contractible (y is terminal)
        Category D; D.objects = 2;
        // morphisms: 0 = 1_x, 1 = 1_y, 2 = f: x -> x, 3 = g: x -> y
        D.src = {0, 1, 0, 0}; D.tgt = {0, 1, 0, 1}; D.identity = {0, 1};
        D.compose = [](int g, int f) { if (g == 0 || g == 1) return f; if (f == 0 || f == 1) return g; if (g == 2 && f == 2) return 0; if (g == 3 && f == 2) return 3; throw std::logic_error("not composable"); };
        expect("category with an involution and a terminal object", homology_of(nerve(D, 6), 5), {{1, {}}, {0, {}}, {0, {}}, {0, {}}, {0, {}}, {0, {}}});
    }
    {   // suspensions: Sigma S^1 = S^2, Sigma RP^2 has H_2 = Z/2
        expect("Sigma S^1", homology_of(suspension(minimal_sphere(1)), 3), {{1, {}}, {0, {}}, {1, {}}, {0, {}}});
        SimplicialSet S; int32_t v = S.add(0, 0, {}); SSimp V{v, {}};
        int32_t a = S.add(1, 0, {V, V}); SSimp A{a, {}};
        S.add(2, 0, {A, make_degenerate(v, {0}), A});
        expect("Sigma RP^2", homology_of(suspension(S), 3), {{1, {}}, {0, {}}, {0, {2}}, {0, {}}});
        expect("Sigma^2 RP^2", homology_of(suspension(suspension(S)), 4), {{1, {}}, {0, {}}, {0, {}}, {0, {2}}, {0, {}}});
    }
    {   // nerve of a poset = order complex, with weights: compare certified barcodes
        std::mt19937_64 g(5);
        for (int t = 0; t < 30; ++t) {
            int n = 4 + (int)(g() % 4);
            std::vector<std::vector<char>> le(n, std::vector<char>(n, 0));
            for (int i = 0; i < n; ++i) le[i][i] = 1;
            for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) if (g() % 3 == 0) le[i][j] = 1;
            for (int k = 0; k < n; ++k) for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) if (le[i][k] && le[k][j]) le[i][j] = 1;
            Category C; C.objects = n; std::map<std::pair<int, int>, int> mid; std::vector<double> wt;
            for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) if (le[i][j]) { mid[{i, j}] = (int)C.src.size(); C.src.push_back(i); C.tgt.push_back(j); wt.push_back(i == j ? 0.0 : (double)(1 + g() % 4)); }
            for (int i = 0; i < n; ++i) C.identity.push_back(mid[{i, i}]);
            C.weight = wt;
            C.compose = [&](int gg, int ff) { return mid.at({C.src[ff], C.tgt[gg]}); };
            SimplicialSet S = nerve(C, n);
            check(S.validate().empty(), "poset nerve valid");
            // order complex: chains x_0 < ... < x_k, birth = max weight over pairs x_i < x_j (runs)
            SimplicialBuilder OB;
            for (size_t x = 0; x < S.size(); ++x) {
                // recover the chain of objects from the nerve simplex: walk the vertices via faces
                std::vector<int> objs; SSimp X{(int32_t)x, {}}; int d = S.dim[x];
                for (int v = 0; v <= d; ++v) {     // vertex v = d_0^{v} d_{d-v}... : use the last-vertex faces
                    SSimp y = X; for (int r = d; r > v; --r) y = S.face(y, r); for (int r = 0; r < v; ++r) y = S.face(y, 0);
                    objs.push_back(std::stoi(S.label[y.nd].substr(1)));
                }
                OB.add(objs, S.birth[x]);
            }
            FilteredComplex A = S.chains(), Bc = OB.build();
            Prepared PA(A), PB(Bc); ExplicitBackend ea(PA, n - 1), eb(PB, n - 1);
            auto ra = certify(ea), rb = certify(eb);
            check(ra.first.bars == rb.first.bars, "poset nerve == order complex (trial " + std::to_string(t) + ")");
        }
    }
    // ---- Morse reduction preserves barcodes (all fields, Q, zero-length bars), bad primes, integral homology ----
    {
        std::mt19937_64 g(17);
        int used_q = 0;
        for (int t = 0; t < trials; ++t) {
            // random filtered complexes: simplicial (few birth levels) or a random simplicial set (BG-like sums)
            FilteredComplex C;
            if (t % 3 != 2) {
                SimplicialBuilder B; int n = 6 + (int)(g() % 5), levels = 1 + (int)(g() % 4);
                int m = 8 + (int)(g() % 20);
                for (int i = 0; i < m; ++i) {
                    int k = 2 + (int)(g() % 3); std::vector<int> s; while ((int)s.size() < k) { int v = (int)(g() % n); if (std::find(s.begin(), s.end(), v) == s.end()) s.push_back(v); }
                    B.add(s, (double)(g() % levels));
                }
                C = B.build();
            } else {
                // RP^2 / Klein / torus pieces glued with random births: use the minimal RP^2 with a filtration
                SimplicialSet S; int32_t v = S.add(0, 0, {}); SSimp V{v, {}};
                int32_t a = S.add(1, (double)(g() % 2), {V, V}); SSimp A{a, {}};
                S.add(2, 1 + (double)(g() % 2), {A, make_degenerate(v, {0}), A});
                S = suspension(S);
                C = S.chains();
            }
            check(validate_chain_complex(C).empty(), "random complex valid");
            ReducedComplex R = algebraic_morse(C);
            Prepared P0(C), P1(R.C);
            int md = std::max(0, P0.D);
            for (bool zero : {false, true}) {
                ExplicitBackend e0(P0, md), e1(P1, std::min(md, std::max(0, P1.D)));
                e0.keep_zero_length = zero; e1.keep_zero_length = zero;
                for (uint64_t p : {2ULL, 3ULL, 5ULL}) {
                    Barcode a = run_field(e0, p, 0, md, nullptr).bars, b = run_field(e1, p, 0, std::min(md, std::max(0, P1.D)), nullptr).bars;
                    b.resize(a.size());
                    if (zero) for (auto& z : R.zero) if (z.first < (int)b.size()) b[z.first].push_back({z.second, z.second});
                    for (auto& v : b) std::sort(v.begin(), v.end());
                    check(a == b, "Morse reduction preserves the F_" + std::to_string(p) + " barcode" + (zero ? " (zero-length)" : "") + " trial " + std::to_string(t));
                }
            }
            ExplicitBackend e0(P0, md), e1(P1, std::min(md, std::max(0, P1.D)));
            auto q0 = certify(e0), q1 = certify(e1);
            Barcode b1 = q1.first.bars; b1.resize(q0.first.bars.size());
            std::string bp0, bp1; for (auto& b : q0.second.bad) bp0 += b.str() + " "; for (auto& b : q1.second.bad) bp1 += b.str() + " ";
            check(q0.first.bars == b1 && bp0 == bp1, "Morse reduction preserves the Q barcode and bad primes, trial " + std::to_string(t) + " [" + bp0 + "] [" + bp1 + "]");
            if (!bp0.empty()) ++used_q;
            IntegralHomology h0 = integral_homology(C, md), h1 = integral_homology(R.C, md);
            check(h0.free == h1.free && show(h0) == show(h1), "Morse reduction preserves integral homology, trial " + std::to_string(t));
            // images for every pair of levels
            std::vector<double> lv(C.birth.begin(), C.birth.end()); std::sort(lv.begin(), lv.end()); lv.erase(std::unique(lv.begin(), lv.end()), lv.end());
            for (size_t i = 0; i < lv.size(); ++i) for (size_t j = i; j < lv.size(); ++j)
                check(show(integral_image(C, lv[i], lv[j], md)) == show(integral_image(R.C, lv[i], lv[j], md)), "Morse reduction preserves images, trial " + std::to_string(t));
        }
        printf("(%d random complexes had bad primes)\n", used_q);
    }
    printf("%d checks, %d failures\n", checks, fails);
    return fails ? 1 : 0;
}
