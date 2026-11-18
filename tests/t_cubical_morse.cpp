// Differential test: the discrete-gradient cubical engine against the classic implicit engine.
// Barcodes over F_2, F_3, F_5, Q (exact int64 and two-lane multi-modular), with and without zero-length
// bars, partial runs from degree k > 0 with clearing data, and the full certificate (bad primes).
//   g++ -std=c++17 -O2 -Iinclude tests/t_cubical_morse.cpp -o t_cubical_morse && ./t_cubical_morse [trials] [seed]
#include "perseus/cubical_engine.hpp"
#include "perseus/certify.hpp"
#include <random>
#include <cstdio>
#include <string>

using namespace perseus;

static int fails = 0, checks = 0;
static void check(bool ok, const std::string& what) { ++checks; if (!ok) { ++fails; if (fails <= 30) printf("FAIL %s\n", what.c_str()); } }

static std::string show(const Barcode& b) {
    std::string s;
    for (size_t k = 0; k < b.size(); ++k) { s += " H" + std::to_string(k) + ":"; for (auto& x : b[k]) s += "[" + std::to_string(x.birth) + "," + std::to_string(x.death) + ")"; }
    return s;
}
template <class R>
static void cmp(const CubicalClassic& c, const CubicalMorse& m, const R& ring, const char* rn, bool keep, const std::string& tag) {
    int kmax = m.maxdim();
    auto a = c.run(ring, 0, kmax, nullptr, keep), b = m.run(ring, 0, kmax, nullptr, keep);
    bool ok = a.bars == b.bars && a.unlucky == b.unlucky;
    check(ok, tag + " " + rn + (keep ? " zero" : "") + "\n  classic" + show(a.bars) + "\n  morse  " + show(b.bars));
    // partial runs from degree k, with each engine's own clearing data
    for (int k = 1; k <= kmax && ok; ++k) {
        auto a2 = c.run(ring, k, kmax, &a.clear_out[k - 1], keep), b2 = m.run(ring, k, kmax, &b.clear_out[k - 1], keep);
        bool same = true;
        for (int j = k; j <= kmax; ++j) same = same && a2.bars[j] == b2.bars[j];
        check(same, tag + " " + rn + " from degree " + std::to_string(k));
    }
}

static void run_case(const CubicalGrid& G, int md, const std::string& tag0, int levels = 0, bool absent = false) {
    int d = (int)G.ext.size(); (void)d;
    CubicalClassic C(G, md); CubicalMorse M(G, md);
    check(M.check_columns() == 0 && M.check_order() == 0 && M.check_dd() == 0, tag0 + ": Morse complex self-check");
    std::string tag = tag0 + " d=" + std::to_string(d) + " ext";
    for (auto e : G.ext) tag += " " + std::to_string(e);
    tag += " maxdim " + std::to_string(md) + " levels " + std::to_string(levels) + (absent ? " absent" : "");
    for (bool keep : {false, true}) {
        cmp(C, M, GF2(), "F2", keep, tag);
        cmp(C, M, Mont<uint32_t, 1>({3}), "F3", keep, tag);
        cmp(C, M, Mont<uint32_t, 1>({5}), "F5", keep, tag);
        { Q64 R1, R2; auto a = C.run(R1, 0, M.maxdim(), nullptr, keep), b = M.run(R2, 0, M.maxdim(), nullptr, keep);
          check(a.bars == b.bars && R1.overflow == R2.overflow, tag + " Q64" + (keep ? " zero" : "") + "\n  classic" + show(a.bars) + "\n  morse  " + show(b.bars)); }
        cmp(C, M, Mont<uint64_t, 2>({4611686018427387847ULL, 4611686018427388039ULL}), "Q(2 lanes)", keep, tag);
    }
    // the full certificate in auto and fast mode
    for (const char* mode : {"auto", "fast"}) {
        CertifyOptions o; o.mode = mode; o.seed = "1234";
        auto rc = certify(C, o); auto rm = certify(M, o);
        std::string bc, bm;
        for (auto& b : rc.second.bad) bc += " " + b.str();
        for (auto& b : rm.second.bad) bm += " " + b.str();
        check(rc.first.bars == rm.first.bars && bc == bm && rc.second.barcode_ok == rm.second.barcode_ok, tag + " certify " + mode + " bad primes:" + bc + " vs" + bm);
    }
}

int main(int argc, char** argv) {
    int trials = argc > 1 ? atoi(argv[1]) : 300;
    std::mt19937_64 g(argc > 2 ? atoll(argv[2]) : 1);
    for (int t = 0; t < trials; ++t) {
        int d = 2 + (int)(g() % 3);
        CubicalGrid G;
        int big = argc > 3 ? atoi(argv[3]) : 1;
        int maxside = (d == 2 ? 12 : d == 3 ? 7 : 4) * big;
        // every 8th trial: a larger grid of nearly distinct values, where Morse path sums of +-2 and long
        // V-paths occur (in 4D, 10^4 white noise already has them)
        bool large = t % 8 == 7;
        if (large) maxside = d == 2 ? 120 : d == 3 ? 24 : 11;
        size_t tot = 1;
        for (int i = 0; i < d; ++i) { int64_t e = (large ? maxside / 2 : 1) + (int64_t)(g() % (large ? maxside / 2 + 1 : maxside)); G.ext.push_back(e); tot *= (size_t)e; }
        G.periodic.assign(d, false);
        int levels = (int)(large ? (g() % 2 ? 1000000 : 64) : g() % 4 == 0 ? 1000000 : 2 + g() % 8);   // heavy ties most of the time
        bool absent = g() % 3 == 0;
        bool reals = large || g() % 5 == 0;
        for (size_t i = 0; i < tot; ++i) {
            double v = reals ? std::ldexp((double)(g() % 1000000), -7) : (double)(g() % levels);
            if (absent && g() % 7 == 0) v = std::numeric_limits<double>::infinity();
            G.top.push_back(v);
        }
        int md = g() % 4 == 0 ? (int)(g() % (d + 1)) : -1;
        run_case(G, md, "trial " + std::to_string(t), levels, absent);
    }

    // degenerate shapes: axes of length 1 and 2, very elongated grids, a single cube, grids that are entirely
    // absent, constant, or a checkerboard of absent cubes
    {
        std::vector<std::vector<int64_t>> shapes = {{1, 1}, {1, 9}, {9, 1}, {2, 2}, {300, 2}, {2, 300}, {1, 1, 1}, {1, 7, 5}, {7, 1, 5}, {2, 2, 2},
                                                    {60, 1, 1}, {1, 1, 60}, {2, 2, 90}, {3, 1, 4}, {1, 1, 1, 1}, {1, 1, 1, 9}, {2, 2, 2, 2},
                                                    {3, 1, 4, 1}, {1, 5, 1, 5}, {40, 2, 2, 1}, {2, 2, 2, 25}, {3, 3, 3, 3}};
        const char* kinds[] = {"uniform", "two levels", "constant", "absent", "one cube", "checkerboard", "reals"};
        int si = 0;
        for (auto& sh : shapes) {
            for (int kind = 0; kind < 7; ++kind) {
                CubicalGrid G; G.ext = sh; G.periodic.assign(sh.size(), false);
                size_t tot = 1; for (auto e : sh) tot *= (size_t)e;
                for (size_t i = 0; i < tot; ++i) {
                    double v;
                    switch (kind) {
                        case 0: v = (double)(g() % 5); break;
                        case 1: v = (double)(g() % 2); break;
                        case 2: v = 3.0; break;
                        case 3: v = std::numeric_limits<double>::infinity(); break;
                        case 4: v = i == tot / 2 ? 1.0 : std::numeric_limits<double>::infinity(); break;
                        case 5: { size_t r = i, par = 0; for (auto e : sh) { par += r % (size_t)e; r /= (size_t)e; } v = par % 2 ? std::numeric_limits<double>::infinity() : (double)(g() % 3); } break;
                        default: v = std::ldexp((double)(g() % 100000), -9) - 50.0; break;
                    }
                    G.top.push_back(v);
                }
                for (int md : {-1, 0, 1}) {
                    std::string tag = std::string("shape ") + std::to_string(si) + " " + kinds[kind] + " md " + std::to_string(md);
                    run_case(G, md, tag);
                    CubicalBackend be(G, md);                      // the default route must pick the gradient engine
                    check(std::string(be.name()) == "morse", tag + ": default engine");
                }
            }
            ++si;
        }
    }
    // axes of length 0: the empty complex (this used to crash both implicit engines)
    for (std::vector<int64_t> sh : {std::vector<int64_t>{0, 5}, std::vector<int64_t>{3, 0, 2}, std::vector<int64_t>{2, 2, 0, 2}}) {
        CubicalGrid G; G.ext = sh; G.periodic.assign(sh.size(), false);
        for (const char* e : {"auto", "classic", "morse"}) {
            CubicalBackend be(G, -1, e);
            CertifyOptions o; auto r = certify(be, o);
            bool none = true; for (auto& v : r.first.bars) none = none && v.empty();
            check(std::string(be.name()) == "empty" && none && r.second.bad.empty(), std::string("empty grid, engine ") + e);
        }
    }
    // int64 guard: with the path bound lowered to 1, any path sum of +-2 must make the gradient engine refuse
    // (std::overflow_error) and the default backend fall back to the classic engine with the same barcodes
    {
        int64_t saved = CubicalMorse::path_limit(); CubicalMorse::path_limit() = 1;
        int triggered = 0;
        for (int seed = 0; seed < 40 && triggered < 3; ++seed) {
            std::mt19937_64 h(1000 + seed);
            CubicalGrid G; G.ext = {10, 10, 10, 10}; G.periodic.assign(4, false);
            for (int i = 0; i < 10000; ++i) G.top.push_back(std::ldexp((double)(h() % (1ULL << 40)), -40));
            bool threw = false;
            try { CubicalMorse M(G, -1); } catch (std::overflow_error&) { threw = true; }
            if (!threw) continue;
            ++triggered;
            CubicalBackend be(G, -1);
            check(std::string(be.name()) == "classic", "overflow fallback: default backend is classic");
            CubicalClassic C(G, -1);
            { Q64 R1, R2; check(be.run(R1, 0, 4, nullptr).bars == C.run(R2, 0, 4, nullptr).bars, "overflow fallback: Q barcodes"); }
            bool forced = false;
            try { CubicalBackend bm(G, -1, "morse"); } catch (std::overflow_error&) { forced = true; }
            check(forced, "overflow with --engine morse is an error, not a silent result");
        }
        check(triggered > 0, "overflow guard was exercised");
        CubicalMorse::path_limit() = saved;
        printf("overflow fallback exercised on %d grids\n", triggered);
    }
    printf("%d checks, %d failures\n", checks, fails);
    return fails ? 1 : 0;
}
