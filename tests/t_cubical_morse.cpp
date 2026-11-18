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
        CubicalClassic C(G, md); CubicalMorse M(G, md);
        check(M.check_columns() == 0 && M.check_order() == 0 && M.check_dd() == 0, "trial " + std::to_string(t) + ": Morse complex self-check");
        std::string tag = "trial " + std::to_string(t) + " d=" + std::to_string(d) + " ext";
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
    printf("%d checks, %d failures\n", checks, fails);
    return fails ? 1 : 0;
}
