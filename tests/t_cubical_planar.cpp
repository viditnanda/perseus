// Differential test: the planar (2D union-find) cubical engine against the classic and discrete-gradient
// engines. Barcodes over F_2, F_3 and Q, with and without zero-length bars, on random grids with heavy ties,
// absent cells, degenerate shapes and structured images; the certificate takes the unimodular path and
// agrees with the classic engine's proved certificate (no bad primes).
//   g++ -std=c++17 -O2 -Iinclude tests/t_cubical_planar.cpp -o t_cubical_planar && ./t_cubical_planar [trials] [seed]
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

static void run_case(const CubicalGrid& G, int md, const std::string& tag) {
    CubicalClassic C(G, md);
    for (bool keep : {false, true}) {
        CubicalPlanar P(G, md, keep);
        auto cmp = [&](auto ring, const char* rn) {
            auto a = C.run(ring, 0, P.maxdim(), nullptr, keep); auto b = P.run(ring, 0, P.maxdim(), nullptr);
            check(a.bars == b.bars, tag + " " + rn + (keep ? " zero" : "") + "\n  classic" + show(a.bars) + "\n  planar " + show(b.bars));
        };
        cmp(GF2(), "F2"); cmp(Mont<uint32_t, 1>({3}), "F3");
        { Q64 R1; auto a = C.run(R1, 0, P.maxdim(), nullptr, keep); auto b = P.run(Q64(), 0, P.maxdim(), nullptr); check(a.bars == b.bars, tag + " Q" + (keep ? " zero" : "")); }
        if (!keep && CubicalMorse::supports(G)) {
            CubicalMorse M(G, md);
            check(M.run(GF2(), 0, M.maxdim(), nullptr).bars == P.run(GF2(), 0, P.maxdim(), nullptr).bars, tag + " planar == morse");
        }
    }
    // certificates: the planar route is unimodular and must agree with the classic proved certificate
    CubicalBackend be(G, md, "auto"); CubicalBackend bc(G, md, "classic");
    CertifyOptions o; auto rp = certify(be, o); auto rc = certify(bc, o);
    check(std::string(be.name()) == "planar" && rp.second.mode == "unimodular" && rp.second.bad.empty() && rp.second.complete(), tag + " certificate path");
    check(rp.first.bars == rc.first.bars && rc.second.bad.empty(), tag + " certificate barcode / classic bad primes " + std::to_string(rc.second.bad.size()));
}

int main(int argc, char** argv) {
    int trials = argc > 1 ? atoi(argv[1]) : 400;
    std::mt19937_64 g(argc > 2 ? atoll(argv[2]) : 1);
    const double INF = std::numeric_limits<double>::infinity();
    for (int t = 0; t < trials; ++t) {
        CubicalGrid G; G.periodic = {false, false};
        bool large = t % 10 == 9;
        int64_t n0 = 1 + g() % (large ? 120 : 14), n1 = 1 + g() % (large ? 120 : 14);
        G.ext = {n0, n1};
        int kind = (int)(g() % 6);
        int levels = 2 + (int)(g() % 6);
        double pabs = g() % 3 == 0 ? (g() % 4) * 0.15 : 0;
        for (int64_t i = 0; i < n0 * n1; ++i) {
            double v;
            int64_t x = i % n0, y = i / n0;
            switch (kind) {
                case 0: v = (double)(g() % levels); break;
                case 1: v = std::ldexp((double)(g() % 1000000), -7); break;
                case 2: v = std::floor(levels * (0.5 + 0.5 * std::sin(0.4 * x) * std::cos(0.3 * y))); break;   // smooth, quantised
                case 3: v = (double)((x / 3 + y / 3) % 2); break;                                                  // checkerboard of blocks
                case 4: v = (double)(std::abs(x - n0 / 2) + std::abs(y - n1 / 2)) * (g() % 2 ? 1 : -1); break;    // rings with ties
                default: v = (g() % 2) ? 0.0 : -0.0; break;                                                        // +-0
            }
            if (pabs > 0 && (double)(g() % 1000) < pabs * 1000) { int a = (int)(g() % 3); v = a == 0 ? INF : a == 1 ? -INF : std::numeric_limits<double>::quiet_NaN(); }
            G.top.push_back(v);
        }
        int md = g() % 4 == 0 ? (int)(g() % 3) : -1;
        run_case(G, md, "trial " + std::to_string(t) + " " + std::to_string(n0) + "x" + std::to_string(n1) + " kind " + std::to_string(kind) + " md " + std::to_string(md));
    }
    // degenerate shapes and fillings
    std::vector<std::pair<int64_t, int64_t>> shapes = {{1, 1}, {1, 2}, {2, 1}, {2, 2}, {1, 30}, {30, 1}, {2, 200}, {200, 2}, {3, 3}};
    for (auto sh : shapes)
        for (int kind = 0; kind < 6; ++kind) {
            CubicalGrid G; G.periodic = {false, false}; G.ext = {sh.first, sh.second};
            int64_t n = sh.first * sh.second;
            for (int64_t i = 0; i < n; ++i) {
                double v = kind == 0 ? (double)(g() % 3) : kind == 1 ? 1.0 : kind == 2 ? INF : kind == 3 ? (i == n / 2 ? 2.0 : INF)
                         : kind == 4 ? ((i % 2) ? INF : (double)(g() % 2)) : std::ldexp((double)(g() % 1000), -3);
                G.top.push_back(v);
            }
            for (int md : {-1, 0, 1}) run_case(G, md, "shape " + std::to_string(sh.first) + "x" + std::to_string(sh.second) + " kind " + std::to_string(kind) + " md " + std::to_string(md));
        }
    // a ring of absent pixels around a present core inside an absent frame: essential H1 classes
    {
        CubicalGrid G; G.periodic = {false, false}; G.ext = {9, 9};
        for (int y = 0; y < 9; ++y) for (int x = 0; x < 9; ++x) {
            int r = std::max(std::abs(x - 4), std::abs(y - 4));
            G.top.push_back(r == 2 || r == 4 ? INF : (double)(r + (x + y) % 2));
        }
        run_case(G, -1, "nested absent rings");
    }
    printf("%d checks, %d failures\n", checks, fails);
    return fails ? 1 : 0;
}
