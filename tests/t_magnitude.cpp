// Tests for integral magnitude homology: against an independent field reduction of the unreduced
// complexes, on trees (diagonal), cycles, random graphs, and the Kaneta-Yoshinaga graph of RP^2 (Z/2 in
// MH_{3,4}).   g++ -std=c++17 -O2 -pthread -Iinclude tests/t_magnitude.cpp -o t_magnitude && ./t_magnitude
#include "perseus/magnitude.hpp"
#include "perseus/explicit_engine.hpp"
#include "perseus/certify.hpp"
#include <cstdio>
#include <random>
#include <set>

using namespace perseus;
static int fails = 0, checks = 0;
static void check(bool ok, const std::string& what) { ++checks; if (!ok) { ++fails; printf("FAIL %s\n", what.c_str()); } }

// Betti numbers over F_p of MC_{*,l} summed over all (s, t), from the unreduced complexes (explicit engine)
static std::vector<std::vector<size_t>> field_ranks(const MagGraph& G, int lmax, uint64_t p) {
    std::vector<std::vector<size_t>> r(lmax + 1, std::vector<size_t>(lmax + 1, 0));
    FlatTuples byl;
    for (int s = 0; s < G.n; ++s) for (int t = 0; t < G.n; ++t) {
        magnitude_tuples(G, s, t, lmax, byl);
        for (int l = 0; l <= lmax; ++l) {
            if (byl[l].empty()) continue;
            FilteredComplex C = magnitude_complex(G, byl[l]);
            Prepared P(C); ExplicitBackend be(P, P.D);
            Barcode b = run_field(be, p, 0, P.D, nullptr).bars;
            for (size_t k = 0; k < b.size(); ++k) for (auto& x : b[k]) if (std::isinf(x.death)) r[l][k]++;
        }
    }
    return r;
}
static void compare(const std::string& name, const MagGraph& G, int lmax, int threads = 2) {
    MagResult R = magnitude_homology(G, lmax, threads);
    for (uint64_t p : {2ULL, 3ULL}) {
        auto f = field_ranks(G, lmax, p);
        bool ok = true;
        for (int l = 0; l <= lmax; ++l) for (int k = 0; k <= lmax; ++k) if (magnitude_rank_mod(R, k, l, p) != f[l][k]) ok = false;
        check(ok, name + ": integral answer reduces to the F_" + std::to_string(p) + " ranks");
    }
}

int main() {
    // trees: MH is diagonal (Hepworth-Willerton)
    {
        std::mt19937_64 g(3);
        for (int t = 0; t < 6; ++t) {
            int n = 4 + (int)(g() % 6); std::vector<std::pair<int, int>> E;
            for (int v = 1; v < n; ++v) E.push_back({v, (int)(g() % v)});
            MagGraph G(n, E, false);
            MagResult R = magnitude_homology(G, 5, 2);
            bool diag = true;
            for (int l = 0; l <= 5; ++l) for (int k = 0; k <= 5; ++k) if (k != l && (R.table[l][k].free || !R.table[l][k].torsion.empty())) diag = false;
            check(diag, "tree " + std::to_string(t) + ": magnitude homology is diagonal");
            compare("tree " + std::to_string(t), G, 4);
        }
    }
    // cycles and random graphs (directed and undirected)
    for (int n = 3; n <= 7; ++n) {
        std::vector<std::pair<int, int>> E; for (int i = 0; i < n; ++i) E.push_back({i, (i + 1) % n});
        compare("C_" + std::to_string(n), MagGraph(n, E, false), 5);
        compare("directed C_" + std::to_string(n), MagGraph(n, E, true), 5);
    }
    {
        std::mt19937_64 g(11);
        for (int t = 0; t < 8; ++t) {
            int n = 5 + (int)(g() % 3); std::vector<std::pair<int, int>> E;
            for (int a = 0; a < n; ++a) for (int b = 0; b < n; ++b) if (a != b && g() % 4 == 0) E.push_back({a, b});
            compare("random graph " + std::to_string(t), MagGraph(n, E, t % 2 == 1), 4);
        }
    }
    // Kaneta-Yoshinaga: the Hasse diagram of the face poset of RP^2_6 with a bottom and a top; Z/2 = H_1(RP^2)
    // embeds in MH_{3,4}
    {
        const int F[10][3] = {{1,2,4},{1,2,6},{1,3,4},{1,3,5},{1,5,6},{2,3,5},{2,3,6},{2,4,5},{3,4,6},{4,5,6}};
        std::map<std::vector<int>, int> id; std::vector<std::pair<int, int>> E;
        int bot = 0; id[{}] = bot; int nid = 1;
        std::set<std::vector<int>> faces;
        for (auto& f : F) { std::vector<int> s(f, f + 3); for (int m = 1; m < 8; ++m) { std::vector<int> sub; for (int i = 0; i < 3; ++i) if (m >> i & 1) sub.push_back(s[i]); faces.insert(sub); } }
        for (auto& s : faces) id[s] = nid++;
        int top = nid++;
        for (auto& s : faces) {
            if (s.size() == 1) E.push_back({bot, id[s]});
            if (s.size() == 3) E.push_back({id[s], top});
            for (size_t i = 0; i < s.size() && s.size() > 1; ++i) { auto t = s; t.erase(t.begin() + i); E.push_back({id[t], id[s]}); }
        }
        MagGraph G(nid, E, false);
        MagResult R = magnitude_homology(G, 4, 2, {{bot, top}});
        const MagEntry& e = R.table[4][3];
        check(e.torsion.count("2") && e.torsion.at("2") >= 1, "Kaneta-Yoshinaga graph of RP^2_6: Z/2 in MH_{3,4}(0, 1)");
        printf("KY(RP^2_6), pair (bottom, top): %d vertices; MH_{3,4}: free %zu, torsion", nid, e.free);
        for (auto& kv : e.torsion) printf(" Z/%s^%zu", kv.first.c_str(), kv.second);
        printf(" (%zu generators, %zu after reduction)\n", e.generators, e.residual);
    }
    printf("%d checks, %d failures\n", checks, fails);
    return fails ? 1 : 0;
}
