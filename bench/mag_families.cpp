// mag_families.cpp -- certified integral magnitude homology of the Kaneta-Yoshinaga / Sazdanovic-Summers
// torsion graphs: the Hasse diagram of P^ = P + {0^, 1^} for P the face poset of a regular CW structure on
//   RP^n   = boundary of the (n+1)-cube modulo the antipodal map        (Z/2 in H_odd, odd < n)
//   L(p,1) = (C_2p * C_2p * ... * C_2p) modulo the diagonal rotation     (Z/p in H_1; lens space of dim 2n+1)
//   RP^2_6 = the 6-vertex triangulation of RP^2
// By Kaneta-Yoshinaga, MH_{k, r+2}(0^, 1^) = H~_{k-2}(Delta(P)) (r = rank of P), so the torsion appears at
// (k, l) = (3, r + 2). The program checks that, then pushes to longer lengths l > r + 2 for the pair (0^, 1^),
// and, for small graphs, to all (s, t) pairs.
//
//   mag_families list                       the cases
//   mag_families <case> <lmax> [all] [edges-file]
// prints one JSON line: graph size, per-(k, l) groups, generators, residual cells, wall and CPU seconds, peak RSS.
#include "perseus/magnitude.hpp"
#include <sys/resource.h>
#include <cstdio>
#include <set>
#include <fstream>

using namespace perseus;

struct Poset {                                   // a free quotient of a face poset, by cover relations
    std::vector<int> rank;                       // per element
    std::set<std::pair<int, int>> cover;         // (lower, upper)
};

// faces of a simplicial complex given by facets, modulo a free vertex action g (applied `order` times)
static Poset simplicial_quotient(const std::vector<std::vector<int>>& facets, const std::vector<int>& g, int order) {
    std::set<std::vector<int>> faces;
    for (auto& F : facets) {
        int m = (int)F.size();
        for (int mask = 1; mask < (1 << m); ++mask) { std::vector<int> s; for (int i = 0; i < m; ++i) if (mask >> i & 1) s.push_back(F[i]); std::sort(s.begin(), s.end()); faces.insert(s); }
    }
    auto act = [&](std::vector<int> s) { for (int& v : s) v = g[v]; std::sort(s.begin(), s.end()); return s; };
    std::map<std::vector<int>, int> orbit;
    Poset P;
    for (auto& s : faces) {
        if (orbit.count(s)) continue;
        int id = (int)P.rank.size(); P.rank.push_back((int)s.size() - 1);
        std::vector<int> t = s;
        for (int i = 0; i < order; ++i) { orbit[t] = id; t = act(t); }
        if (t != s) throw std::runtime_error("the action does not have the stated order");
    }
    for (auto& s : faces) if (s.size() > 1)
        for (size_t i = 0; i < s.size(); ++i) { std::vector<int> f = s; f.erase(f.begin() + i); P.cover.insert({orbit[f], orbit[s]}); }
    return P;
}

// faces of the boundary of the N-cube ({-1,0,1}^N minus 0) modulo v -> -v
static Poset cube_rp(int N) {
    std::vector<std::vector<int>> F; std::map<std::vector<int>, int> id;
    int total = 1; for (int i = 0; i < N; ++i) total *= 3;
    Poset P;
    auto canon = [](std::vector<int> v) { std::vector<int> w = v; for (int& x : w) x = -x; return std::min(v, w); };
    for (int c = 0; c < total; ++c) {
        std::vector<int> v(N); int x = c, zeros = 0;
        for (int i = 0; i < N; ++i) { v[i] = x % 3 - 1; x /= 3; zeros += v[i] == 0; }
        if (zeros == N) continue;
        auto cv = canon(v);
        if (!id.count(cv)) { id[cv] = (int)P.rank.size(); P.rank.push_back(zeros); }
    }
    for (auto& kv : id) {                        // faces of v: set one free coordinate to +-1
        const std::vector<int>& v = kv.first;
        for (int i = 0; i < N; ++i) if (v[i] == 0) for (int s : {-1, 1}) { std::vector<int> w = v; w[i] = s; P.cover.insert({id[canon(w)], kv.second}); }
    }
    return P;
}

// join of m copies of C_N with the rotation by N/p on every copy (a free Z/p action when N/p >= 2)
static Poset lens(int p, int m, int N) {
    std::vector<std::vector<int>> facets{{}};
    for (int c = 0; c < m; ++c) {
        std::vector<std::vector<int>> nx;
        for (auto& f : facets) for (int i = 0; i < N; ++i) { auto g = f; g.push_back(c * N + i); g.push_back(c * N + (i + 1) % N); nx.push_back(g); }
        facets = nx;
    }
    std::vector<int> g(m * N);
    for (int c = 0; c < m; ++c) for (int i = 0; i < N; ++i) g[c * N + i] = c * N + (i + N / p) % N;
    return simplicial_quotient(facets, g, p);
}

static Poset rp2_6() {
    std::vector<std::vector<int>> T{{0,1,2},{0,2,3},{0,3,4},{0,4,5},{0,1,5},{1,2,4},{2,3,5},{1,3,4},{2,4,5},{1,3,5}};
    std::vector<int> id(6); for (int i = 0; i < 6; ++i) id[i] = i;
    return simplicial_quotient(T, id, 1);
}

// the Kaneta-Yoshinaga graph: Hasse diagram of P + {0^, 1^}; returns (n, edges, bottom, top, rank)
struct KY { int n, bot, top, r; std::vector<std::pair<int, int>> E; };
static KY ky_graph(const Poset& P) {
    KY G; int m = (int)P.rank.size(); G.r = 0;
    for (int x : P.rank) G.r = std::max(G.r, x + 1);           // rank of P = number of levels
    G.bot = m; G.top = m + 1; G.n = m + 2;
    std::vector<char> has_up(m, 0);
    for (auto& c : P.cover) { G.E.push_back(c); has_up[c.first] = 1; }
    for (int x = 0; x < m; ++x) { if (P.rank[x] == 0) G.E.push_back({G.bot, x}); if (!has_up[x]) G.E.push_back({x, G.top}); }
    return G;
}

static Poset make_case(const std::string& c) {
    if (c == "rp2_6") return rp2_6();
    if (c.rfind("rp", 0) == 0) return cube_rp(std::stoi(c.substr(2)) + 1);           // rpN
    if (c.rfind("lens", 0) == 0) {                                                       // lensP_M: L(P,1) of dim 2M-1
        auto u = c.find('_'); int p = std::stoi(c.substr(4, u - 4)), m = std::stoi(c.substr(u + 1));
        return lens(p, m, 2 * p);
    }
    throw std::runtime_error("unknown case " + c);
}

int main(int argc, char** argv) {
    if (argc < 2 || std::string(argv[1]) == "list") {
        puts("rp2_6 rp2 rp3 rp4 rp5 rp6 lens2_2 lens3_2 lens5_2 lens7_2 lens2_3 lens3_3");
        return 0;
    }
    std::string c = argv[1];
    Poset P = make_case(c);
    KY G = ky_graph(P);
    int lmax = argc > 2 ? std::stoi(argv[2]) : G.r + 2;
    bool all = argc > 3 && std::string(argv[3]) == "all";
    if (argc > 4) { std::ofstream f(argv[4]); f << G.n << "\n"; for (auto& e : G.E) f << e.first << " " << e.second << "\n"; }
    MagGraph M(G.n, G.E, false);
    std::vector<std::pair<int, int>> pairs; if (!all) pairs.push_back({G.bot, G.top});
    int threads = 1; if (const char* t = getenv("THREADS")) threads = atoi(t);
    MagResult R = magnitude_homology(M, lmax, threads, pairs);
    struct rusage ru; getrusage(RUSAGE_SELF, &ru);
    printf("{\"case\": \"%s\", \"vertices\": %d, \"edges\": %zu, \"poset_rank\": %d, \"distance\": %d, \"lmax\": %d, \"pairs\": %zu, "
           "\"seconds\": %.3f, \"cpu_seconds\": %.3f, \"peak_mb\": %.1f, \"largest_complex\": %zu, \"groups\": [",
           c.c_str(), G.n, G.E.size(), G.r, M.d(G.bot, G.top), lmax, R.pairs, R.seconds,
           ru.ru_utime.tv_sec + ru.ru_stime.tv_sec + 1e-6 * (ru.ru_utime.tv_usec + ru.ru_stime.tv_usec), ru.ru_maxrss / 1024.0, R.max_complex);
    bool first = true;
    for (int l = 0; l <= lmax; ++l) for (int k = 0; k <= l; ++k) {
        const MagEntry& e = R.table[l][k];
        if (!e.generators) continue;
        printf("%s{\"k\": %d, \"l\": %d, \"free\": %zu, \"torsion\": {", first ? "" : ", ", k, l, e.free);
        bool f2 = true; for (auto& kv : e.torsion) { printf("%s\"%s\": %zu", f2 ? "" : ", ", kv.first.c_str(), kv.second); f2 = false; }
        printf("}, \"generators\": %zu, \"residual\": %zu}", e.generators, e.residual);
        first = false;
    }
    puts("]}");
    return 0;
}
