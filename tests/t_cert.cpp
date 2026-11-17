#include "perseus/explicit_engine.hpp"
#include "perseus/certify.hpp"
#include <cstdio>
#include <random>
using namespace perseus;

static void report(const char* name, const FilteredComplex& C, int maxdim = 99) {
    Prepared P(C); ExplicitBackend be(P, maxdim);
    CertifyOptions o; o.threads = 2;
    auto res = certify(be, o);
    auto& q = res.first; auto& c = res.second;
    printf("%-22s", name);
    for (size_t k = 0; k < q.bars.size(); ++k) printf(" |H%zu|=%zu", k, q.bars[k].size());
    printf("  nonunit=%zu lanes=%d retries=%d unrec=%d unresolved=%zu cand={", c.values.size(), c.lanes, c.retries, c.unreconstructed, c.unresolved.size());
    for (auto& kv : c.candidates) printf("%llu,", (unsigned long long)kv.first);
    printf("} bad={");
    for (auto& b : c.bad) printf("%llu,", (unsigned long long)b.p);
    printf("}  tQ=%.3f trec=%.3f tcls=%.3f\n", c.t_q, c.t_recon, c.t_classify);
}
int main(int argc, char** argv) {
    {   int T[10][3] = {{0,1,2},{0,2,3},{0,3,4},{0,4,5},{0,5,1},{1,2,4},{2,3,5},{3,4,1},{4,5,2},{5,1,3}};
        SimplicialBuilder B; for (auto& t : T) B.add({t[0], t[1], t[2]}, 0); report("RP2", B.build()); }
    {   // Moore space M(Z/m,1): circle (one vertex v, one edge e) and a 2-cell with boundary m*e
        for (int m : {6, 35, 1001}) {
            FilteredComplex C; int v = C.add_cell(0, 0, {}); int e = C.add_cell(1, 1, {}); C.add_cell(2, 2, {{e, m}}); (void)v;
            char nm[64]; snprintf(nm, 64, "Moore(Z/%d)", m); report(nm, C);
        }
    }
    {   // the counterexample: a (dim0,b0), c1 (b1, d = p a)... as a chain complex in degrees 1,2
        FilteredComplex C; int a = C.add_cell(0, 0, {});
        int e = C.add_cell(1, 0, {}); (void)a;
        C.add_cell(2, 1, {{e, 3}}); C.add_cell(2, 2, {{e, 1}});
        report("counterexample p=3", C);
    }
    {   // big torsion: 2-cell with boundary (2^61-1)*e  (prime > 2^31)
        FilteredComplex C; C.add_cell(0, 0, {}); int e = C.add_cell(1, 0, {});
        C.add_cell(2, 1, {{e, 2305843009213693951LL}}); report("Moore(Z/(2^61-1))", C);
    }
    {   // product of two huge primes -> needs bigint factoring of a 120-bit number
        FilteredComplex C; C.add_cell(0, 0, {}); int e = C.add_cell(1, 0, {}); int f = C.add_cell(1, 0, {});
        C.add_cell(2, 1, {{e, 1000000007LL}, {f, 998244353LL}});
        C.add_cell(2, 2, {{e, 1000000009LL}, {f, 1000000007LL}});
        report("2x2 det", C);
    }
    // Linial-Meshulam Y_2(n, m): all edges, m random triangles in random order (filtration = insertion)
    int n = argc > 1 ? atoi(argv[1]) : 40;
    std::mt19937_64 g(argc > 2 ? atoi(argv[2]) : 7);
    std::vector<std::vector<int>> tri;
    for (int a = 0; a < n; ++a) for (int b = a + 1; b < n; ++b) for (int c = b + 1; c < n; ++c) tri.push_back({a, b, c});
    std::shuffle(tri.begin(), tri.end(), g);
    SimplicialBuilder B;
    for (int a = 0; a < n; ++a) for (int b = a + 1; b < n; ++b) B.add({a, b}, 0);
    for (size_t t = 0; t < tri.size(); ++t) B.add(tri[t], 1 + (double)t);
    char nm[64]; snprintf(nm, 64, "LM full n=%d", n);
    report(nm, B.build(), 1);
}
