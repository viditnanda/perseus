#include "perseus/explicit_engine.hpp"
#include "perseus/bigint.hpp"
#include <cstdio>
#include <random>
using namespace perseus;

static FilteredComplex rp2() {
    // 6-vertex RP^2, all born at 0
    int T[10][3] = {{0,1,2},{0,2,3},{0,3,4},{0,4,5},{0,5,1},{1,2,4},{2,3,5},{3,4,1},{4,5,2},{5,1,3}};
    SimplicialBuilder B; for (auto& t : T) B.add({t[0], t[1], t[2]}, 0); return B.build();
}
template <class R> void show(const char* name, const FilteredComplex& C, const R& ring) {
    Prepared P(C); RunResult r = run_explicit(P, ring, 0, P.D);
    auto bc = barcode(P, r);
    printf("%s:", name);
    for (size_t k = 0; k < bc.size(); ++k) printf(" b%zu=%zu", k, bc[k].size());
    printf(" nonunit=%zu unlucky=%d", r.nonunit.size(), (int)r.unlucky);
    for (auto& nu : r.nonunit) { printf(" [k=%d", nu.k); for (auto x : nu.res) printf(" %llu", (unsigned long long)x); printf("]"); }
    printf("\n");
}
int main() {
    FilteredComplex C = rp2();
    printf("validate: '%s'\n", C.validate(true).c_str());
    show("RP2/F2", C, GF2());
    show("RP2/F3", C, Mont<uint32_t, 1>({3}));
    show("RP2/Q", C, Mont<uint32_t, 2>(big_primes(2)));
    show("RP2/Q64", C, Mont<uint64_t, 2>({(1ULL << 61) - 1, 4611686018427387847ULL}));
    // torus via cubical periodic grid 3x3
    FilteredComplex T = build_cubical({3, 3}, {true, true}, std::vector<double>(9, 1.0));
    printf("torus validate '%s' cells=%zu\n", T.validate(true).c_str(), T.size());
    show("T2/Q", T, Mont<uint32_t, 2>(big_primes(2)));
    // bigint sanity
    BigInt a = BigInt::from_u64(2147483647ULL) * BigInt::from_u64(2147483629ULL) * BigInt::from_u64(4294967291ULL);
    BigInt q, r; divmod(a, BigInt::from_u64(2147483629ULL), q, r);
    printf("bigint %s / p = %s rem %s\n", a.str().c_str(), q.str().c_str(), r.str().c_str());
    BigInt x("" == nullptr ? 0 : 123456789); (void)x;
    std::mt19937_64 g(1);
    int bad = 0;
    for (int it = 0; it < 20000; ++it) {
        int la = 1 + g() % 6, lb = 1 + g() % 4;
        BigInt u, v; u.a.resize(la); v.a.resize(lb);
        for (auto& w : u.a) w = (uint32_t)g(); for (auto& w : v.a) w = (uint32_t)g();
        if (g() % 3 == 0) v.a.back() = 1 + (uint32_t)(g() % 7);
        u.trim(); v.trim(); if (v.is_zero()) continue;
        BigInt qq, rr; divmod(u, v, qq, rr);
        if (!(qq * v + rr == u) || !(rr < v)) ++bad;
    }
    printf("divmod random failures: %d\n", bad);
}
