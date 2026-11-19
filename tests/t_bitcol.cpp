// t_bitcol.cpp -- the bit-tree working column (bitcol.hpp) against a std::map reference, over F_2, F_p
// (Montgomery), multi-modular Q and exact int64 rationals: random additions with heavy cancellation, the
// pivot after every step, and the drained column.
#include "perseus/bitcol.hpp"
#include <cstdio>
#include <map>
#include <random>
#include <string>

using namespace perseus;

static int fails = 0, checks = 0;
static void check(bool ok, const std::string& what) { ++checks; if (!ok) { ++fails; if (fails < 20) printf("FAIL %s\n", what.c_str()); } }

template <class R>
static void trial(const R& ring, const char* name, std::mt19937_64& g, size_t n, int steps, int64_t span) {
    typedef typename R::T T;
    BitCol<R> c; c.init(n);
    std::map<uint32_t, T> ref;
    auto eq = [&](T a, T b) { return ring.is_zero(ring.add(a, ring.neg(b))); };
    for (int s = 0; s < steps; ++s) {
        // rows cluster in a window so that duplicates and cancellations are frequent
        const size_t win = (g() % 3 == 0) ? n : std::min<size_t>(n, 200);
        uint32_t r = (uint32_t)(g() % win);
        int64_t v = (int64_t)(g() % (2 * span + 1)) - span; if (v == 0) v = 1;
        T x = ring.from_int(v);
        if (ring.is_zero(x)) continue;
        c.add(ring, r, x);
        auto it = ref.find(r);
        if (it == ref.end()) ref[r] = x;
        else { T y = ring.add(it->second, x); if (ring.is_zero(y)) ref.erase(it); else it->second = y; }
        check(c.count == ref.size(), std::string(name) + ": count");
        if (!ref.empty()) {
            uint32_t m = c.min_row();
            check(m == ref.begin()->first, std::string(name) + ": min row");
            check(eq(c.coef(ring, m), ref.begin()->second), std::string(name) + ": pivot coefficient");
        } else check(c.empty(), std::string(name) + ": empty");
        if (s % 97 == 0 && !ref.empty()) {          // erase the pivot, as the reduction does
            uint32_t m = c.min_row(); c.erase(m); ref.erase(ref.begin());
        }
    }
    std::vector<std::pair<uint32_t, T>> out;
    c.drain(ring, [&](uint32_t r, T x) { out.push_back({r, x}); });
    bool same = out.size() == ref.size(); size_t i = 0;
    for (auto& kv : ref) { if (!same) break; same = out[i].first == kv.first && eq(out[i].second, kv.second); ++i; }
    check(same, std::string(name) + ": drained column");
    check(c.empty() && c.count == 0, std::string(name) + ": empty after drain");
    for (auto& w : c.lev) for (uint64_t x : w) if (x) { check(false, std::string(name) + ": stale bits"); return; }
}

int main(int argc, char** argv) {
    int trials = argc > 1 ? atoi(argv[1]) : 200;
    std::mt19937_64 g(argc > 2 ? strtoull(argv[2], nullptr, 10) : 7);
    GF2 f2; Mont<uint32_t, 1> f3({3}), f7({7}); Mont<uint64_t, 2> qm({4611686018427387847ULL, 4611686018427387817ULL}); Q64 q;
    const size_t sizes[] = {1, 2, 63, 64, 65, 4095, 4096, 4097, 262145, 300000};
    for (int t = 0; t < trials; ++t) {
        size_t n = sizes[t % 10];
        trial(f2, "F2", g, n, 2000, 1);
        trial(f3, "F3", g, n, 2000, 5);
        trial(f7, "F7", g, n, 2000, 20);
        trial(qm, "Q multimodular", g, n, 1000, 1000);
        trial(q, "Q64", g, n, 500, 30);
    }
    printf("%d checks, %d failures\n", checks, fails);
    return fails != 0;
}
