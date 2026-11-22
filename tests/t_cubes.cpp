// The packed sparse cubical builder against the ordered-map builder: identical complexes (dimensions,
// births, numbering, faces and signs) on random cube sets in dimensions 1-5, including repeated
// anchors, negative and far-apart coordinates, and boxes too large to pack (map fallback).
#include "perseus/complex.hpp"

#include <cstdio>
#include <cstdlib>
#include <random>

using namespace perseus;

static bool same(const FilteredComplex& a, const FilteredComplex& b) {
    return a.dim == b.dim && a.birth == b.birth && a.off == b.off && a.face == b.face && a.coef == b.coef;
}

int main(int argc, char** argv) {
    int trials = argc > 1 ? std::atoi(argv[1]) : 300;
    std::mt19937_64 rng(20261122);
    int fails = 0, fallback = 0;
    for (int t = 0; t < trials; ++t) {
        int d = 1 + (int)(rng() % 5);
        int extent = 2 + (int)(rng() % (d <= 2 ? 12 : d == 3 ? 6 : 4));
        int64_t offset = (int64_t)(rng() % 2001) - 1000;
        int64_t far = (t % 7 == 0 && d >= 2) ? (int64_t)1 << 40 : 0;   // a far cube on every axis forces the map fallback
        size_t n = 1 + rng() % (d <= 2 ? 80 : 40);
        std::vector<std::vector<int64_t>> an; std::vector<double> bs;
        for (size_t i = 0; i < n; ++i) {
            std::vector<int64_t> a(d);
            for (int k = 0; k < d; ++k) a[k] = offset + (int64_t)(rng() % extent);
            if (far && i == 0) for (int k = 0; k < d; ++k) a[k] += far;
            an.push_back(a); bs.push_back((double)(rng() % 5) * 0.5);
        }
        // the builders take distinct anchors (read_cubes merges repeats); merge here the same way
        std::vector<size_t> ord(an.size()); for (size_t i = 0; i < ord.size(); ++i) ord[i] = i;
        std::sort(ord.begin(), ord.end(), [&](size_t x, size_t y) { return an[x] != an[y] ? an[x] < an[y] : bs[x] < bs[y]; });
        std::vector<std::vector<int64_t>> A; std::vector<double> B;
        for (size_t i = 0; i < ord.size(); ++i) {
            if (i && an[ord[i]] == an[ord[i - 1]]) continue;
            A.push_back(an[ord[i]]); B.push_back(bs[ord[i]]);
        }
        fallback += far != 0;
        FilteredComplex p = build_cubical_sparse(d, A, B), m = build_cubical_sparse_map(d, A, B);
        std::string err = p.validate(true);
        if (!same(p, m) || !err.empty()) {
            ++fails;
            fprintf(stderr, "trial %d (d=%d, %zu cubes): %s\n", t, d, A.size(), err.empty() ? "complexes differ" : err.c_str());
        }
    }
    printf("t_cubes: %d trials (%d through the map fallback), %d failures\n", trials, fallback, fails);
    return fails ? 1 : 0;
}
