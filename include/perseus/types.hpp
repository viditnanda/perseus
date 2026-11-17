// types.hpp -- barcodes and run records shared by all backends.
#pragma once
#include <vector>
#include <cstdint>
#include <limits>
#include "bigint.hpp"

namespace perseus {

struct Bar { double birth, death; };       // death = +inf for essential classes
inline bool operator<(const Bar& a, const Bar& b) { return a.birth != b.birth ? a.birth < b.birth : a.death < b.death; }
inline bool operator==(const Bar& a, const Bar& b) { return a.birth == b.birth && a.death == b.death; }
typedef std::vector<std::vector<Bar>> Barcode;   // per degree, sorted

// a pivot that is not +-1: residues (multi-modular rings) or its exact value (exact rings)
struct Pivot { int k; uint64_t key; double birth, death; std::vector<uint64_t> res; bool exact = false; BigInt num, den; };

template <class Clear>
struct BRun {
    int kmin = 0;
    Barcode bars;                         // degrees kmin..kmax filled
    std::vector<Pivot> nonunit;
    std::vector<Clear> clear_out;         // clear_out[k]: produced by degree k (index valid for k >= kmin)
    bool unlucky = false;
    std::vector<uint64_t> moduli;
    long long additions = 0, columns = 0;
    // a priori data for the probabilistic error bound, per degree k (pairs between k- and (k+1)-cells):
    //   hbits[k]  >= log2 of every minor of the matrix the engine reduces in degree k (Hadamard bound)
    //   pairs[k]  >= the number of pairs in degree k
    std::vector<double> hbits, pairs;
    void note_bound(int k, double h, double p) {
        if ((int)hbits.size() <= k) { hbits.resize(k + 1, 0.0); pairs.resize(k + 1, 0.0); }
        hbits[k] = std::max(hbits[k], h); pairs[k] = std::max(pairs[k], p);
    }
};

} // namespace perseus
