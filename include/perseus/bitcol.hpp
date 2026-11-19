// bitcol.hpp -- a working column for long column reductions: a 64-ary tree of bit words over the row
// indices, plus a dense coefficient array for rings other than F_2.
//
// A reduction adds multiples of stored columns to a working column and repeatedly needs its pivot (here: the
// smallest row). A sorted vector costs a full merge per addition; a binary heap keeps duplicate rows until
// they surface and needs a sort to store the result. The tree cancels a duplicate the moment it arrives
// (F_2: a bit flip; other rings: an addition in the coefficient array, clearing the bit on zero) and finds
// the smallest row by descending log_64(n) words. This is the "bit tree" pivot column of PHAT (Bauer,
// Kerber, Reininghaus, Wagner, J. Symbolic Comput. 78, 2017), extended with coefficients.
//
// Memory: n/64 + n/64^2 + ... words, and n coefficients for rings other than F_2, allocated once per degree;
// a column is emptied by draining it, which touches only its own entries.
#pragma once
#include <cstdint>
#include <cstddef>
#include <vector>
#include <type_traits>
#include <cstdlib>
#include "ring.hpp"

namespace perseus {

template <class R>
struct BitCol {
    typedef typename R::T T;
    static constexpr bool F2 = std::is_same<R, GF2>::value;
    std::vector<std::vector<uint64_t>> lev;     // lev[0]: one bit per row; lev[l+1]: one bit per nonzero word of lev[l]
    std::vector<T> val;                         // coefficients (unused over F_2)
    size_t count = 0;                           // number of nonzero entries

    void init(size_t n) {
        lev.clear(); count = 0;
        size_t w = (n + 63) / 64; if (w == 0) w = 1;
        while (true) { lev.emplace_back(w, 0); if (w == 1) break; w = (w + 63) / 64; }
        if (!F2) val.assign(n, T());
    }
    bool empty() const { return count == 0; }
    bool has(uint32_t r) const { return (lev[0][r >> 6] >> (r & 63)) & 1; }
    void set_bit(uint32_t r) {
        for (size_t l = 0; l < lev.size(); ++l) {
            uint64_t& w = lev[l][r >> 6]; bool was = w != 0; w |= 1ULL << (r & 63);
            if (was) return;
            r >>= 6;
        }
    }
    void clear_bit(uint32_t r) {
        for (size_t l = 0; l < lev.size(); ++l) {
            uint64_t& w = lev[l][r >> 6]; w &= ~(1ULL << (r & 63));
            if (w) return;
            r >>= 6;
        }
    }
    // add c (nonzero) at row r
    void add(const R& ring, uint32_t r, T c) {
        if (has(r)) {
            if (F2) { clear_bit(r); --count; return; }
            T s = ring.add(val[r], c);
            if (ring.is_zero(s)) { clear_bit(r); --count; } else val[r] = s;
        } else {
            set_bit(r); ++count;
            if (!F2) val[r] = c;
        }
    }
    T coef(const R& ring, uint32_t r) const { if (F2) return ring.one(); return val[r]; }
    void erase(uint32_t r) { clear_bit(r); --count; }
    // the smallest row; the column must not be empty
    uint32_t min_row() const {
        size_t i = 0;
        for (size_t l = lev.size(); l-- > 0;) i = (i << 6) | (size_t)__builtin_ctzll(lev[l][i]);
        return (uint32_t)i;
    }
    // remove every entry in increasing row order, passing (row, coefficient) to f
    template <class F> void drain(const R& ring, F f) {
        while (count) { uint32_t r = min_row(); f(r, coef(ring, r)); erase(r); }
    }
};

// Working columns longer than this move from a sorted vector into a BitCol. PERSEUS_LONG_COLUMN overrides it
// (for tests: 0 sends every column through the tree).
inline size_t long_column_threshold() {
    static const size_t v = [] { const char* s = std::getenv("PERSEUS_LONG_COLUMN"); return s ? (size_t)std::strtoull(s, nullptr, 10) : (size_t)512; }();
    return v;
}

} // namespace perseus
