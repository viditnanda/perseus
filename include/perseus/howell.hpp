// howell.hpp -- column sweeps of a boundary matrix over Z/m, m arbitrary, in Howell form.
//
// The columns of d_j arrive in filtration order. The span of the columns seen so far is kept as a strong
// Groebner basis over the principal ideal ring Z/m (a Howell basis, Howell 1986; Storjohann's thesis,
// ch. 4): one basis vector per row r, the row of its lowest nonzero entry (the latest-born row, as in a
// column reduction), whose entry there generates the ideal (d_r), d_r | m. A vector arriving on a row
// already taken is reduced if d_r divides its entry; otherwise the two are combined by an extended gcd,
// the pivot ideal shrinks to (gcd), and both remainders go on. Every pivot vector u with d_r != 1 also
// sends its annihilator multiple (m / d_r) u down, which keeps the basis a Howell basis.
//
// What this gives, for every prime power p^k exactly dividing m (Chinese remainder theorem: Z/m is the
// product of the Z/p^k, and the Howell form of a product is the product of the Howell forms):
//   * k = 1: the rows r with p not dividing d_r are exactly the pivot rows (lows) of the column reduction
//     over F_p, so each event "p stops dividing d_r at column c" is the F_p persistence pair (r, c);
//   * any k: log_p |span mod p^k| = sum over rows of (k - min(k, v_p(d_r))).
// The sweep records every change of every d_r with the column that caused it; callers read off the above.
#pragma once
#include "bigint.hpp"
#include <vector>
#include <cstdint>
#include <numeric>
#include <algorithm>

namespace perseus {
namespace howell {

typedef unsigned __int128 u128;
typedef __int128 i128;

// arithmetic modulo m < 2^62, representatives in [0, m)
struct ModU64 {
    typedef uint64_t T;
    uint64_t m;
    explicit ModU64(uint64_t mod) : m(mod) {}
    T zero() const { return 0; }
    T modulus() const { return m; }
    T from_int(int64_t c) const { int64_t r = c % (int64_t)m; return r < 0 ? (T)(r + (int64_t)m) : (T)r; }
    T mul(T a, T b) const { return (T)((u128)a * b % m); }
    T add(T a, T b) const { T s = a + b; return s >= m ? s - m : s; }
    T sub(T a, T b) const { return a >= b ? a - b : a + m - b; }
    static bool is_zero(T a) { return a == 0; }
    static bool is_one(T a) { return a == 1; }
    static bool eq(T a, T b) { return a == b; }
    T gcd_m(T a) const { return std::gcd(a, m); }                // gcd(0, m) = m
    static T gcd(T a, T b) { return std::gcd(a, b); }
    static bool divides(T d, T a) { return a % d == 0; }
    static T div(T a, T d) { return a / d; }
    // inverse of a modulo n, gcd(a, n) = 1 (n = 1 gives 0)
    static T inv_mod(T a, T n) {
        if (n == 1) return 0;
        i128 r0 = n, r1 = a % n, s0 = 0, s1 = 1;
        while (r1) { i128 qq = r0 / r1, t = r0 - qq * r1; r0 = r1; r1 = t; t = s0 - qq * s1; s0 = s1; s1 = t; }
        i128 x = s0 % (i128)n; if (x < 0) x += n; return (T)x;
    }
    // integers a, b >= 0, not both 0: g = gcd(a, b) = s a + t b, with s, t reduced into [0, m)
    void egcd(T a, T b, T& g, T& s, T& t) const {
        i128 r0 = a, r1 = b, s0 = 1, s1 = 0, t0 = 0, t1 = 1;
        while (r1) { i128 qq = r0 / r1, x;
            x = r0 - qq * r1; r0 = r1; r1 = x; x = s0 - qq * s1; s0 = s1; s1 = x; x = t0 - qq * t1; t0 = t1; t1 = x; }
        g = (T)r0;
        i128 ss = s0 % (i128)m; if (ss < 0) ss += m; s = (T)ss;
        i128 tt = t0 % (i128)m; if (tt < 0) tt += m; t = (T)tt;
    }
    // valuation of a (a | m^k, a > 0) at a prime p, capped at cap
    static int val(T a, uint64_t p, int cap) { int v = 0; while (v < cap && a % p == 0) { a /= p; ++v; } return v; }
};

// arithmetic modulo a BigInt m >= 2, representatives in [0, m)
struct ModBig {
    typedef BigInt T;
    BigInt m;
    explicit ModBig(const BigInt& mod) : m(mod) {}
    T zero() const { return BigInt(0); }
    T modulus() const { return m; }
    static BigInt nonneg_mod(const BigInt& a, const BigInt& n) { BigInt r = a % n; if (r.neg) r = r + n; return r; }
    T from_int(int64_t c) const { BigInt b = BigInt::from_u64((uint64_t)(c < 0 ? -c : c)); if (c < 0) b = -b; return nonneg_mod(b, m); }
    T mul(const T& a, const T& b) const { return (a * b) % m; }
    T add(const T& a, const T& b) const { BigInt s = a + b; if (!(s < m)) s = s - m; return s; }
    T sub(const T& a, const T& b) const { BigInt s = a - b; if (s.neg) s = s + m; return s; }
    static bool is_zero(const T& a) { return a.is_zero(); }
    static bool is_one(const T& a) { return a == BigInt(1); }
    static bool eq(const T& a, const T& b) { return a == b; }
    T gcd_m(const T& a) const { return a.is_zero() ? m : perseus::gcd(a, m); }
    static T gcd(const T& a, const T& b) { return perseus::gcd(a, b); }
    static bool divides(const T& d, const T& a) { return (a % d).is_zero(); }
    static T div(const T& a, const T& d) { return a / d; }
    static void egcd_signed(const BigInt& a, const BigInt& b, BigInt& g, BigInt& s, BigInt& t) {
        BigInt r0 = a, r1 = b, s0(1), s1(0), t0(0), t1(1);
        while (!r1.is_zero()) { BigInt qq, rr; divmod(r0, r1, qq, rr);
            r0 = r1; r1 = rr; BigInt x = s0 - qq * s1; s0 = s1; s1 = x; x = t0 - qq * t1; t0 = t1; t1 = x; }
        g = r0; s = s0; t = t0;
    }
    static T inv_mod(const T& a, const T& n) {
        if (n == BigInt(1)) return BigInt(0);
        BigInt g, s, t; egcd_signed(nonneg_mod(a, n), n, g, s, t);
        return nonneg_mod(s, n);
    }
    void egcd(const T& a, const T& b, T& g, T& s, T& t) const {
        BigInt ss, tt; egcd_signed(a, b, g, ss, tt); s = nonneg_mod(ss, m); t = nonneg_mod(tt, m);
    }
    static int val(T a, const BigInt& p, int cap) {
        int v = 0; while (v < cap) { BigInt qq, rr; divmod(a, p, qq, rr); if (!rr.is_zero()) break; a = qq; ++v; } return v;
    }
};

// One sweep of a sparse integer matrix given column by column (rows ascending within a column).
template <class M>
struct Sweep {
    typedef typename M::T T;
    struct Ent { int32_t row; T c; };
    struct Event { int32_t col, row; T dold, dnew; };   // dold == m: the row had no pivot before

    const M& R;
    std::vector<std::vector<Ent>> piv;   // pivot vector per row (empty: none)
    std::vector<T> a, d, inv;            // pivot entry, d = gcd(a, m), inverse of a/d modulo m/d
    std::vector<Event> events;
    size_t ops = 0;                      // vector additions and scalings, for diagnostics

    Sweep(const M& ring, size_t nrows) : R(ring), piv(nrows), a(nrows), d(nrows), inv(nrows) {}

    // x + c y
    void axpy(std::vector<Ent>& x, const T& c, const std::vector<Ent>& y, std::vector<Ent>& out) const {
        out.clear(); size_t i = 0, j = 0;
        while (i < x.size() || j < y.size()) {
            if (j == y.size() || (i < x.size() && x[i].row < y[j].row)) out.push_back(x[i++]);
            else if (i == x.size() || y[j].row < x[i].row) { T s = R.mul(c, y[j].c); if (!M::is_zero(s)) out.push_back({y[j].row, s}); ++j; }
            else { T s = R.add(x[i].c, R.mul(c, y[j].c)); if (!M::is_zero(s)) out.push_back({x[i].row, s}); ++i; ++j; }
        }
    }
    std::vector<Ent> scale(const T& c, const std::vector<Ent>& y) const {
        std::vector<Ent> out; for (auto& e : y) { T s = R.mul(c, e.c); if (!M::is_zero(s)) out.push_back({e.row, s}); } return out;
    }

    void insert(std::vector<Ent> x0, int32_t col) {
        const T m = R.modulus();
        std::vector<std::vector<Ent>> queue; queue.push_back(std::move(x0));
        std::vector<Ent> x, tmp;
        while (!queue.empty()) {
            x.swap(queue.back()); queue.pop_back();
            while (!x.empty()) {
                const int32_t r = x.back().row; const T av = x.back().c;
                if (piv[r].empty()) {
                    T dd = R.gcd_m(av);
                    events.push_back({col, r, m, dd});
                    a[r] = av; d[r] = dd; inv[r] = M::inv_mod(M::div(av, dd), M::div(m, dd));
                    if (M::is_one(dd)) { piv[r].swap(x); x.clear(); break; }
                    std::vector<Ent> ann = scale(M::div(m, dd), x); ++ops;
                    piv[r].swap(x); x.swap(ann);       // the annihilator multiple has a zero on row r
                    continue;
                }
                if (M::divides(d[r], av)) {            // reduce: c a_r = av (mod m)
                    T c = R.mul(M::div(av, d[r]), inv[r]);
                    axpy(x, R.sub(R.zero(), c), piv[r], tmp); x.swap(tmp); ++ops;
                    continue;
                }
                // combine: U = s x + t w has entry g = gcd(av, a_r) on row r; both remainders lose row r
                T g, s, t; R.egcd(av, a[r], g, s, t);
                std::vector<Ent> U, sx = scale(s, x); axpy(sx, t, piv[r], U);
                std::vector<Ent> x1, w1;
                axpy(x, R.sub(R.zero(), M::div(av, g)), U, x1);
                axpy(piv[r], R.sub(R.zero(), M::div(a[r], g)), U, w1);
                ops += 4;
                T dn = R.gcd_m(g);
                events.push_back({col, r, d[r], dn});
                a[r] = g; d[r] = dn; inv[r] = M::inv_mod(M::div(g, dn), M::div(m, dn));
                if (!w1.empty()) queue.push_back(std::move(w1));
                if (!M::is_one(dn)) { std::vector<Ent> ann = scale(M::div(m, dn), U); if (!ann.empty()) queue.push_back(std::move(ann)); }
                piv[r].swap(U);
                x.swap(x1);
            }
        }
    }
};

} // namespace howell
} // namespace perseus
