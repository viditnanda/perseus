// bigint.hpp -- a small arbitrary-precision signed integer, enough for multi-modular
// reconstruction (CRT, Wang's rational reconstruction) and for factoring pivot values.
// Base 2^32 limbs, schoolbook multiplication, Knuth algorithm D for division.
#pragma once
#include <cstdint>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <cstdio>
#include <cstring>

namespace perseus {

// Limbs: a vector of 32-bit limbs that keeps up to 8 of them (256 bits) inline, so the arithmetic on the
// numbers a reduction meets (pivot values of a few hundred bits) does not touch the heap.
class Limbs {
    static constexpr size_t N = 8;
    uint32_t buf_[N];
    uint32_t* p_ = buf_;
    size_t n_ = 0, cap_ = N;
    void grow(size_t c) {
        if (c <= cap_) return;
        size_t nc = std::max(c, 2 * cap_);
        uint32_t* q = new uint32_t[nc];
        if (n_) std::memcpy(q, p_, n_ * sizeof(uint32_t));
        if (p_ != buf_) delete[] p_;
        p_ = q; cap_ = nc;
    }
public:
    typedef uint32_t value_type;
    Limbs() {}
    explicit Limbs(size_t n, uint32_t v = 0) { assign(n, v); }
    Limbs(const Limbs& o) { grow(o.n_); std::copy(o.p_, o.p_ + o.n_, p_); n_ = o.n_; }
    Limbs(Limbs&& o) noexcept {
        if (o.p_ == o.buf_) { std::copy(o.buf_, o.buf_ + o.n_, buf_); n_ = o.n_; }
        else { p_ = o.p_; n_ = o.n_; cap_ = o.cap_; o.p_ = o.buf_; o.cap_ = N; }
        o.n_ = 0;
    }
    Limbs& operator=(const Limbs& o) { if (this != &o) { n_ = 0; grow(o.n_); std::copy(o.p_, o.p_ + o.n_, p_); n_ = o.n_; } return *this; }
    Limbs& operator=(Limbs&& o) noexcept {
        if (this == &o) return *this;
        if (o.p_ == o.buf_) { n_ = 0; grow(o.n_); std::copy(o.buf_, o.buf_ + o.n_, p_); n_ = o.n_; o.n_ = 0; return *this; }
        if (p_ != buf_) delete[] p_;
        p_ = o.p_; n_ = o.n_; cap_ = o.cap_; o.p_ = o.buf_; o.cap_ = N; o.n_ = 0;
        return *this;
    }
    ~Limbs() { if (p_ != buf_) delete[] p_; }
    size_t size() const { return n_; }
    bool empty() const { return n_ == 0; }
    uint32_t& operator[](size_t i) { return p_[i]; }
    const uint32_t& operator[](size_t i) const { return p_[i]; }
    uint32_t& back() { return p_[n_ - 1]; }
    const uint32_t& back() const { return p_[n_ - 1]; }
    void pop_back() { --n_; }
    void push_back(uint32_t v) { grow(n_ + 1); p_[n_++] = v; }
    void resize(size_t n, uint32_t v = 0) { grow(n); for (size_t i = n_; i < n; ++i) p_[i] = v; n_ = n; }
    void assign(size_t n, uint32_t v) { n_ = 0; resize(n, v); }
    void clear() { n_ = 0; }
    uint32_t* begin() { return p_; }
    uint32_t* end() { return p_ + n_; }
    const uint32_t* begin() const { return p_; }
    const uint32_t* end() const { return p_ + n_; }
};

struct BigInt {
    Limbs a;                   // magnitude, little endian, no leading zeros
    bool neg = false;

    BigInt() {}
    BigInt(int64_t x) { neg = x < 0; uint64_t u = neg ? (uint64_t)0 - (uint64_t)x : (uint64_t)x; set_u64(u); }
    static BigInt from_u64(uint64_t u) { BigInt r; r.set_u64(u); return r; }
    void set_u64(uint64_t u) { a.clear(); while (u) { a.push_back((uint32_t)u); u >>= 32; } }
    bool is_zero() const { return a.empty(); }
    void trim() { while (!a.empty() && a.back() == 0) a.pop_back(); if (a.empty()) neg = false; }
    size_t bits() const { if (a.empty()) return 0; return 32 * (a.size() - 1) + (32 - __builtin_clz(a.back())); }
    bool fits_u64() const { return a.size() <= 2; }
    uint64_t to_u64() const { uint64_t r = 0; for (size_t i = a.size(); i-- > 0;) r = (r << 32) | a[i]; return r; }

    static int cmp_mag(const Limbs& x, const Limbs& y) {
        if (x.size() != y.size()) return x.size() < y.size() ? -1 : 1;
        for (size_t i = x.size(); i-- > 0;) if (x[i] != y[i]) return x[i] < y[i] ? -1 : 1;
        return 0;
    }
    static Limbs add_mag(const Limbs& x, const Limbs& y) {
        Limbs r(std::max(x.size(), y.size()) + 1, 0u); uint64_t c = 0;
        for (size_t i = 0; i < r.size(); ++i) {
            uint64_t s = c + (i < x.size() ? x[i] : 0) + (i < y.size() ? y[i] : 0);
            r[i] = (uint32_t)s; c = s >> 32;
        }
        while (!r.empty() && r.back() == 0) r.pop_back();
        return r;
    }
    static Limbs sub_mag(const Limbs& x, const Limbs& y) { // |x| >= |y|
        Limbs r(x.size()); int64_t b = 0;
        for (size_t i = 0; i < x.size(); ++i) {
            int64_t s = (int64_t)x[i] - b - (i < y.size() ? (int64_t)y[i] : 0);
            b = s < 0; if (s < 0) s += ((int64_t)1 << 32);
            r[i] = (uint32_t)s;
        }
        while (!r.empty() && r.back() == 0) r.pop_back();
        return r;
    }
    friend int cmp(const BigInt& x, const BigInt& y) {
        if (x.neg != y.neg) return x.neg ? -1 : 1;
        int c = cmp_mag(x.a, y.a); return x.neg ? -c : c;
    }
    friend bool operator<(const BigInt& x, const BigInt& y) { return cmp(x, y) < 0; }
    friend bool operator==(const BigInt& x, const BigInt& y) { return cmp(x, y) == 0; }
    friend bool operator!=(const BigInt& x, const BigInt& y) { return cmp(x, y) != 0; }
    BigInt operator-() const { BigInt r = *this; if (!r.is_zero()) r.neg = !r.neg; return r; }
    friend BigInt operator+(const BigInt& x, const BigInt& y) {
        BigInt r;
        if (x.neg == y.neg) { r.a = add_mag(x.a, y.a); r.neg = x.neg; }
        else if (cmp_mag(x.a, y.a) >= 0) { r.a = sub_mag(x.a, y.a); r.neg = x.neg; }
        else { r.a = sub_mag(y.a, x.a); r.neg = y.neg; }
        r.trim(); return r;
    }
    friend BigInt operator-(const BigInt& x, const BigInt& y) { return x + (-y); }
    friend BigInt operator*(const BigInt& x, const BigInt& y) {
        BigInt r; if (x.is_zero() || y.is_zero()) return r;
        r.a.assign(x.a.size() + y.a.size(), 0);
        for (size_t i = 0; i < x.a.size(); ++i) {
            uint64_t c = 0;
            for (size_t j = 0; j < y.a.size(); ++j) {
                uint64_t t = (uint64_t)x.a[i] * y.a[j] + r.a[i + j] + c;
                r.a[i + j] = (uint32_t)t; c = t >> 32;
            }
            r.a[i + y.a.size()] += (uint32_t)c;
        }
        r.neg = x.neg != y.neg; r.trim(); return r;
    }
    // remainder of |x| by a small positive m
    uint64_t mod_small(uint64_t m) const {
        u128_t r = 0;
        for (size_t i = a.size(); i-- > 0;) r = ((r << 32) | a[i]) % m;
        return (uint64_t)r;
    }
    BigInt div_small(uint32_t m, uint32_t* rem = nullptr) const {
        BigInt q; q.a.resize(a.size()); uint64_t r = 0;
        for (size_t i = a.size(); i-- > 0;) { uint64_t cur = (r << 32) | a[i]; q.a[i] = (uint32_t)(cur / m); r = cur % m; }
        q.neg = neg; q.trim(); if (rem) *rem = (uint32_t)r; return q;
    }
    typedef unsigned __int128 u128_t;
    // truncated division of magnitudes (Knuth D); returns quotient, sets remainder
    static void divmod_mag(const Limbs& u0, const Limbs& v0,
                           Limbs& q, Limbs& r) {
        if (v0.empty()) throw std::runtime_error("BigInt: division by zero");
        if (cmp_mag(u0, v0) < 0) { q.clear(); r = u0; return; }
        if (v0.size() == 1) {
            q.assign(u0.size(), 0); uint64_t rem = 0;
            for (size_t i = u0.size(); i-- > 0;) { uint64_t cur = (rem << 32) | u0[i]; q[i] = (uint32_t)(cur / v0[0]); rem = cur % v0[0]; }
            while (!q.empty() && q.back() == 0) q.pop_back();
            r.clear(); if (rem) r.push_back((uint32_t)rem);
            return;
        }
        int s = __builtin_clz(v0.back());
        size_t n = v0.size(), m = u0.size() - n;
        Limbs v(n), u(u0.size() + 1);
        for (size_t i = n; i-- > 0;) v[i] = (v0[i] << s) | (s && i ? (uint32_t)((uint64_t)v0[i - 1] >> (32 - s)) : 0);
        u[u0.size()] = s ? (uint32_t)((uint64_t)u0.back() >> (32 - s)) : 0;
        for (size_t i = u0.size(); i-- > 0;) u[i] = (u0[i] << s) | (s && i ? (uint32_t)((uint64_t)u0[i - 1] >> (32 - s)) : 0);
        q.assign(m + 1, 0);
        const uint64_t B = (uint64_t)1 << 32;
        for (size_t j = m + 1; j-- > 0;) {
            uint64_t num = ((uint64_t)u[j + n] << 32) | u[j + n - 1];
            uint64_t qhat = num / v[n - 1], rhat = num % v[n - 1];
            while (qhat >= B || qhat * v[n - 2] > ((rhat << 32) | u[j + n - 2])) { --qhat; rhat += v[n - 1]; if (rhat >= B) break; }
            int64_t borrow = 0; uint64_t carry = 0;
            for (size_t i = 0; i < n; ++i) {
                uint64_t p = qhat * v[i] + carry; carry = p >> 32;
                int64_t t = (int64_t)u[i + j] - borrow - (int64_t)(uint32_t)p;
                borrow = t < 0; u[i + j] = (uint32_t)(t + (borrow ? (int64_t)B : 0));
            }
            int64_t t = (int64_t)u[j + n] - borrow - (int64_t)carry;
            borrow = t < 0; u[j + n] = (uint32_t)(t + (borrow ? (int64_t)B : 0));
            if (borrow) {   // add back
                --qhat; uint64_t c = 0;
                for (size_t i = 0; i < n; ++i) { uint64_t s2 = (uint64_t)u[i + j] + v[i] + c; u[i + j] = (uint32_t)s2; c = s2 >> 32; }
                u[j + n] += (uint32_t)c;
            }
            q[j] = (uint32_t)qhat;
        }
        while (!q.empty() && q.back() == 0) q.pop_back();
        r.assign(n, 0);
        for (size_t i = 0; i < n; ++i) r[i] = (u[i] >> s) | (s ? (uint32_t)((uint64_t)u[i + 1] << (32 - s)) : 0);
        while (!r.empty() && r.back() == 0) r.pop_back();
    }
    // floor division / non-negative remainder for nonnegative operands; truncation otherwise
    friend void divmod(const BigInt& x, const BigInt& y, BigInt& q, BigInt& r) {
        divmod_mag(x.a, y.a, q.a, r.a);
        q.neg = x.neg != y.neg; r.neg = x.neg; q.trim(); r.trim();
    }
    friend BigInt operator/(const BigInt& x, const BigInt& y) {
        if (y.a.size() == 1) { BigInt q = x.div_small(y.a[0]); if (y.neg && !q.is_zero()) q.neg = !q.neg; return q; }
        BigInt q, r; divmod(x, y, q, r); return q;
    }
    friend BigInt operator%(const BigInt& x, const BigInt& y) { BigInt q, r; divmod(x, y, q, r); return r; }
    BigInt abs() const { BigInt r = *this; r.neg = false; return r; }

    std::string str() const {
        if (a.empty()) return "0";
        std::string s; BigInt t = abs();
        while (!t.is_zero()) { uint32_t rem; t = t.div_small(1000000000u, &rem);
            char buf[16]; snprintf(buf, sizeof buf, t.is_zero() ? "%u" : "%09u", rem);
            std::string chunk(buf); std::reverse(chunk.begin(), chunk.end()); s += chunk; }
        if (neg) s += '-';
        std::reverse(s.begin(), s.end()); return s;
    }
};

// gcd: Euclid's remainders while the numbers are longer than 128 bits, then binary gcd on two unsigned
// 128-bit words (the pivot values and the entries of a reduction are mostly in that range)
namespace detail_gcd {
typedef unsigned __int128 u128;
inline int ctz128(u128 x) { uint64_t lo = (uint64_t)x; return lo ? __builtin_ctzll(lo) : 64 + __builtin_ctzll((uint64_t)(x >> 64)); }
inline u128 to128(const BigInt& x) { u128 r = 0; for (size_t i = x.a.size(); i-- > 0;) r = (r << 32) | x.a[i]; return r; }
inline BigInt from128(u128 v) { BigInt r; while (v) { r.a.push_back((uint32_t)v); v >>= 32; } return r; }
inline u128 binary(u128 a, u128 b) {
    if (a == 0) return b;
    if (b == 0) return a;
    int s = ctz128(a | b); a >>= ctz128(a);
    do { b >>= ctz128(b); if (a > b) { u128 t = a; a = b; b = t; } b -= a; } while (b);
    return a << s;
}
}
inline BigInt gcd(BigInt x, BigInt y) {
    x.neg = y.neg = false;
    if (x.is_zero()) return y;
    if (y.is_zero()) return x;
    if (x.a.size() <= 2 || y.a.size() <= 2) {          // one word: reduce the other modulo it, then binary gcd
        if (x.a.size() > 2) std::swap(x, y);            // x fits in 64 bits
        uint64_t s = x.to_u64(), r = y.a.size() <= 2 ? y.to_u64() : y.mod_small(s);
        return detail_gcd::from128(detail_gcd::binary(s, r));
    }
    while (!y.is_zero()) {
        if (x.a.size() <= 4 && y.a.size() <= 4) return detail_gcd::from128(detail_gcd::binary(detail_gcd::to128(x), detail_gcd::to128(y)));
        BigInt r = x % y; x = y; y = r;
    }
    return x;
}

} // namespace perseus
