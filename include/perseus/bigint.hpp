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

namespace perseus {

struct BigInt {
    std::vector<uint32_t> a;   // magnitude, little endian, no leading zeros
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

    static int cmp_mag(const std::vector<uint32_t>& x, const std::vector<uint32_t>& y) {
        if (x.size() != y.size()) return x.size() < y.size() ? -1 : 1;
        for (size_t i = x.size(); i-- > 0;) if (x[i] != y[i]) return x[i] < y[i] ? -1 : 1;
        return 0;
    }
    static std::vector<uint32_t> add_mag(const std::vector<uint32_t>& x, const std::vector<uint32_t>& y) {
        std::vector<uint32_t> r(std::max(x.size(), y.size()) + 1, 0); uint64_t c = 0;
        for (size_t i = 0; i < r.size(); ++i) {
            uint64_t s = c + (i < x.size() ? x[i] : 0) + (i < y.size() ? y[i] : 0);
            r[i] = (uint32_t)s; c = s >> 32;
        }
        while (!r.empty() && r.back() == 0) r.pop_back();
        return r;
    }
    static std::vector<uint32_t> sub_mag(const std::vector<uint32_t>& x, const std::vector<uint32_t>& y) { // |x| >= |y|
        std::vector<uint32_t> r(x.size()); int64_t b = 0;
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
    static void divmod_mag(const std::vector<uint32_t>& u0, const std::vector<uint32_t>& v0,
                           std::vector<uint32_t>& q, std::vector<uint32_t>& r) {
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
        std::vector<uint32_t> v(n), u(u0.size() + 1);
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
    friend BigInt operator/(const BigInt& x, const BigInt& y) { BigInt q, r; divmod(x, y, q, r); return q; }
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

inline BigInt gcd(BigInt x, BigInt y) { x.neg = y.neg = false; while (!y.is_zero()) { BigInt r = x % y; x = y; y = r; } return x; }

} // namespace perseus
