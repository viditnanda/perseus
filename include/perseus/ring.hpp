// ring.hpp -- coefficient rings for the Perseus engines.
//
// Every engine is templated on a Ring with the interface
//   using T;                       element type (trivially copyable)
//   T zero(), one();  T from_int(int64_t);
//   T add(T,T), neg(T), mul(T,T), inv(T);
//   bool is_zero(T)                (zero in every lane)
//   bool partly_zero(T)            (zero in some lane but not all: the rational is
//                                   nonzero but one modulus divides it -> unlucky run)
//   bool is_pm1(T)                 (+-1 in every lane)
//   void residues(T, uint64_t* out)  standard representatives, one per lane
//   void store(T, Pivot&)          record a non-unit pivot (residues, or the exact value)
//   static constexpr int lanes;
//
// Rings provided:
//   GF2            -- F_2, one byte, xor/and.
//   Q64            -- exact rationals with int64 numerator/denominator; on overflow it sets a flag and
//                     degenerates to zero (the caller reruns with QBig).
//   QBig           -- exact rationals with arbitrary precision.
//   Mont<W,K>      -- K independent odd moduli (W = uint32_t for moduli < 2^31,
//                     W = uint64_t for moduli < 2^63), Montgomery form.
//                     K = 1 is the prime field F_p; K >= 2 carries Q through K large
//                     primes at once (multi-modular), which is how we compute over Q.
#pragma once
#include <cstdint>
#include <cstddef>
#include <vector>
#include <stdexcept>
#include <random>
#include <set>
#include <array>
#include <string>
#include <cstdio>
#include <stdexcept>
#include <cmath>
#include "types.hpp"

namespace perseus {

typedef unsigned __int128 u128;

struct GF2 {
    using T = uint8_t;
    static constexpr int lanes = 1;
    T zero() const { return 0; }
    T one() const { return 1; }
    T from_int(int64_t x) const { return (T)(x & 1); }
    T add(T a, T b) const { return a ^ b; }
    T neg(T a) const { return a; }
    T mul(T a, T b) const { return a & b; }
    T inv(T a) const { return a; }
    bool is_zero(T a) const { return a == 0; }
    bool partly_zero(T) const { return false; }
    bool is_pm1(T a) const { return a == 1; }
    void residues(T a, uint64_t* out) const { out[0] = a; }
    void store(T a, Pivot& pv) const { pv.res.assign(1, a); }
    uint64_t modulus(int) const { return 2; }
};

template <class W> struct DW;
template <> struct DW<uint32_t> { typedef uint64_t type; };
template <> struct DW<uint64_t> { typedef u128 type; };

template <class W, int K>
struct Mont {
    typedef typename DW<W>::type D;
    static constexpr int lanes = K;
    static constexpr int BITS = 8 * sizeof(W);
    struct T { W v[K]; };
    W m[K], mi[K], r2[K], one_[K];

    Mont() {}
    explicit Mont(const std::vector<uint64_t>& mods) {
        if ((int)mods.size() < K) throw std::runtime_error("Mont: not enough moduli");
        for (int l = 0; l < K; ++l) {
            uint64_t q = mods[l];
            if (!(q & 1) || q < 3) throw std::runtime_error("Mont: modulus must be odd and >= 3");
            if (BITS == 32 && q >= (1ULL << 31)) throw std::runtime_error("Mont<32>: modulus must be < 2^31");
            if (BITS == 64 && q >= (1ULL << 63)) throw std::runtime_error("Mont<64>: modulus must be < 2^63");
            m[l] = (W)q;
            W x = 1;                                     // Newton: x = q^{-1} mod 2^BITS
            for (int i = 0; i < 7; ++i) x *= (W)2 - (W)q * x;
            mi[l] = (W)(0 - x);                          // -q^{-1}
            // r2 = 2^(2*BITS) mod q
            u128 r = ((u128)1 << BITS) % q;              // R mod q  (BITS <= 64)
            r = (r * r) % q;
            r2[l] = (W)r;
            one_[l] = redc((D)r2[l], l);                 // R mod q = mont(1)
        }
    }
    inline W redc(D t, int l) const {
        W u = (W)t * mi[l];
        // t < m*2^BITS and u*m < m*2^BITS, so t + u*m < 2m*2^BITS fits in D (m < 2^(BITS-1))
        W r = (W)((t + (D)u * m[l]) >> BITS);
        return r >= m[l] ? r - m[l] : r;
    }
    inline W to_m(uint64_t x, int l) const { return redc((D)(W)(x % m[l]) * r2[l], l); }
    inline uint64_t from_m(W a, int l) const { return (uint64_t)redc((D)a, l); }

    T zero() const { T z; for (int l = 0; l < K; ++l) z.v[l] = 0; return z; }
    T one() const { T z; for (int l = 0; l < K; ++l) z.v[l] = one_[l]; return z; }
    T from_int(int64_t x) const {
        T z;
        for (int l = 0; l < K; ++l) {
            uint64_t a = x >= 0 ? (uint64_t)x % m[l] : (m[l] - (uint64_t)(-(x + 1)) % m[l] - 1) % m[l];
            z.v[l] = to_m(a, l);
        }
        return z;
    }
    T add(T a, T b) const {
        T z;
        for (int l = 0; l < K; ++l) { W s = a.v[l] + b.v[l]; z.v[l] = s >= m[l] ? s - m[l] : s; }
        return z;
    }
    T neg(T a) const { T z; for (int l = 0; l < K; ++l) z.v[l] = a.v[l] ? m[l] - a.v[l] : 0; return z; }
    T mul(T a, T b) const { T z; for (int l = 0; l < K; ++l) z.v[l] = redc((D)a.v[l] * b.v[l], l); return z; }
    W pw(W b, uint64_t e, int l) const { W r = one_[l]; while (e) { if (e & 1) r = redc((D)r * b, l); b = redc((D)b * b, l); e >>= 1; } return r; }
    T inv(T a) const { T z; for (int l = 0; l < K; ++l) z.v[l] = a.v[l] ? pw(a.v[l], (uint64_t)m[l] - 2, l) : 0; return z; }
    bool is_zero(T a) const { for (int l = 0; l < K; ++l) if (a.v[l]) return false; return true; }
    bool partly_zero(T a) const {
        bool z = false, nz = false;
        for (int l = 0; l < K; ++l) { if (a.v[l]) nz = true; else z = true; }
        return z && nz;
    }
    bool is_pm1(T a) const {
        bool p = true, n = true;
        for (int l = 0; l < K; ++l) { if (a.v[l] != one_[l]) p = false; if (a.v[l] != m[l] - one_[l]) n = false; }
        return p || n;
    }
    void residues(T a, uint64_t* out) const { for (int l = 0; l < K; ++l) out[l] = from_m(a.v[l], l); }
    void store(T a, Pivot& pv) const { pv.res.resize(K); residues(a, pv.res.data()); }
    uint64_t modulus(int l) const { return m[l]; }
};

// --- exact rationals -----------------------------------------------------------------
inline int64_t gcd_i64(int64_t a, int64_t b) { if (a < 0) a = -a; if (b < 0) b = -b; while (b) { int64_t t = a % b; a = b; b = t; } return a; }
struct Q64 {
    struct T { int64_t n, d; };                   // d > 0, gcd(n, d) = 1, |n|, d <= LIM
    static constexpr int lanes = 1;
    static constexpr int64_t LIM = (int64_t)1 << 62;
    mutable bool overflow = false;
    T zero() const { return {0, 1}; }
    T one() const { return {1, 1}; }
    T bad() const { overflow = true; return {0, 1}; }
    T chk(int64_t n, int64_t d) const { if (n > LIM || n < -LIM || d > LIM) return bad(); return {n, d}; }
    T from_int(int64_t x) const { return chk(x, 1); }
    T norm(int64_t n, int64_t d) const { if (n == 0) return {0, 1}; int64_t g = gcd_i64(n, d); return chk(n / g, d / g); }
    T add(T a, T b) const {
        if (overflow) return {0, 1};
        int64_t r;
        if (a.d == 1 && b.d == 1) { if (__builtin_add_overflow(a.n, b.n, &r)) return bad(); return chk(r, 1); }
        int64_t g = gcd_i64(a.d, b.d), da = a.d / g, db = b.d / g, x, y, d;
        if (__builtin_mul_overflow(a.n, db, &x) || __builtin_mul_overflow(b.n, da, &y) || __builtin_add_overflow(x, y, &r) ||
            __builtin_mul_overflow(a.d, db, &d)) return bad();
        return norm(r, d);
    }
    T neg(T a) const { return {-a.n, a.d}; }
    T mul(T a, T b) const {
        if (overflow) return {0, 1};
        int64_t r, d;
        if (a.d == 1 && b.d == 1) { if (__builtin_mul_overflow(a.n, b.n, &r)) return bad(); return chk(r, 1); }
        int64_t g1 = gcd_i64(a.n, b.d), g2 = gcd_i64(b.n, a.d);
        if (g1 == 0) g1 = 1; if (g2 == 0) g2 = 1;
        if (__builtin_mul_overflow(a.n / g1, b.n / g2, &r) || __builtin_mul_overflow(a.d / g2, b.d / g1, &d)) return bad();
        if (r == 0) return {0, 1};
        return chk(r, d);
    }
    T inv(T a) const { if (a.n == 0) return {0, 1}; return a.n > 0 ? T{a.d, a.n} : T{-a.d, -a.n}; }
    bool is_zero(T a) const { return overflow || a.n == 0; }
    bool partly_zero(T) const { return false; }
    bool is_pm1(T a) const { return a.d == 1 && (a.n == 1 || a.n == -1); }
    void residues(T, uint64_t*) const {}
    void store(T a, Pivot& pv) const { pv.exact = true; pv.num = BigInt(a.n); pv.den = BigInt(a.d); }
    uint64_t modulus(int) const { return 0; }
};
struct QBig {
    struct T { BigInt n, d; };                    // d > 0, gcd(n, d) = 1
    static constexpr int lanes = 1;
    mutable bool overflow = false;                // never set
    static T norm(BigInt n, BigInt d) {
        if (n.is_zero()) return {BigInt(0), BigInt(1)};
        if (d.neg) { n = -n; d = -d; }
        BigInt g = gcd(n, d);
        if (g != BigInt(1)) { n = n / g; d = d / g; }
        return {n, d};
    }
    T zero() const { return {BigInt(0), BigInt(1)}; }
    T one() const { return {BigInt(1), BigInt(1)}; }
    T from_int(int64_t x) const { return {BigInt(x), BigInt(1)}; }
    T add(const T& a, const T& b) const {
        if (a.d == BigInt(1) && b.d == BigInt(1)) return {a.n + b.n, BigInt(1)};
        return norm(a.n * b.d + b.n * a.d, a.d * b.d);
    }
    T neg(const T& a) const { return {-a.n, a.d}; }
    T mul(const T& a, const T& b) const {
        if (a.d == BigInt(1) && b.d == BigInt(1)) return {a.n * b.n, BigInt(1)};
        return norm(a.n * b.n, a.d * b.d);
    }
    T inv(const T& a) const { if (a.n.is_zero()) return zero(); return norm(a.d, a.n); }
    bool is_zero(const T& a) const { return a.n.is_zero(); }
    bool partly_zero(const T&) const { return false; }
    bool is_pm1(const T& a) const { return a.d == BigInt(1) && a.n.a.size() == 1 && a.n.a[0] == 1; }
    void residues(const T&, uint64_t*) const {}
    void store(const T& a, Pivot& pv) const { pv.exact = true; pv.num = a.n; pv.den = a.d; }
    uint64_t modulus(int) const { return 0; }
};

// --- primes ---------------------------------------------------------------------
inline uint64_t mulmod64(uint64_t a, uint64_t b, uint64_t m) { return (uint64_t)((u128)a * b % m); }
inline uint64_t powmod64(uint64_t b, uint64_t e, uint64_t m) {
    uint64_t r = 1 % m; b %= m;
    while (e) { if (e & 1) r = mulmod64(r, b, m); b = mulmod64(b, b, m); e >>= 1; }
    return r;
}
// deterministic Miller-Rabin for 64-bit integers
inline bool is_prime64(uint64_t n) {
    if (n < 2) return false;
    static const uint64_t sp[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (uint64_t p : sp) { if (n % p == 0) return n == p; }
    uint64_t d = n - 1; int s = 0; while (!(d & 1)) { d >>= 1; ++s; }
    for (uint64_t a : sp) {
        uint64_t x = powmod64(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool comp = true;
        for (int r = 1; r < s; ++r) { x = mulmod64(x, x, n); if (x == n - 1) { comp = false; break; } }
        if (comp) return false;
    }
    return true;
}
// the j-th prime below 2^31 counting down from 2^31 - 1 (offset lets a run pick a fresh set)
inline std::vector<uint64_t> big_primes(int count, int offset = 0) {
    std::vector<uint64_t> out; uint64_t x = (1ULL << 31) - 1; int skipped = 0;
    while ((int)out.size() < count) {
        if (is_prime64(x)) { if (skipped < offset) ++skipped; else out.push_back(x); }
        x -= 2;
    }
    return out;
}

// K distinct primes drawn uniformly at random from [2^(b-1), 2^b), b <= 62, excluding `avoid`
// ChaCha20 (RFC 8439) as a stream of 64-bit words, keyed by a 256-bit seed. Fast mode draws its primes
// from this stream, so the error bound rests on the seed being random and on ChaCha20 being a secure
// generator, not on the statistical quality of a non-cryptographic generator.
struct ChaCha20 {
    using result_type = uint64_t;
    uint32_t key[8], ctr = 0; uint32_t buf[16]; int pos = 16;
    static constexpr result_type min() { return 0; }
    static constexpr result_type max() { return ~0ULL; }
    explicit ChaCha20(const std::array<uint64_t, 4>& seed) {
        for (int i = 0; i < 4; ++i) { key[2 * i] = (uint32_t)seed[i]; key[2 * i + 1] = (uint32_t)(seed[i] >> 32); }
    }
    static uint32_t rotl(uint32_t x, int r) { return (x << r) | (x >> (32 - r)); }
    static void qr(uint32_t& a, uint32_t& b, uint32_t& c, uint32_t& d) {
        a += b; d ^= a; d = rotl(d, 16); c += d; b ^= c; b = rotl(b, 12);
        a += b; d ^= a; d = rotl(d, 8);  c += d; b ^= c; b = rotl(b, 7);
    }
    static void block(const uint32_t key[8], uint32_t counter, const uint32_t nonce[3], uint32_t out[16]) {
        uint32_t x[16] = {0x61707865, 0x3320646e, 0x79622d32, 0x6b206574,
                          key[0], key[1], key[2], key[3], key[4], key[5], key[6], key[7], counter, nonce[0], nonce[1], nonce[2]};
        uint32_t s0[16]; for (int i = 0; i < 16; ++i) s0[i] = x[i];
        for (int i = 0; i < 10; ++i) {
            qr(x[0], x[4], x[8], x[12]); qr(x[1], x[5], x[9], x[13]); qr(x[2], x[6], x[10], x[14]); qr(x[3], x[7], x[11], x[15]);
            qr(x[0], x[5], x[10], x[15]); qr(x[1], x[6], x[11], x[12]); qr(x[2], x[7], x[8], x[13]); qr(x[3], x[4], x[9], x[14]);
        }
        for (int i = 0; i < 16; ++i) out[i] = x[i] + s0[i];
    }
    result_type operator()() {
        if (pos >= 16) { const uint32_t nonce[3] = {0, 0, 0}; block(key, ctr++, nonce, buf); pos = 0; }
        uint64_t v = (uint64_t)buf[pos] | ((uint64_t)buf[pos + 1] << 32); pos += 2; return v;
    }
};
// 256-bit seeds, written as 64 hex digits (shorter strings are read as numbers and zero-padded)
inline std::string seed_hex(const std::array<uint64_t, 4>& s) {
    char b[65]; snprintf(b, sizeof b, "%016llx%016llx%016llx%016llx", (unsigned long long)s[3], (unsigned long long)s[2], (unsigned long long)s[1], (unsigned long long)s[0]);
    return b;
}
inline std::array<uint64_t, 4> parse_seed(const std::string& h) {
    std::string t = h; if (t.size() > 2 && t[0] == '0' && (t[1] == 'x' || t[1] == 'X')) t = t.substr(2);
    if (t.empty() || t.size() > 64) throw std::runtime_error("--seed: expected 1 to 64 hex digits");
    std::array<uint64_t, 4> s{0, 0, 0, 0};
    for (char c : t) {
        int d = (c >= '0' && c <= '9') ? c - '0' : (c >= 'a' && c <= 'f') ? c - 'a' + 10 : (c >= 'A' && c <= 'F') ? c - 'A' + 10 : -1;
        if (d < 0) throw std::runtime_error("--seed: expected 1 to 64 hex digits");
        for (int i = 3; i > 0; --i) s[i] = (s[i] << 4) | (s[i - 1] >> 60);
        s[0] = (s[0] << 4) | (uint64_t)d;
    }
    return s;
}
inline std::array<uint64_t, 4> fresh_seed() {
    std::random_device rd; std::array<uint64_t, 4> s;
    for (auto& w : s) w = ((uint64_t)rd() << 32) | (uint64_t)rd();
    return s;
}
template <class RNG>
inline std::vector<uint64_t> random_primes(int K, int b, RNG& rng, const std::set<uint64_t>& avoid = {}) {
    const uint64_t lo = 1ULL << (b - 1), mask = lo - 1;
    std::set<uint64_t> got; std::vector<uint64_t> out;
    while ((int)out.size() < K) {
        uint64_t x = (lo | (rng() & mask)) | 1;   // a uniform odd number in [2^(b-1), 2^b), then rejection: uniform over primes
        if (!is_prime64(x) || got.count(x) || avoid.count(x)) continue;
        got.insert(x); out.push_back(x);
    }
    return out;
}
// a rigorous lower bound for the number of primes in [2^(b-1), 2^b), b >= 20 (Dusart 2010:
// pi(x) >= x/ln x (1 + 1/ln x) for x >= 599, pi(x) <= x/ln x (1 + 1/ln x + 2.51/ln^2 x) for x >= 355991)
inline double primes_in_window(int b) {
    double hi = std::ldexp(1.0, b), lo = std::ldexp(1.0, b - 1), Lh = std::log(hi), Ll = std::log(lo);
    return hi / Lh * (1 + 1 / Lh) - lo / Ll * (1 + 1 / Ll + 2.51 / (Ll * Ll));
}

} // namespace perseus
