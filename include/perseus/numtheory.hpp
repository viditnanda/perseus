// numtheory.hpp -- CRT, rational reconstruction and factoring on BigInt.
#pragma once
#include "bigint.hpp"
#include "ring.hpp"
#include <vector>
#include <map>
#include <random>

namespace perseus {

inline BigInt isqrt(const BigInt& n) {        // floor(sqrt(n)), n >= 0
    if (n.is_zero()) return BigInt(0);
    BigInt x = BigInt(1);
    { BigInt t = BigInt(1); size_t sh = (n.bits() + 1) / 2; // x = 2^sh >= sqrt(n)
      std::vector<uint32_t> a(sh / 32 + 1, 0); a[sh / 32] = 1u << (sh % 32); x.a = a; x.trim(); (void)t; }
    while (true) { BigInt y = (x + n / x).div_small(2); if (!(y < x)) return x; x = y; }
}

// CRT: residues r[i] mod pairwise coprime m[i] -> x mod M (0 <= x < M)
inline void crt(const std::vector<uint64_t>& r, const std::vector<uint64_t>& m, BigInt& x, BigInt& M) {
    x = BigInt(0); M = BigInt(1);
    for (size_t i = 0; i < r.size(); ++i) {
        uint64_t xm = x.mod_small(m[i]), Mm = M.mod_small(m[i]);
        uint64_t d = (r[i] % m[i] + m[i] - xm) % m[i];
        uint64_t t = mulmod64(d, powmod64(Mm, m[i] - 2, m[i]), m[i]);   // m[i] prime
        x = x + M * BigInt::from_u64(t);
        M = M * BigInt::from_u64(m[i]);
    }
}

// Wang's rational reconstruction: find n/d = x mod M with |n|, d <= sqrt(M/2). false if none.
inline bool ratrecon(const BigInt& x, const BigInt& M, BigInt& num, BigInt& den) {
    BigInt B = isqrt(M.div_small(2));
    BigInt r0 = M, r1 = x, t0(0), t1(1);
    while (B < r1) {
        BigInt q = r0 / r1;
        BigInt r2 = r0 - q * r1; r0 = r1; r1 = r2;
        BigInt t2 = t0 - q * t1; t0 = t1; t1 = t2;
    }
    if (t1.is_zero() || B < t1.abs()) return false;
    if (gcd(r1, t1) != BigInt(1)) return false;
    num = r1; den = t1;
    if (den.neg) { num = -num; den = -den; }
    return true;
}

inline BigInt powmod(BigInt b, BigInt e, const BigInt& m) {
    BigInt r(1); b = b % m;
    while (!e.is_zero()) {
        uint32_t bit; e = e.div_small(2, &bit);
        if (bit) r = (r * b) % m;
        b = (b * b) % m;
    }
    return r;
}
inline bool is_probable_prime(const BigInt& n) {
    if (n.fits_u64()) return is_prime64(n.to_u64());
    static const uint32_t sp[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};
    for (uint32_t p : sp) if (n.mod_small(p) == 0) return false;
    BigInt nm1 = n - BigInt(1), d = nm1; int s = 0;
    while (d.mod_small(2) == 0) { d = d.div_small(2); ++s; }
    for (uint32_t a : sp) {
        BigInt x = powmod(BigInt(a), d, n);
        if (x == BigInt(1) || x == nm1) continue;
        bool comp = true;
        for (int i = 1; i < s; ++i) { x = (x * x) % n; if (x == nm1) { comp = false; break; } }
        if (comp) return false;
    }
    return true;
}
inline uint64_t gcd64(uint64_t a, uint64_t b) { while (b) { uint64_t t = a % b; a = b; b = t; } return a; }
inline uint64_t rho64(uint64_t n, uint64_t c) {           // Pollard-Brent, returns a nontrivial factor or n
    if (n % 2 == 0) return 2;
    uint64_t y = 2, x, g = 1, q = 1, ys = 0; uint64_t r = 1, m = 128;
    auto f = [&](uint64_t v) { return (mulmod64(v, v, n) + c) % n; };
    do {
        x = y; for (uint64_t i = 0; i < r; ++i) y = f(y);
        uint64_t k = 0;
        do {
            ys = y;
            for (uint64_t i = 0; i < std::min(m, r - k); ++i) { y = f(y); q = mulmod64(q, x > y ? x - y : y - x, n); }
            g = gcd64(q, n); k += m;
        } while (k < r && g == 1);
        r <<= 1;
    } while (g == 1 && r < (1ULL << 40));
    if (g == n) { do { ys = f(ys); g = gcd64(x > ys ? x - ys : ys - x, n); } while (g == 1); }
    return g;
}
inline void factor64(uint64_t n, std::map<uint64_t, int>& out) {
    if (n < 2) return;
    if (is_prime64(n)) { out[n]++; return; }
    for (uint64_t c = 1;; ++c) { uint64_t d = rho64(n, c); if (d != n && d != 1) { factor64(d, out); factor64(n / d, out); return; } }
}
// factor |n|: primes fitting in 64 bits into `primes`, probable primes >= 2^64 into `bigprimes`, and
// composites resisting Pollard rho into `unresolved`
inline void factor_big(BigInt n, std::map<uint64_t, int>& primes, std::vector<BigInt>& unresolved, std::vector<BigInt>& bigprimes, uint64_t rho_budget = 1u << 16) {
    n.neg = false;
    if (n.is_zero()) return;
    static std::vector<uint32_t> small;
    if (small.empty()) { std::vector<bool> c(1 << 16, false); for (uint32_t i = 2; i < (1u << 16); ++i) if (!c[i]) { small.push_back(i); for (uint32_t j = i * i; j < (1u << 16); j += i) c[j] = true; } }
    for (uint32_t p : small) {
        if (n.fits_u64() && n.to_u64() < (uint64_t)p * p) break;
        while (n.mod_small(p) == 0) { primes[p]++; n = n.div_small(p); }
    }
    if (n == BigInt(1) || n.is_zero()) return;
    if (n.fits_u64()) { factor64(n.to_u64(), primes); return; }
    if (is_probable_prime(n)) { bigprimes.push_back(n); return; }
    // Pollard-Brent rho on BigInt, gcds batched over blocks of m steps (backtracking when a block overshoots)
    BigInt one(1);
    for (uint64_t c = 1; c < 3; ++c) {
        BigInt C = BigInt::from_u64(c), y = BigInt::from_u64(2 + c), x, ys, q(1), g(1);
        auto f = [&](const BigInt& v) { return (v * v + C) % n; };
        const uint64_t m = 128; uint64_t r = 1, steps = 0;
        do {
            x = y; for (uint64_t i = 0; i < r; ++i) y = f(y);
            uint64_t k = 0;
            do {
                ys = y;
                for (uint64_t i = 0; i < std::min(m, r - k); ++i) { y = f(y); q = (q * (x < y ? y - x : x - y)) % n; }
                g = gcd(q, n); k += m; steps += m;
            } while (k < r && g == one && steps < rho_budget);
            r <<= 1;
        } while (g == one && steps < rho_budget);
        if (g == n) { do { ys = f(ys); g = gcd(x < ys ? ys - x : x - ys, n); } while (g == one); }
        if (g != one && g != n) { factor_big(g, primes, unresolved, bigprimes, rho_budget); factor_big(n / g, primes, unresolved, bigprimes, rho_budget); return; }
    }
    unresolved.push_back(n);
}

// the field F_p for a prime p of any size (BigInt arithmetic): classification of large candidate primes
struct ModBig {
    struct T { BigInt v; };
    static constexpr int lanes = 1;
    BigInt p;
    explicit ModBig(const BigInt& q) : p(q) {}
    T red(BigInt x) const { x = x % p; if (x.neg) x = x + p; return {x}; }
    T zero() const { return {BigInt(0)}; }
    T one() const { return {BigInt(1)}; }
    T from_int(int64_t x) const { return red(BigInt(x)); }
    T add(const T& a, const T& b) const { BigInt s = a.v + b.v; if (!(s < p)) s = s - p; return {s}; }
    T neg(const T& a) const { return a.v.is_zero() ? a : T{p - a.v}; }
    T mul(const T& a, const T& b) const { return red(a.v * b.v); }
    T inv(const T& a) const {                         // extended Euclid
        BigInt r0 = p, r1 = a.v, t0(0), t1(1);
        while (!r1.is_zero()) { BigInt q = r0 / r1, r2 = r0 - q * r1; r0 = r1; r1 = r2; BigInt t2 = t0 - q * t1; t0 = t1; t1 = t2; }
        return red(t0);
    }
    bool is_zero(const T& a) const { return a.v.is_zero(); }
    bool partly_zero(const T&) const { return false; }
    bool is_pm1(const T& a) const { return a.v == BigInt(1) || a.v == p - BigInt(1); }
    void residues(const T&, uint64_t*) const {}
    void store(const T&, Pivot&) const {}
    uint64_t modulus(int) const { return 0; }
};

// Z/N for an odd N >= 3 of up to 64L - 2 bits, possibly composite (multi-limb Montgomery arithmetic).
// Used to classify all prime factors of N at once: if a pivot x has 1 < gcd(x, N) < N, the factor is
// recorded in `split` and the ring degenerates to zero so the run ends quickly; otherwise every pivot is a
// unit or zero modulo every prime factor of N, so the run is a valid F_p run for each of them, with the
// same pivots -- and they all have the same barcode.
template <int L>
struct MontN {
    struct T { uint64_t w[L]; };
    static constexpr int lanes = 1;
    uint64_t n[L], ninv; T r2, one_;
    BigInt N;
    mutable bool split_found = false; mutable BigInt split;
    static void from_big(const BigInt& x, uint64_t* w) {
        for (int i = 0; i < L; ++i) w[i] = 0;
        for (size_t i = 0; i < x.a.size() && i / 2 < (size_t)L; ++i) w[i / 2] |= (uint64_t)x.a[i] << (32 * (i & 1));
    }
    static BigInt to_big(const uint64_t* w) {
        BigInt r; r.a.resize(2 * L); for (int i = 0; i < L; ++i) { r.a[2 * i] = (uint32_t)w[i]; r.a[2 * i + 1] = (uint32_t)(w[i] >> 32); }
        r.trim(); return r;
    }
    explicit MontN(const BigInt& NN) : N(NN) {
        if ((int)NN.bits() > 64 * L - 2) throw std::runtime_error("MontN: modulus too large");
        from_big(NN, n);
        uint64_t x = 1; for (int i = 0; i < 7; ++i) x *= 2 - n[0] * x;
        ninv = 0 - x;
        // R^2 mod N with R = 2^(64L)
        BigInt R2(1); { std::vector<uint32_t> a(4 * L + 1, 0); a[4 * L] = 1; R2.a = a; R2.trim(); }
        R2 = R2 % NN; from_big(R2, r2.w);
        T o; for (int i = 0; i < L; ++i) o.w[i] = 0; o.w[0] = 1; one_ = mul(o, r2);
    }
    bool geq(const uint64_t* a, const uint64_t* b) const { for (int i = L; i-- > 0;) if (a[i] != b[i]) return a[i] > b[i]; return true; }
    void subn(uint64_t* a, const uint64_t* b) const { unsigned __int128 br = 0; for (int i = 0; i < L; ++i) { unsigned __int128 d = (unsigned __int128)a[i] - b[i] - br; a[i] = (uint64_t)d; br = (d >> 64) & 1; } }
    T mul(const T& a, const T& b) const {                      // CIOS Montgomery multiplication
        uint64_t t[L + 2] = {0};
        for (int i = 0; i < L; ++i) {
            unsigned __int128 c = 0;
            for (int j = 0; j < L; ++j) { c += (unsigned __int128)a.w[j] * b.w[i] + t[j]; t[j] = (uint64_t)c; c >>= 64; }
            c += t[L]; t[L] = (uint64_t)c; t[L + 1] = (uint64_t)(c >> 64);
            uint64_t m = t[0] * ninv;
            c = (unsigned __int128)m * n[0] + t[0]; c >>= 64;
            for (int j = 1; j < L; ++j) { c += (unsigned __int128)m * n[j] + t[j]; t[j - 1] = (uint64_t)c; c >>= 64; }
            c += t[L]; t[L - 1] = (uint64_t)c; c >>= 64; t[L] = t[L + 1] + (uint64_t)c;
        }
        T r; for (int i = 0; i < L; ++i) r.w[i] = t[i];
        if (t[L] || geq(r.w, n)) subn(r.w, n);
        return r;
    }
    T zero() const { T z; for (int i = 0; i < L; ++i) z.w[i] = 0; return z; }
    T one() const { return one_; }
    T from_int(int64_t x) const {
        BigInt b(x); b = b % N; if (b.neg) b = b + N;
        T t; from_big(b, t.w); return mul(t, r2);
    }
    T add(const T& a, const T& b) const {
        if (split_found) return zero();
        T r; unsigned __int128 c = 0;
        for (int i = 0; i < L; ++i) { c += (unsigned __int128)a.w[i] + b.w[i]; r.w[i] = (uint64_t)c; c >>= 64; }
        if (c || geq(r.w, n)) subn(r.w, n);
        return r;
    }
    T neg(const T& a) const { if (is_zero(a)) return a; T r; for (int i = 0; i < L; ++i) r.w[i] = n[i]; subn(r.w, a.w); return r; }
    BigInt value(const T& a) const { T o = zero(); o.w[0] = 1; T v = mul(a, o); return to_big(v.w); }
    T inv(const T& a) const {
        BigInt x = value(a);
        BigInt r0 = N, r1 = x, t0(0), t1(1);
        while (!r1.is_zero()) { BigInt q = r0 / r1, r2_ = r0 - q * r1; r0 = r1; r1 = r2_; BigInt t2 = t0 - q * t1; t0 = t1; t1 = t2; }
        if (r0 != BigInt(1)) { if (!split_found) { split_found = true; split = r0; } return zero(); }
        t0 = t0 % N; if (t0.neg) t0 = t0 + N;
        T t; from_big(t0, t.w); return mul(t, r2);
    }
    bool is_zero(const T& a) const { if (split_found) return true; for (int i = 0; i < L; ++i) if (a.w[i]) return false; return true; }
    bool partly_zero(const T&) const { return false; }
    bool is_pm1(const T& a) const { BigInt v = value(a); return v == BigInt(1) || v == N - BigInt(1); }
    void residues(const T&, uint64_t*) const {}
    void store(const T&, Pivot&) const {}
    uint64_t modulus(int) const { return 0; }
};

} // namespace perseus
