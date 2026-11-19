// t_exact.cpp -- the exact arithmetic behind proved mode: BigInt (inline limbs, 128-bit binary gcd, one-limb
// division) against identities, and the hybrid rational ring QHyb against QBig (plain big rationals) on
// random operation sequences whose numbers cross 62, 128 and 256 bits.
#include "perseus/ring.hpp"
#include <cstdio>
#include <random>
#include <string>

using namespace perseus;

static int fails = 0, checks = 0;
static void check(bool ok, const std::string& what) { ++checks; if (!ok) { ++fails; if (fails < 20) printf("FAIL %s\n", what.c_str()); } }

static BigInt rnd_big(std::mt19937_64& g, int maxbits) {
    int bits = 1 + (int)(g() % (uint64_t)maxbits);
    BigInt x(0);
    for (int b = 0; b < bits; b += 32) { x = x * BigInt((int64_t)1 << 32) + BigInt((int64_t)(g() & 0xffffffffULL)); }
    if (g() % 2) x = -x;
    return x;
}

static bool same(const QHyb::T& a, const QBig::T& b) { return QHyb::bn(a) == b.n && QHyb::bd(a) == b.d; }

int main(int argc, char** argv) {
    int trials = argc > 1 ? atoi(argv[1]) : 2000;
    std::mt19937_64 g(argc > 2 ? strtoull(argv[2], nullptr, 10) : 11);
    // BigInt identities
    for (int t = 0; t < trials; ++t) {
        int mb = (t % 4 == 0) ? 600 : (t % 4 == 1 ? 300 : (t % 4 == 2 ? 128 : 64));
        BigInt a = rnd_big(g, mb), b = rnd_big(g, mb);
        if (b.is_zero()) continue;
        BigInt q, r; divmod(a, b, q, r);
        check(q * b + r == a, "a = q b + r");
        check(cmp(r.abs(), b.abs()) < 0, "|r| < |b|");
        check(a * b / b == a, "(a b) / b = a");
        BigInt d = gcd(a, b);
        check(!d.is_zero() && (a % d).is_zero() && (b % d).is_zero(), "gcd divides both");
        check(gcd(a / d, b / d) == BigInt(1), "cofactors coprime");
        BigInt s = BigInt((int64_t)(1 + g() % 4000000000ULL));
        check((a * s) / s == a, "one-limb division");
        check(gcd(a * s, s) == s, "gcd with a one-word number");
    }
    // QHyb against QBig
    QHyb H; QBig B;
    for (int t = 0; t < trials / 4; ++t) {
        QHyb::T x = H.from_int((int64_t)(g() % 2000) - 1000); QBig::T y = B.from_int(QHyb::to_i64(QHyb::bn(x)));
        for (int s = 0; s < 60; ++s) {
            int64_t p = (int64_t)(g() % 2000000) - 1000000, q = 1 + (int64_t)(g() % 1000000);
            if (s % 7 == 3) { p = (int64_t)(g() >> 2) * (g() % 2 ? 1 : -1); q = 1 + (int64_t)(g() >> 3); }   // near 2^62
            QHyb::T u = H.mul(H.from_int(p), H.inv(H.from_int(q))); QBig::T v = B.mul(B.from_int(p), B.inv(B.from_int(q)));
            switch (g() % 4) {
                case 0: x = H.add(x, u); y = B.add(y, v); break;
                case 1: x = H.mul(x, u); y = B.mul(y, v); break;
                case 2: x = H.add(x, H.neg(u)); y = B.add(y, B.neg(v)); break;
                default: if (!H.is_zero(x)) { x = H.inv(x); y = B.inv(y); } break;
            }
            check(same(x, y), "QHyb == QBig after an operation");
            check(H.is_zero(x) == B.is_zero(y) && H.is_pm1(x) == B.is_pm1(y), "QHyb zero and unit tests");
            if (QHyb::bn(x).bits() > 900) { x = H.from_int(1); y = B.from_int(1); }
        }
    }
    printf("%d checks, %d failures\n", checks, fails);
    return fails != 0;
}
