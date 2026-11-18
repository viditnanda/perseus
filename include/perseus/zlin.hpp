// zlin.hpp -- exact integer linear algebra for the integral outputs: Smith normal form (invariant
// factors), kernel bases and quotients of lattices. Every routine runs first in int64 with overflow
// checks and is repeated in BigInt when a check fails, so results are exact for any input size the
// matrices allow (they are dense: these routines are meant for residual complexes after reduction).
#pragma once
#include "bigint.hpp"
#include <vector>
#include <cstdint>
#include <stdexcept>
#include <algorithm>
#include <utility>

namespace perseus {
namespace zlin {

struct Overflow {};

// ---- arithmetic on int64 (checked) and BigInt, under one set of names ---------------------------------
inline int64_t zadd(int64_t a, int64_t b) { int64_t r; if (__builtin_add_overflow(a, b, &r)) throw Overflow(); return r; }
inline int64_t zsub(int64_t a, int64_t b) { int64_t r; if (__builtin_sub_overflow(a, b, &r)) throw Overflow(); return r; }
inline int64_t zmul(int64_t a, int64_t b) { int64_t r; if (__builtin_mul_overflow(a, b, &r)) throw Overflow(); return r; }
inline int64_t zquo(int64_t a, int64_t b) { if (b == -1 && a == INT64_MIN) throw Overflow(); return a / b; }   // truncated
inline bool zzero(int64_t a) { return a == 0; }
inline bool zabsless(int64_t a, int64_t b) {   // |a| < |b|
    if (a == INT64_MIN || b == INT64_MIN) throw Overflow();
    return (a < 0 ? -a : a) < (b < 0 ? -b : b);
}
inline int64_t zneg(int64_t a) { if (a == INT64_MIN) throw Overflow(); return -a; }
inline bool zneg_p(int64_t a) { return a < 0; }
inline BigInt zbig(int64_t a) { return BigInt(a); }

inline BigInt zadd(const BigInt& a, const BigInt& b) { return a + b; }
inline BigInt zsub(const BigInt& a, const BigInt& b) { return a - b; }
inline BigInt zmul(const BigInt& a, const BigInt& b) { return a * b; }
inline BigInt zquo(const BigInt& a, const BigInt& b) {   // truncated quotient
    BigInt q = a.abs() / b.abs();
    return (a.neg != b.neg) ? -q : q;
}
inline bool zzero(const BigInt& a) { return a.is_zero(); }
inline bool zabsless(const BigInt& a, const BigInt& b) { return BigInt::cmp_mag(a.a, b.a) < 0; }
inline BigInt zneg(const BigInt& a) { return -a; }
inline bool zneg_p(const BigInt& a) { return a.neg; }
inline BigInt zbig(const BigInt& a) { return a; }

template <class T> using Mat = std::vector<std::vector<T>>;   // row-major, m[r][c]

// ---- Smith normal form: rank and invariant factors d_1 | d_2 | ... (positive, nonzero) ------------------
inline std::vector<BigInt> finish_snf(std::vector<BigInt> d) { return d; }
template <class T>
std::vector<BigInt> snf_impl(Mat<T> A) {
    // the classical algorithm with a global least-modulus pivot at every round and rounded quotients, which
    // keeps the entries small (pivoting on the first nonzero entry blows up even on 8 x 9 matrices)
    size_t R = A.size(), C = R ? A[0].size() : 0;
    std::vector<BigInt> d;
    auto rquo = [](const T& a, const T& b) {   // a / b rounded to a nearest integer
        T q = zquo(a, b), r = zsub(a, zmul(q, b));
        T r2 = zadd(r, r);
        if (zabsless(b, r2)) {                 // |2r| > |b|: step q one towards a / b
            bool same = zneg_p(r) == zneg_p(b);   // r and b have the same sign
            q = same ? zadd(q, T(1)) : zsub(q, T(1));
        }
        return q;
    };
    size_t t = 0;
    while (t < R && t < C) {
        while (true) {
            size_t pr = R, pc = C;
            for (size_t i = t; i < R; ++i) for (size_t j = t; j < C; ++j)
                if (!zzero(A[i][j]) && (pr == R || zabsless(A[i][j], A[pr][pc]))) { pr = i; pc = j; }
            if (pr == R) return finish_snf(d);
            std::swap(A[t], A[pr]);
            if (pc != t) for (size_t i = 0; i < R; ++i) std::swap(A[i][t], A[i][pc]);
            bool clear = true;
            for (size_t i = t + 1; i < R; ++i) {
                if (zzero(A[i][t])) continue;
                T q = rquo(A[i][t], A[t][t]);
                for (size_t j = t; j < C; ++j) if (!zzero(A[t][j])) A[i][j] = zsub(A[i][j], zmul(q, A[t][j]));
                if (!zzero(A[i][t])) clear = false;
            }
            for (size_t j = t + 1; j < C; ++j) {
                if (zzero(A[t][j])) continue;
                T q = rquo(A[t][j], A[t][t]);
                for (size_t i = t; i < R; ++i) if (!zzero(A[i][t])) A[i][j] = zsub(A[i][j], zmul(q, A[i][t]));
                if (!zzero(A[t][j])) clear = false;
            }
            if (!clear) continue;
            // row and column t are clear; the pivot must divide every remaining entry
            size_t bad = R;
            for (size_t i = t + 1; i < R && bad == R; ++i)
                for (size_t j = t + 1; j < C; ++j) {
                    if (zzero(A[i][j])) continue;
                    if (!zzero(zsub(A[i][j], zmul(zquo(A[i][j], A[t][t]), A[t][t])))) { bad = i; break; }
                }
            if (bad == R) break;
            for (size_t j = t; j < C; ++j) A[t][j] = zadd(A[t][j], A[bad][j]);
        }
        BigInt v = zbig(A[t][t]); v.neg = false; d.push_back(v);
        ++t;
    }
    return finish_snf(d);
}

// invariant factors of an integer matrix (given as int64 entries)
inline std::vector<BigInt> invariant_factors(const Mat<int64_t>& A) {
    try { return snf_impl<int64_t>(A); }
    catch (Overflow&) {
        Mat<BigInt> B(A.size());
        for (size_t i = 0; i < A.size(); ++i) for (int64_t x : A[i]) B[i].push_back(BigInt(x));
        return snf_impl<BigInt>(B);
    }
}

// ---- column echelon form by column operations (extended gcd), with optional transformation -------------
// On return A (r x c) has its nonzero columns first, in echelon form: column j has its first nonzero entry
// (the pivot) in row piv[j], strictly increasing in j; the remaining columns are zero. If V is given (c x c,
// start with the identity), the same column operations are applied to it, so A_in * V = A_out and the
// columns of V matching zero columns of A_out form a basis of the kernel of A_in.
template <class T>
size_t col_echelon(Mat<T>& A, Mat<T>* V, std::vector<size_t>& piv) {
    size_t R = A.size(), C = R ? A[0].size() : (V ? V->size() : 0);
    piv.clear();
    size_t c0 = 0;
    auto colop = [&](size_t j, size_t k, const T& a, const T& b, const T& cc, const T& dd) {
        // (col j, col k) <- (a col j + b col k, cc col j + dd col k), a unimodular 2x2 transformation
        for (size_t i = 0; i < R; ++i) { T x = A[i][j], y = A[i][k]; A[i][j] = zadd(zmul(a, x), zmul(b, y)); A[i][k] = zadd(zmul(cc, x), zmul(dd, y)); }
        if (V) for (size_t i = 0; i < V->size(); ++i) { T x = (*V)[i][j], y = (*V)[i][k]; (*V)[i][j] = zadd(zmul(a, x), zmul(b, y)); (*V)[i][k] = zadd(zmul(cc, x), zmul(dd, y)); }
    };
    auto swapcol = [&](size_t j, size_t k) {
        if (j == k) return;
        for (size_t i = 0; i < R; ++i) std::swap(A[i][j], A[i][k]);
        if (V) for (size_t i = 0; i < V->size(); ++i) std::swap((*V)[i][j], (*V)[i][k]);
    };
    T one = T(1), zero = T(0);
    for (size_t r = 0; r < R && c0 < C; ++r) {
        // gather the gcd of row r over columns c0.. into column c0 (Euclid by column operations)
        while (true) {
            size_t best = C;
            for (size_t j = c0; j < C; ++j) if (!zzero(A[r][j]) && (best == C || zabsless(A[r][j], A[r][best]))) best = j;
            if (best == C) break;
            swapcol(c0, best);
            bool rest = false;
            for (size_t j = c0 + 1; j < C; ++j) {
                if (zzero(A[r][j])) continue;
                T q = zquo(A[r][j], A[r][c0]);
                colop(j, c0, one, zneg(q), zero, one);    // col j <- col j - q col c0 (col c0 unchanged)
                if (!zzero(A[r][j])) rest = true;
            }
            if (!rest) break;
        }
        if (zzero(A[r][c0])) continue;
        piv.push_back(r); ++c0;
    }
    return c0;   // rank
}

template <class T>
Mat<T> to_T(const Mat<int64_t>& A) { Mat<T> B(A.size()); for (size_t i = 0; i < A.size(); ++i) for (int64_t x : A[i]) B[i].push_back(T(x)); return B; }

// kernel of A (r x c, int64 entries): a Z-basis, as BigInt vectors of length c
template <class T>
std::vector<std::vector<BigInt>> kernel_impl(const Mat<int64_t>& A0, size_t c) {
    Mat<T> A = to_T<T>(A0);
    Mat<T> V(c, std::vector<T>(c, T(0))); for (size_t i = 0; i < c; ++i) V[i][i] = T(1);
    if (A.empty()) A.assign(0, std::vector<T>(c));
    std::vector<size_t> piv;
    size_t rank = A.empty() ? 0 : col_echelon<T>(A, &V, piv);
    std::vector<std::vector<BigInt>> out;
    for (size_t j = rank; j < c; ++j) { std::vector<BigInt> v(c); for (size_t i = 0; i < c; ++i) v[i] = zbig(V[i][j]); out.push_back(v); }
    return out;
}
inline std::vector<std::vector<BigInt>> kernel_basis(const Mat<int64_t>& A, size_t c) {
    try { return kernel_impl<int64_t>(A, c); } catch (Overflow&) { return kernel_impl<BigInt>(A, c); }
}

// S / B for lattices B <= S <= Z^n, given by generators (columns, BigInt). Returns the free rank and the
// invariant factors > 1 of the quotient. Throws if B is not contained in S.
struct GroupStructure { size_t free_rank = 0; std::vector<BigInt> torsion; };
inline GroupStructure quotient(const std::vector<std::vector<BigInt>>& S, const std::vector<std::vector<BigInt>>& Bg, size_t n) {
    // basis of S: column echelon form of the generators
    Mat<BigInt> A(n, std::vector<BigInt>(S.size()));
    for (size_t j = 0; j < S.size(); ++j) for (size_t i = 0; i < n; ++i) A[i][j] = S[j][i];
    std::vector<size_t> piv;
    size_t m = S.empty() ? 0 : col_echelon<BigInt>(A, nullptr, piv);
    // coordinates of each generator of B in that basis (forward substitution on the echelon basis)
    Mat<BigInt> Y(m, std::vector<BigInt>(Bg.size()));
    for (size_t g = 0; g < Bg.size(); ++g) {
        std::vector<BigInt> x = Bg[g];
        for (size_t j = 0; j < m; ++j) {
            size_t r = piv[j];
            if (x[r].is_zero()) continue;
            BigInt q, rem; divmod(x[r], A[r][j], q, rem);
            if (!rem.is_zero()) throw std::runtime_error("zlin::quotient: B is not contained in S");
            Y[j][g] = q;
            for (size_t i = r; i < n; ++i) if (!A[i][j].is_zero()) x[i] = x[i] - q * A[i][j];
        }
        for (auto& e : x) if (!e.is_zero()) throw std::runtime_error("zlin::quotient: B is not contained in S");
    }
    GroupStructure G;
    std::vector<BigInt> d = m && !Bg.empty() ? snf_impl<BigInt>(Y) : std::vector<BigInt>();
    G.free_rank = m - d.size();
    for (auto& x : d) if (!(x == BigInt(1))) G.torsion.push_back(x);
    return G;
}

} // namespace zlin
} // namespace perseus
