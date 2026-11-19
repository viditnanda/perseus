// laurent.hpp -- filtered chain complexes with coefficients in the Laurent ring Z[t, t^-1] (twisted / Novikov
// coefficients: the chains of an infinite cyclic cover, Milnor fibres, Alexander modules).
//
// Z[t, t^-1] is not a principal ideal domain, so its homology has no barcode or Smith form. What is computed is
// what the principal ideal domains around it see (docs/TECHNIQUES.md, "Laurent coefficients"):
//
//   generic fibre    Q(t)               the barcode (column reduction over the field Q(t))
//   PID              Q[t, t^-1]         the homology of the whole complex: free rank and invariant factors
//                                       (Alexander polynomials), by a Smith form over Q[t]
//   special fibres   Q[t]/(f)           the BAD IRREDUCIBLES f, whose residue field gives another barcode
//   integral content F_p(t), F_p[t,t^-1] the BAD PRIMES p, whose barcode over F_p(t), or whose homology over
//                                       F_p[t, t^-1], is not the reduction of the rational answer
//
// Everything works on integer polynomials (fraction-free): a nonzero integer is a unit of every ring above
// except F_p, and t is a unit of every one of them, so rows and columns are scaled by integers and powers of t
// freely and the scalars that F_p could see are recorded. The certificates mirror docs/CERTIFICATE.md: the
// pivot values of the reduction over Q(t) (ratios of minors, independent of the reduction) are units in the
// discrete valuation ring of every irreducible f and every prime p that divides none of their numerators or
// denominators, so only those are candidates; each candidate is settled by a rerun over its residue field,
// and an unfactored polynomial g is settled for all its irreducible factors at once by a run over Q[t]/(g) as
// if it were a field (a non-invertible pivot splits g).
#pragma once
#include "bigint.hpp"
#include "numtheory.hpp"
#include <vector>
#include <map>
#include <set>
#include <string>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <cmath>
#include <cctype>

namespace perseus {
namespace laurent {

typedef std::vector<BigInt> Poly;     // c[0] + c[1] t + ... ; no trailing zeros; empty = 0

// ---- coefficient domain: Z (p = 0, standing for Q by fraction-free arithmetic) or F_p ---------------------
struct Dom {
    uint64_t p = 0;
    Dom() {}
    explicit Dom(uint64_t p_) : p(p_) {}
    BigInt red(const BigInt& x) const {
        if (!p) return x;
        uint64_t r = x.mod_small(p); if (x.neg && r) r = p - r;
        return BigInt::from_u64(r);
    }
    void trim(Poly& a) const { if (p) for (auto& x : a) x = red(x); while (!a.empty() && a.back().is_zero()) a.pop_back(); }
    uint64_t inv(uint64_t a) const {       // a^-1 mod p
        unsigned __int128 r = 1, b = a % p; uint64_t e = p - 2;
        while (e) { if (e & 1) r = r * b % p; b = b * b % p; e >>= 1; }
        return (uint64_t)r;
    }
};

inline int deg(const Poly& a) { return (int)a.size() - 1; }
inline const BigInt& lc(const Poly& a) { return a.back(); }
inline Poly padd(const Poly& a, const Poly& b, const Dom& D) {
    Poly r(std::max(a.size(), b.size()));
    for (size_t i = 0; i < r.size(); ++i) r[i] = (i < a.size() ? a[i] : BigInt()) + (i < b.size() ? b[i] : BigInt());
    D.trim(r); return r;
}
inline Poly pneg(Poly a) { for (auto& x : a) x = -x; return a; }
inline Poly psub(const Poly& a, const Poly& b, const Dom& D) { return padd(a, pneg(b), D); }
inline Poly pmul(const Poly& a, const Poly& b, const Dom& D) {
    if (a.empty() || b.empty()) return Poly();
    Poly r(a.size() + b.size() - 1);
    for (size_t i = 0; i < a.size(); ++i) { if (a[i].is_zero()) continue; for (size_t j = 0; j < b.size(); ++j) if (!b[j].is_zero()) r[i + j] = r[i + j] + a[i] * b[j]; }
    D.trim(r); return r;
}
inline Poly pscale(const BigInt& c, Poly a, const Dom& D) { for (auto& x : a) x = c * x; D.trim(a); return a; }
inline Poly pshift(const Poly& a, int k) { if (a.empty()) return a; Poly r(k, BigInt()); r.insert(r.end(), a.begin(), a.end()); return r; }   // t^k a, k >= 0
inline int pval(const Poly& a) { int v = 0; while (v < (int)a.size() && a[v].is_zero()) ++v; return v; }
inline Poly strip_t(const Poly& a) { int v = pval(a); return Poly(a.begin() + v, a.end()); }
inline Poly pconst(const BigInt& c, const Dom& D) { Poly r{c}; D.trim(r); return r; }

// content: over Z the gcd of the coefficients, with the sign of the leading one; over F_p the leading coefficient
inline BigInt content(const Poly& a, const Dom& D) {
    if (a.empty()) return BigInt(1);
    if (D.p) return lc(a);
    BigInt g; for (auto& x : a) { g = gcd(g, x); if (g == BigInt(1)) break; }
    if (lc(a).neg) g = -g;
    return g;
}
// a / c for a scalar c dividing every coefficient (Z) or c != 0 (F_p)
inline Poly pdiv_scalar(Poly a, const BigInt& c, const Dom& D) {
    if (D.p) { BigInt ci = BigInt::from_u64(D.inv(D.red(c).to_u64())); return pscale(ci, a, D); }
    for (auto& x : a) { BigInt q, r; divmod(x, c, q, r); if (!r.is_zero()) throw std::logic_error("laurent: inexact scalar division"); x = q; }
    return a;
}
inline Poly prim(const Poly& a, const Dom& D) { if (a.empty()) return a; return pdiv_scalar(a, content(a, D), D); }   // lc > 0 (Z), monic (F_p)

// s * b = q * a + r with deg r < deg a, s = lc(a)^e (pseudo-division; exact division over F_p up to the scalar)
inline void pquorem(const Poly& b, const Poly& a, const Dom& D, BigInt& s, Poly& q, Poly& r) {
    if (a.empty()) throw std::logic_error("laurent: division by zero polynomial");
    s = BigInt(1); q.clear(); r = b;
    const BigInt& la = lc(a);
    while (!r.empty() && deg(r) >= deg(a)) {
        int k = deg(r) - deg(a); BigInt c = lc(r);
        r = psub(pscale(la, r, D), pmul(pshift(Poly{c}, k), a, D), D);
        q = padd(pscale(la, q, D), pshift(pconst(c, D), k), D);
        s = D.p ? D.red(s * la) : s * la;
    }
}
inline Poly prem(const Poly& b, const Poly& a, const Dom& D) { BigInt s; Poly q, r; pquorem(b, a, D, s, q, r); return r; }
// lc(a)^E * b mod a, for E >= deg b - deg a + 1 (exactly E pseudo-division steps, so every b gets the same scalar)
inline Poly prem_fixed(const Poly& b, const Poly& a, int E, const Dom& D) {
    Poly r = b; const BigInt& la = lc(a);
    for (int i = E - 1; i >= 0; --i) {
        size_t at = (size_t)(deg(a) + i);
        BigInt c = at < r.size() ? r[at] : BigInt();
        r = pscale(la, r, D);
        if (!c.is_zero()) r = psub(r, pmul(pshift(Poly{c}, i), a, D), D);
    }
    if (!r.empty() && deg(r) >= deg(a)) throw std::logic_error("laurent: prem_fixed exponent too small");
    return r;
}

// gcd, normalised (Z: content gcd times the primitive gcd, lc > 0; F_p: monic). gcd(0, 0) = 0.
inline Poly pgcd(const Poly& a0, const Poly& b0, const Dom& D) {
    if (a0.empty()) return prim(b0, D);
    if (b0.empty()) return prim(a0, D);
    BigInt c = D.p ? BigInt(1) : gcd(content(a0, D), content(b0, D));
    Poly a = prim(a0, D), b = prim(b0, D);
    if (deg(a) < deg(b)) std::swap(a, b);
    while (!b.empty()) { Poly r = prem(a, b, D); a = b; b = prim(r, D); }
    return D.p ? a : pscale(c, a, D);
}
// a / b, b dividing a exactly (over Z the quotient is integral when b is primitive)
inline Poly pexact(const Poly& a, const Poly& b, const Dom& D) {
    if (b.empty()) throw std::logic_error("laurent: exact division by 0");
    Poly r = a, q(std::max(0, deg(a) - deg(b) + 1));
    BigInt li = D.p ? BigInt::from_u64(D.inv(D.red(lc(b)).to_u64())) : BigInt(1);
    while (!r.empty() && deg(r) >= deg(b)) {
        int k = deg(r) - deg(b); BigInt c;
        if (D.p) c = D.red(lc(r) * li);
        else { BigInt rr; divmod(lc(r), lc(b), c, rr); if (!rr.is_zero()) throw std::logic_error("laurent: inexact polynomial division"); }
        q[k] = c;
        r = psub(r, pmul(pshift(Poly{c}, k), b, D), D);
    }
    if (!r.empty()) throw std::logic_error("laurent: inexact polynomial division (remainder)");
    D.trim(q); return q;
}
inline Poly pderiv(const Poly& a, const Dom& D) { Poly r; for (size_t i = 1; i < a.size(); ++i) r.push_back(BigInt((int64_t)i) * a[i]); D.trim(r); return r; }
// squarefree part over Q (characteristic 0), primitive, t stripped
inline Poly sqfree(const Poly& a0) {
    Dom Z; Poly a = prim(strip_t(a0), Z);
    if (deg(a) < 1) return a;
    Poly g = pgcd(a, pderiv(a, Z), Z);
    return deg(g) < 1 ? a : prim(pexact(a, prim(g, Z), Z), Z);
}
// normal form up to units of K[t, t^-1] (K = Q or F_p): t stripped, primitive (lc > 0) or monic
inline Poly unit_normal(const Poly& a, const Dom& D) { return prim(strip_t(a), D); }
inline Poly reduce_mod_p(const Poly& a, uint64_t p) { Dom D(p); Poly r = a; D.trim(r); return r; }

inline std::string to_string(const Poly& a, int shift = 0) {
    if (a.empty()) return "0";
    std::string s;
    for (int i = deg(a); i >= 0; --i) {
        if (a[i].is_zero()) continue;
        BigInt c = a[i]; bool neg = c.neg; c.neg = false; int e = i + shift;
        s += s.empty() ? (neg ? "-" : "") : (neg ? " - " : " + ");
        bool one = c == BigInt(1);
        if (e == 0) s += c.str();
        else {
            if (!one) s += c.str();
            s += "t"; if (e != 1) s += "^" + std::to_string(e);
        }
    }
    return s;
}

// ---- Laurent polynomials (input/output only) -------------------------------------------------------------
struct LPoly {                                 // t^v * c(t), c(0) != 0 unless zero
    int v = 0; Poly c;
    bool zero() const { return c.empty(); }
    void norm() { int k = pval(c); if (c.empty()) { v = 0; return; } c.erase(c.begin(), c.begin() + k); v += k; }
    std::string str() const { return to_string(c, v); }
};
inline LPoly lmul(const LPoly& a, const LPoly& b) { LPoly r; Dom Z; r.c = pmul(a.c, b.c, Z); r.v = a.v + b.v; r.norm(); return r; }
inline LPoly ladd(LPoly a, LPoly b) {
    if (a.zero()) return b; if (b.zero()) return a;
    Dom Z; int v = std::min(a.v, b.v);
    LPoly r; r.v = v; r.c = padd(pshift(a.c, a.v - v), pshift(b.c, b.v - v), Z); r.norm(); return r;
}
inline LPoly lmonomial(int64_t c, int e) { LPoly r; if (c) { r.c = {BigInt(c)}; r.v = e; } return r; }

// parse "3", "t", "-t^-1 + 2t^3", "1-t", "2*t^2" (integers, t, ^exponent, + and -)
inline LPoly parse_lpoly(const std::string& s0) {
    std::string s; for (char ch : s0) if (!isspace((unsigned char)ch)) s += ch;
    if (s.empty()) throw std::invalid_argument("laurent: empty polynomial");
    LPoly r; size_t i = 0;
    auto err = [&]() { return std::invalid_argument("laurent: cannot parse polynomial '" + s0 + "'"); };
    while (i < s.size()) {
        int sign = 1;
        if (s[i] == '+' || s[i] == '-') { if (s[i] == '-') sign = -1; ++i; }
        else if (i > 0) throw err();
        size_t j = i; while (j < s.size() && isdigit((unsigned char)s[j])) ++j;
        BigInt c(1);
        if (j > i) { c = BigInt(0); for (size_t k = i; k < j; ++k) c = c * BigInt(10) + BigInt(s[k] - '0'); }
        bool hasnum = j > i; i = j;
        if (i < s.size() && s[i] == '*') { if (!hasnum) throw err(); ++i; }
        int e = 0;
        if (i < s.size() && s[i] == 't') {
            ++i; e = 1;
            if (i < s.size() && s[i] == '^') {
                ++i; int es = 1; if (i < s.size() && (s[i] == '-' || s[i] == '+')) { if (s[i] == '-') es = -1; ++i; }
                size_t k = i; while (k < s.size() && isdigit((unsigned char)s[k])) ++k;
                if (k == i) throw err();
                e = es * std::stoi(s.substr(i, k - i)); i = k;
            }
        } else if (!hasnum) throw err();
        if (sign < 0) c = -c;
        LPoly m; if (!c.is_zero()) { m.c = {c}; m.v = e; }
        r = ladd(r, m);
    }
    return r;
}

// ---- filtered based complexes over Z[t, t^-1] ---------------------------------------------------------------
struct LaurentComplex {
    std::vector<int> dim;
    std::vector<double> birth;
    std::vector<std::vector<std::pair<int32_t, LPoly>>> bd;
    int32_t add_cell(int d, double b, std::vector<std::pair<int32_t, LPoly>> f) {
        std::map<int32_t, LPoly> m;
        for (auto& x : f) m[x.first] = ladd(m[x.first], x.second);
        std::vector<std::pair<int32_t, LPoly>> g; for (auto& kv : m) if (!kv.second.zero()) g.push_back(kv);
        dim.push_back(d); birth.push_back(b); bd.push_back(g);
        return (int32_t)dim.size() - 1;
    }
    size_t size() const { return dim.size(); }
    int max_dim() const { int m = -1; for (int d : dim) m = std::max(m, d); return m; }
    // faces exist, are earlier, one dimension lower, born no later; d d = 0 exactly over Z[t, t^-1]
    std::string validate() const {
        for (size_t c = 0; c < size(); ++c) {
            if (dim[c] < 0) return "cell " + std::to_string(c) + ": negative dimension";
            if (std::isnan(birth[c])) return "cell " + std::to_string(c) + ": birth is NaN";
            for (auto& f : bd[c]) {
                if (f.first < 0 || (size_t)f.first >= c) return "cell " + std::to_string(c) + ": a face must be an earlier cell";
                if (dim[f.first] != dim[c] - 1) return "cell " + std::to_string(c) + ": face " + std::to_string(f.first) + " has the wrong dimension";
                if (birth[f.first] > birth[c]) return "cell " + std::to_string(c) + ": born before its face " + std::to_string(f.first);
            }
            std::map<int32_t, LPoly> acc;
            for (auto& f : bd[c]) for (auto& g : bd[f.first]) acc[g.first] = ladd(acc[g.first], lmul(f.second, g.second));
            for (auto& kv : acc) if (!kv.second.zero()) return "boundary of boundary is nonzero at cell " + std::to_string(c) + " (on cell " + std::to_string(kv.first) + ": " + kv.second.str() + ")";
        }
        return "";
    }
};

// ---- algebraic Morse reduction with unit incidences +-t^k ----------------------------------------------------
// A pair sigma < tau with [tau : sigma] = u a unit of Z[t, t^-1] (+-t^k) is cancelled by Gaussian elimination:
// d rho <- d rho - [rho : sigma] u^-1 d tau for the other cofaces rho of sigma; sigma and tau leave the complex.
// This is a chain homotopy equivalence over Z[t, t^-1] (filtered, when only pairs with equal births are taken),
// so it preserves the homology over every ring above and every barcode, and with them the bad lists. The
// residual cells keep their dimension and birth; faces are no longer required to precede their cofaces by
// index (the reductions order cells by (birth, dimension)).
inline bool is_unit_lp(const LPoly& a) { return a.c.size() == 1 && a.c[0].abs() == BigInt(1); }
inline LaurentComplex laurent_reduce(const LaurentComplex& X, bool equal_births, size_t* pairs = nullptr) {
    size_t n = X.size();
    std::vector<std::map<int32_t, LPoly>> bd(n);
    std::vector<std::vector<int32_t>> cob(n);
    for (size_t c = 0; c < n; ++c) for (auto& f : X.bd[c]) { bd[c][f.first] = f.second; cob[f.first].push_back((int32_t)c); }
    std::vector<char> alive(n, 1);
    size_t np = 0;
    auto live_cob = [&](int32_t s) {
        auto& L = cob[s]; std::sort(L.begin(), L.end()); L.erase(std::unique(L.begin(), L.end()), L.end());
        size_t w = 0; for (int32_t r : L) if (alive[r] && bd[r].count(s)) L[w++] = r;
        L.resize(w); return w;
    };
    std::vector<int32_t> ord(n); for (size_t i = 0; i < n; ++i) ord[i] = (int32_t)i;
    std::stable_sort(ord.begin(), ord.end(), [&](int32_t a, int32_t b) { return X.dim[a] > X.dim[b]; });
    for (int pass = 0; pass < 8; ++pass) {
        size_t before = np;
        for (int32_t tau : ord) {
            if (!alive[tau]) continue;
            int32_t sigma = -1; size_t best = SIZE_MAX; LPoly u;
            for (auto& f : bd[tau]) {
                if (!alive[f.first] || !is_unit_lp(f.second)) continue;
                if (equal_births && X.birth[f.first] != X.birth[tau]) continue;
                size_t c = live_cob(f.first);
                if (c < best) { best = c; sigma = f.first; u = f.second; }
            }
            if (sigma < 0) continue;
            // u^-1 = (sign) t^-v
            LPoly uinv; uinv.c = {u.c[0]}; uinv.v = -u.v;
            std::vector<int32_t> rhos = cob[sigma];
            for (int32_t rho : rhos) {
                if (rho == tau || !alive[rho]) continue;
                auto it = bd[rho].find(sigma); if (it == bd[rho].end()) continue;
                LPoly m = lmul(it->second, uinv); m.c = pneg(m.c);       // d rho += m d tau
                for (auto& f : bd[tau]) {
                    LPoly v = ladd(bd[rho][f.first], lmul(m, f.second));
                    if (v.zero()) bd[rho].erase(f.first);
                    else { bd[rho][f.first] = v; cob[f.first].push_back(rho); }
                }
                if (bd[rho].count(sigma)) throw std::logic_error("laurent_reduce: elimination left the face");
            }
            alive[sigma] = alive[tau] = 0; ++np;
            for (int32_t pi : cob[tau]) if (alive[pi]) bd[pi].erase(tau);
        }
        if (np == before) break;
    }
    if (pairs) *pairs = np;
    LaurentComplex R; std::vector<int32_t> nid(n, -1);
    for (size_t c = 0; c < n; ++c) if (alive[c]) { nid[c] = (int32_t)R.size(); R.dim.push_back(X.dim[c]); R.birth.push_back(X.birth[c]); R.bd.emplace_back(); }
    for (size_t c = 0; c < n; ++c) if (alive[c]) for (auto& f : bd[c]) {
        if (!alive[f.first]) throw std::logic_error("laurent_reduce: dead face in a residual boundary");
        R.bd[nid[c]].push_back({nid[f.first], f.second});
    }
    return R;
}

typedef std::vector<std::pair<double, double>> Bars;      // per degree
typedef std::vector<Bars> LBarcode;

// ---- Smith form over K[t, t^-1] (K = Q fraction-free, or F_p) ------------------------------------------------
typedef std::vector<std::vector<Poly>> PMat;
struct SmithResult {
    size_t rank = 0;
    std::vector<Poly> factors;          // all rank invariant factors, unit-normal (units are the constant 1)
    std::vector<BigInt> scalars;        // integers the rows and columns were multiplied or divided by (Q only)
};
inline SmithResult smith(PMat A, const Dom& D) {
    SmithResult S;
    size_t R = A.size(), C = R ? A[0].size() : 0;
    for (auto& row : A) for (auto& x : row) D.trim(x);
    auto note = [&](const BigInt& x) { if (!D.p && !(x.abs() == BigInt(1))) S.scalars.push_back(x.abs()); };
    auto better = [](const Poly& x, const Poly& y) {    // a preferred pivot: lower degree, then a smaller leading coefficient
        if (deg(x) != deg(y)) return deg(x) < deg(y);
        return BigInt::cmp_mag(lc(x).a, lc(y).a) < 0;
    };
    auto rowop = [&](size_t i, size_t t, const BigInt& s, const Poly& q) {   // row i <- s row i - q row t, then / content
        for (size_t j = t; j < C; ++j) A[i][j] = psub(pscale(s, A[i][j], D), pmul(q, A[t][j], D), D);
        if (!D.p) { BigInt g; for (size_t j = t; j < C; ++j) for (auto& x : A[i][j]) g = gcd(g, x);
            if (!g.is_zero() && !(g == BigInt(1))) { note(g); for (size_t j = t; j < C; ++j) A[i][j] = pdiv_scalar(A[i][j], g, D); } }
    };
    auto colop = [&](size_t j, size_t t, const BigInt& s, const Poly& q) {
        for (size_t i = t; i < R; ++i) A[i][j] = psub(pscale(s, A[i][j], D), pmul(q, A[i][t], D), D);
        if (!D.p) { BigInt g; for (size_t i = t; i < R; ++i) for (auto& x : A[i][j]) g = gcd(g, x);
            if (!g.is_zero() && !(g == BigInt(1))) { note(g); for (size_t i = t; i < R; ++i) A[i][j] = pdiv_scalar(A[i][j], g, D); } }
    };
    for (size_t t = 0; t < R && t < C; ++t) {
        size_t pr = R, pc = C;
        for (size_t i = t; i < R; ++i) for (size_t j = t; j < C; ++j)
            if (!A[i][j].empty() && (pr == R || better(A[i][j], A[pr][pc]))) { pr = i; pc = j; }
        if (pr == R) break;
        std::swap(A[t], A[pr]);
        if (pc != t) for (size_t i = 0; i < R; ++i) std::swap(A[i][t], A[i][pc]);
        while (true) {
            for (size_t i = t + 1; i < R; ++i) {
                if (A[i][t].empty()) continue;
                BigInt s; Poly q, r; pquorem(A[i][t], A[t][t], D, s, q, r);
                if (!(s.abs() == BigInt(1))) note(lc(A[t][t]));
                rowop(i, t, s, q);
            }
            for (size_t j = t + 1; j < C; ++j) {
                if (A[t][j].empty()) continue;
                BigInt s; Poly q, r; pquorem(A[t][j], A[t][t], D, s, q, r);
                if (!(s.abs() == BigInt(1))) note(lc(A[t][t]));
                colop(j, t, s, q);
            }
            // a nonzero remainder in row or column t: it has lower degree than the pivot; move it there
            size_t bi = t, bj = t;
            for (size_t i = t + 1; i < R; ++i) if (!A[i][t].empty() && better(A[i][t], A[bi][bj])) { bi = i; bj = t; }
            for (size_t j = t + 1; j < C; ++j) if (!A[t][j].empty() && better(A[t][j], A[bi][bj])) { bi = t; bj = j; }
            if (bi != t || bj != t) { std::swap(A[t], A[bi]); if (bj != t) for (size_t i = 0; i < R; ++i) std::swap(A[i][t], A[i][bj]); continue; }
            bool clear = true;
            for (size_t i = t + 1; i < R && clear; ++i) if (!A[i][t].empty()) clear = false;
            for (size_t j = t + 1; j < C && clear; ++j) if (!A[t][j].empty()) clear = false;
            if (!clear) continue;
            // the pivot must divide every remaining entry; otherwise add that row to row t and repeat
            size_t bad = R;
            for (size_t i = t + 1; i < R && bad == R; ++i) for (size_t j = t + 1; j < C; ++j)
                if (!A[i][j].empty() && !prem(A[i][j], A[t][t], D).empty()) { bad = i; break; }
            if (bad == R) break;
            for (size_t j = t; j < C; ++j) A[t][j] = padd(A[t][j], A[bad][j], D);
        }
        note(content(A[t][t], D));                    // the content of the pivot is a scalar F_p can see
        S.factors.push_back(unit_normal(A[t][t], D));
        ++S.rank;
    }
    return S;
}
inline bool is_unit(const Poly& a) { return deg(a) == 0; }

// homology over K[t, t^-1] of the whole complex, by degree: free rank and the non-unit invariant factors
struct LHomology {
    std::vector<size_t> free;
    std::vector<std::vector<Poly>> torsion;
    std::vector<BigInt> scalars;
};
namespace detail {
// the matrix of the boundary from the given (k+1)-cells to the given k-cells, as polynomials: column j is
// multiplied by t^m so that no negative powers remain (a unit column scaling)
inline PMat boundary_poly(const LaurentComplex& X, const std::vector<int32_t>& rows, const std::vector<int32_t>& cols, const Dom& D) {
    std::map<int32_t, size_t> rix; for (size_t i = 0; i < rows.size(); ++i) rix[rows[i]] = i;
    PMat M(rows.size(), std::vector<Poly>(cols.size()));
    for (size_t j = 0; j < cols.size(); ++j) {
        int m = 0; for (auto& f : X.bd[cols[j]]) m = std::min(m, f.second.v);
        for (auto& f : X.bd[cols[j]]) {
            auto it = rix.find(f.first); if (it == rix.end()) continue;
            Poly p = pshift(f.second.c, f.second.v - m); D.trim(p);
            M[it->second][j] = p;
        }
    }
    return M;
}
}
inline LHomology laurent_homology(const LaurentComplex& X, int maxdim, const Dom& D) {
    int top = std::max(X.max_dim(), maxdim + 1);
    std::vector<std::vector<int32_t>> cells(top + 2);
    for (size_t c = 0; c < X.size(); ++c) cells[X.dim[c]].push_back((int32_t)c);
    std::vector<size_t> rank(top + 3, 0); std::vector<std::vector<Poly>> tors(top + 3);
    LHomology H;
    for (int k = 1; k <= maxdim + 1 && k <= top; ++k) {
        if (cells[k].empty() || cells[k - 1].empty()) continue;
        SmithResult S = smith(detail::boundary_poly(X, cells[k - 1], cells[k], D), D);
        rank[k] = S.rank;
        for (auto& f : S.factors) if (!is_unit(f)) tors[k].push_back(f);
        H.scalars.insert(H.scalars.end(), S.scalars.begin(), S.scalars.end());
    }
    for (int k = 0; k <= maxdim; ++k) {
        size_t nk = k <= top ? cells[k].size() : 0;
        H.free.push_back(nk - rank[k] - rank[k + 1]);
        H.torsion.push_back(tors[k + 1]);
    }
    return H;
}
// the answer over F_p[t, t^-1] predicted by the rational one: the same free ranks, the invariant factors reduced
// mod p (primitive integer polynomials, so none vanishes), units dropped
inline LHomology reduce_homology(const LHomology& H, uint64_t p) {
    Dom F(p); LHomology R; R.free = H.free;
    for (auto& T : H.torsion) {
        std::vector<Poly> v;
        for (auto& f : T) { Poly g = unit_normal(reduce_mod_p(f, p), F); if (!is_unit(g)) v.push_back(g); }
        R.torsion.push_back(v);
    }
    return R;
}
inline bool same_homology(const LHomology& a, const LHomology& b) { return a.free == b.free && a.torsion == b.torsion; }

// ---- column reduction over a field: Q(t), F_p(t), or Q[t]/(g) treated as a field -------------------------------
struct Split { Poly h; };                      // a pivot shares the factor h with g

struct FieldRun {
    LBarcode bars;                             // zero-length bars omitted
    std::vector<std::pair<Poly, Poly>> pivots; // Q(t) only: pivot values num / den, up to powers of t
};
namespace detail {
struct Col { std::vector<std::pair<int32_t, Poly>> e; };   // sorted by row (filtration position)
}

// mode: 0 = frac field (Q(t) with D.p = 0, or F_p(t)), 1 = residue field Q[t]/(g) (D.p = 0)
inline FieldRun field_reduce(const LaurentComplex& X, int maxdim, const Dom& D, const Poly* g, bool want_pivots) {
    size_t n = X.size();
    std::vector<int32_t> ord(n); for (size_t i = 0; i < n; ++i) ord[i] = (int32_t)i;
    std::stable_sort(ord.begin(), ord.end(), [&](int32_t a, int32_t b) { return X.birth[a] != X.birth[b] ? X.birth[a] < X.birth[b] : X.dim[a] < X.dim[b]; });
    std::vector<int32_t> pos(n); for (size_t i = 0; i < n; ++i) pos[ord[i]] = (int32_t)i;
    int top = X.max_dim();
    std::vector<detail::Col> R(n);
    std::vector<int32_t> low_of(n, -1);        // row position -> column position whose pivot it is
    std::vector<char> cleared(n, 0);
    std::vector<Poly> snum(n), sden(n);        // column scale s_j = snum / sden (Q(t), want_pivots)
    FieldRun out;
    auto reduce_mod_g = [&](detail::Col& c) {  // entries mod g, all scaled by the same lc(g)^E, then / content
        if (!g) return;
        int E = 0; for (auto& x : c.e) E = std::max(E, deg(x.second) - deg(*g) + 1);
        if (E == 0) return;
        std::vector<std::pair<int32_t, Poly>> ne;
        for (auto& x : c.e) {
            Poly r = prem_fixed(x.second, *g, E, D);
            if (!r.empty()) ne.push_back({x.first, r});
        }
        c.e.swap(ne);
    };
    auto normalise = [&](detail::Col& c, int32_t j) {
        if (c.e.empty()) return;
        if (g) {
            BigInt h; for (auto& x : c.e) for (auto& y : x.second) h = gcd(h, y);
            if (!(h == BigInt(1))) for (auto& x : c.e) x.second = pdiv_scalar(x.second, h, D);
            return;
        }
        Poly h; for (auto& x : c.e) { h = pgcd(h, x.second, D); if (deg(h) == 0 && (D.p || lc(h).abs() == BigInt(1))) break; }
        if (deg(h) == 0 && (D.p || h[0].abs() == BigInt(1))) return;
        for (auto& x : c.e) x.second = pexact(x.second, h, D);
        if (want_pivots) { Dom Z; sden[j] = pmul(sden[j], h, Z); }
    };
    auto cancel_scale = [&](int32_t j) {
        Dom Z; Poly h = pgcd(snum[j], sden[j], Z);
        if (deg(h) > 0 || !(h[0].abs() == BigInt(1))) { snum[j] = pexact(snum[j], h, Z); sden[j] = pexact(sden[j], h, Z); }
    };
    for (int d = top; d >= 1; --d) {
        for (size_t jj = 0; jj < n; ++jj) {
            int32_t c = ord[jj]; if (X.dim[c] != d) continue;
            if (d - 1 > maxdim) continue;
            int32_t j = (int32_t)jj;
            if (cleared[j]) continue;
            detail::Col col;
            int m = 0; for (auto& f : X.bd[c]) m = std::min(m, f.second.v);
            for (auto& f : X.bd[c]) { Poly p = pshift(f.second.c, f.second.v - m); D.trim(p); if (!p.empty()) col.e.push_back({pos[f.first], p}); }
            std::sort(col.e.begin(), col.e.end(), [](const std::pair<int32_t, Poly>& a, const std::pair<int32_t, Poly>& b) { return a.first < b.first; });
            if (want_pivots) { snum[j] = Poly{BigInt(1)}; sden[j] = Poly{BigInt(1)}; }
            reduce_mod_g(col);
            normalise(col, j);
            while (!col.e.empty()) {
                int32_t lo = col.e.back().first;
                int32_t i = low_of[lo];
                if (i < 0) break;
                Poly a = R[i].e.back().second, b = col.e.back().second;
                // col <- a col - b R[i]
                std::vector<std::pair<int32_t, Poly>> ne; size_t x = 0, y = 0; const auto& P = col.e; const auto& Q = R[i].e;
                while (x < P.size() || y < Q.size()) {
                    if (y == Q.size() || (x < P.size() && P[x].first < Q[y].first)) { Poly v = pmul(a, P[x].second, D); if (!v.empty()) ne.push_back({P[x].first, v}); ++x; }
                    else if (x == P.size() || Q[y].first < P[x].first) { Poly v = pneg(pmul(b, Q[y].second, D)); D.trim(v); if (!v.empty()) ne.push_back({Q[y].first, v}); ++y; }
                    else { Poly v = psub(pmul(a, P[x].second, D), pmul(b, Q[y].second, D), D); if (!v.empty()) ne.push_back({P[x].first, v}); ++x; ++y; }
                }
                col.e.swap(ne);
                if (want_pivots) { Dom Z; snum[j] = pmul(snum[j], a, Z); }
                reduce_mod_g(col);
                normalise(col, j);
                if (want_pivots) cancel_scale(j);
            }
            if (col.e.empty()) continue;
            int32_t lo = col.e.back().first;
            if (g) {                                       // the pivot must be invertible mod g
                Poly h = pgcd(col.e.back().second, *g, D);
                if (deg(h) > 0) throw Split{prim(h, D)};
            }
            low_of[lo] = j; cleared[lo] = 1;
            if (want_pivots) {
                Dom Z; Poly num = pmul(col.e.back().second, sden[j], Z), den = snum[j];
                Poly h = pgcd(num, den, Z); if (deg(h) > 0) { num = pexact(num, h, Z); den = pexact(den, h, Z); }
                out.pivots.push_back({num, den});   // contents kept: they carry the candidate primes
            }
            R[j] = std::move(col);
        }
    }
    out.bars.assign(std::max(0, maxdim) + 1, Bars());
    for (size_t jj = 0; jj < n; ++jj) {
        int32_t c = ord[jj];
        int k = X.dim[c]; if (k > maxdim) continue;
        if (!R[jj].e.empty()) continue;                       // a death column (or a column of a degree above)
        if (low_of[jj] >= 0) {
            double b = X.birth[c], e = X.birth[ord[low_of[jj]]];
            if (e > b) out.bars[k].push_back({b, e});
        } else out.bars[k].push_back({X.birth[c], std::numeric_limits<double>::infinity()});   // essential
    }
    for (auto& b : out.bars) std::sort(b.begin(), b.end());
    return out;
}

// ---- the certified answer ----------------------------------------------------------------------------------
struct BadIrreducible {
    Poly g;                     // squarefree, primitive, coprime to t
    bool irreducible = false;   // proved irreducible over Q (degree 1, or degree 2 with a non-square discriminant)
    LBarcode bars;              // the barcode over Q[t]/(f) for every irreducible factor f of g
};
struct LaurentResult {
    LBarcode bars;                                    // over Q(t)
    std::vector<BadIrreducible> bad_irreducibles;
    size_t candidate_polys = 0;                       // after the gcd-free basis
    std::map<uint64_t, LBarcode> bad_primes;          // p -> barcode over F_p(t)
    std::vector<uint64_t> candidate_primes;
    std::vector<std::string> unresolved;              // integers not factored (their primes are not settled)
    LHomology homology;                               // over Q[t, t^-1], the whole complex
    std::map<uint64_t, LHomology> homology_bad_primes;// p -> over F_p[t, t^-1], where it is not the reduction
    std::vector<uint64_t> homology_candidate_primes;
    bool complete = true;
    size_t cells = 0, pairs_filtered = 0, after_filtered = 0, after_unfiltered = 0;   // unit-pair reduction
};

inline bool proved_irreducible(const Poly& g) {
    if (deg(g) == 1) return true;
    if (deg(g) != 2) return false;
    BigInt disc = g[1] * g[1] - BigInt(4) * g[2] * g[0];
    if (disc.neg) return true;
    BigInt r = isqrt(disc);
    return !(r * r == disc);
}

// a gcd-free basis of the given polynomials (squarefree, primitive, t stripped, degree >= 1): every input is a
// product of powers of basis elements, and the basis elements are pairwise coprime
inline std::vector<Poly> coprime_basis(const std::vector<Poly>& in) {
    Dom Z; std::vector<Poly> B;
    for (auto& p0 : in) { Poly p = sqfree(p0); if (deg(p) >= 1) B.push_back(p); }
    std::sort(B.begin(), B.end(), [](const Poly& a, const Poly& b) { if (a.size() != b.size()) return a.size() < b.size(); for (size_t i = 0; i < a.size(); ++i) if (a[i] != b[i]) return a[i] < b[i]; return false; });
    B.erase(std::unique(B.begin(), B.end()), B.end());
    bool changed = true;
    while (changed) {
        changed = false;
        for (size_t i = 0; i < B.size() && !changed; ++i) for (size_t j = i + 1; j < B.size() && !changed; ++j) {
            Poly h = prim(pgcd(B[i], B[j], Z), Z);
            if (deg(h) < 1) continue;
            Poly a = B[i], b = B[j];
            std::vector<Poly> nb;
            for (size_t k = 0; k < B.size(); ++k) if (k != i && k != j) nb.push_back(B[k]);
            nb.push_back(h);
            Poly a2 = prim(pexact(a, h, Z), Z), b2 = prim(pexact(b, h, Z), Z);
            if (deg(a2) >= 1) nb.push_back(sqfree(a2));
            if (deg(b2) >= 1) nb.push_back(sqfree(b2));
            std::sort(nb.begin(), nb.end(), [](const Poly& x, const Poly& y) { if (x.size() != y.size()) return x.size() < y.size(); for (size_t k = 0; k < x.size(); ++k) if (x[k] != y[k]) return x[k] < y[k]; return false; });
            nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
            B.swap(nb); changed = true;
        }
    }
    return B;
}

// split off the linear factors b t - a of g (rational roots a / b) when the leading and trailing coefficients are
// small enough to enumerate their divisors; returns the factors found and the cofactor (in out.back())
inline std::vector<Poly> split_rational_roots(Poly g) {
    Dom Z; std::vector<Poly> out;
    auto divisors = [](const BigInt& x, std::vector<uint64_t>& d) {
        if (!x.fits_u64() || x.abs().to_u64() > (uint64_t)1e12) return false;
        std::map<uint64_t, int> f; factor64(x.abs().to_u64(), f);
        d = {1};
        for (auto& kv : f) { size_t m = d.size(); uint64_t pw = 1; for (int e = 1; e <= kv.second; ++e) { pw *= kv.first; for (size_t i = 0; i < m; ++i) d.push_back(d[i] * pw); } }
        return d.size() <= 2000;
    };
    while (deg(g) >= 2) {
        std::vector<uint64_t> A, B;
        if (!divisors(g[0], A) || !divisors(lc(g), B)) break;
        bool found = false;
        for (uint64_t b : B) { for (uint64_t a : A) for (int sg : {1, -1}) {
            if (gcd64(a, b) != 1) continue;
            BigInt num = BigInt((int64_t)a) * BigInt(sg), den = BigInt::from_u64(b), v, pa(1);
            std::vector<BigInt> pb(g.size(), BigInt(1)); for (size_t i = 1; i < g.size(); ++i) pb[i] = pb[i - 1] * den;
            for (size_t i = 0; i < g.size(); ++i) { v = v + g[i] * pa * pb[g.size() - 1 - i]; pa = pa * num; }
            if (!v.is_zero()) continue;
            Poly lin{-num, den};
            out.push_back(prim(lin, Z)); g = prim(pexact(g, prim(lin, Z), Z), Z); found = true; break;
        } if (found) break; }
        if (!found) break;
    }
    out.push_back(g);
    return out;
}

inline void factor_into(const BigInt& x, std::set<uint64_t>& primes, std::vector<std::string>& unresolved) {
    std::map<uint64_t, int> pr; std::vector<BigInt> un, big;
    factor_big(x, pr, un, big);
    for (auto& kv : pr) primes.insert(kv.first);
    for (auto& u : un) unresolved.push_back(u.str());
    for (auto& b : big) unresolved.push_back(b.str());
}

// The barcode over Q(t), the bad irreducibles with their barcodes, the bad primes of the generic fibre (the
// barcode over F_p(t)) and, with homology = true, the homology of the whole complex over Q[t, t^-1] and the
// primes p where the homology over F_p[t, t^-1] is not the reduction mod p of the rational answer.
inline LaurentResult laurent_persistence(const LaurentComplex& X0, int maxdim, bool classify = true, bool homology = true, bool reduce = true) {
    LaurentResult L;
    Dom Q;
    LaurentComplex Xr; if (reduce) Xr = laurent_reduce(X0, true, &L.pairs_filtered);
    const LaurentComplex& X = reduce ? Xr : X0;
    L.cells = X0.size(); L.after_filtered = X.size();
    FieldRun base = field_reduce(X, maxdim, Q, nullptr, true);
    L.bars = base.bars;
    if (classify) {
        // candidates: a gcd-free basis of the primitive parts of the numerators and denominators of the pivot
        // values, and the prime factors of their contents
        std::vector<Poly> polys; std::set<uint64_t> primes;
        for (auto& pv : base.pivots) for (const Poly* p : {&pv.first, &pv.second}) {
            Dom Z; Poly s = strip_t(prim(*p, Z));
            if (deg(s) >= 1) polys.push_back(s);
            BigInt c = content(*p, Z).abs();
            if (!(c == BigInt(1))) factor_into(c, primes, L.unresolved);
        }
        std::vector<Poly> basis0 = coprime_basis(polys), basis;
        for (auto& b : basis0) for (auto& f : split_rational_roots(b)) if (deg(f) >= 1) basis.push_back(f);
        L.candidate_polys = basis.size();
        std::vector<Poly> queue(basis.rbegin(), basis.rend());
        std::set<std::vector<std::string>> seen;
        while (!queue.empty()) {
            Poly gq = queue.back(); queue.pop_back();
            std::vector<std::string> key; for (auto& x : gq) key.push_back(x.str());
            if (seen.count(key)) continue;
            seen.insert(key);
            try {
                FieldRun fr = field_reduce(X, maxdim, Q, &gq, false);
                if (fr.bars != L.bars) { BadIrreducible b; b.g = gq; b.irreducible = proved_irreducible(gq); b.bars = fr.bars; L.bad_irreducibles.push_back(b); }
            } catch (Split& s) {
                // gq = h * (gq / h): both parts are settled separately (a further split refines them again)
                Dom Z; Poly h = sqfree(s.h), r = sqfree(pexact(gq, prim(s.h, Z), Z));
                if (deg(h) >= 1) queue.push_back(h);
                if (deg(r) >= 1) queue.push_back(r);
            }
        }
        std::sort(L.bad_irreducibles.begin(), L.bad_irreducibles.end(), [](const BadIrreducible& a, const BadIrreducible& b) { return a.g.size() != b.g.size() ? a.g.size() < b.g.size() : to_string(a.g) < to_string(b.g); });
        L.candidate_primes.assign(primes.begin(), primes.end());
        for (uint64_t p : primes) {
            FieldRun fr = field_reduce(X, maxdim, Dom(p), nullptr, false);
            if (fr.bars != L.bars) L.bad_primes[p] = fr.bars;
        }
        if (!L.unresolved.empty()) L.complete = false;
    }
    if (homology) {
        // the homology of the whole complex ignores the filtration: every unit pair may be cancelled
        LaurentComplex Xh; if (reduce) Xh = laurent_reduce(X, false);
        const LaurentComplex& Y = reduce ? Xh : X;
        L.after_unfiltered = Y.size();
        L.homology = laurent_homology(Y, maxdim, Q);
        if (classify) {
            std::set<uint64_t> primes; std::vector<std::string> un;
            for (auto& s : L.homology.scalars) factor_into(s, primes, un);
            L.homology_candidate_primes.assign(primes.begin(), primes.end());
            for (auto& u : un) L.unresolved.push_back(u);
            if (!un.empty()) L.complete = false;
            for (uint64_t p : primes) {
                LHomology Hp = laurent_homology(Y, maxdim, Dom(p));
                if (!same_homology(Hp, reduce_homology(L.homology, p))) L.homology_bad_primes[p] = Hp;
            }
        }
    }
    return L;
}

// ---- input and output -----------------------------------------------------------------------------------------
// text format (perseus laurent): one cell per line, "dim birth face:poly face:poly ...", faces are 0-based line
// numbers of earlier cells, polynomials without spaces ("1-t", "-t^-1+2t^3", "3"); # comments, blank lines skipped
inline LaurentComplex read_laurent(std::istream& is) {
    LaurentComplex X; std::string line; int ln = 0;
    while (std::getline(is, line)) {
        ++ln;
        size_t h = line.find('#'); if (h != std::string::npos) line = line.substr(0, h);
        std::istringstream ss(line); std::vector<std::string> t; std::string w; while (ss >> w) t.push_back(w);
        if (t.empty()) continue;
        if (t.size() < 2) throw std::runtime_error("laurent: line " + std::to_string(ln) + ": expected dim birth face:poly ...");
        int d = std::stoi(t[0]); double b = std::stod(t[1]);
        std::vector<std::pair<int32_t, LPoly>> f;
        for (size_t i = 2; i < t.size(); ++i) {
            size_t c = t[i].find(':');
            if (c == std::string::npos) throw std::runtime_error("laurent: line " + std::to_string(ln) + ": a face is face:poly");
            long long y = std::stoll(t[i].substr(0, c));
            if (y < 0 || y >= (long long)X.size()) throw std::runtime_error("laurent: line " + std::to_string(ln) + ": a face must be an earlier cell");
            f.push_back({(int32_t)y, parse_lpoly(t[i].substr(c + 1))});
        }
        X.add_cell(d, b, f);
    }
    return X;
}
inline std::string num_json(double x) {
    if (std::isinf(x)) return "null";
    char b[64]; snprintf(b, sizeof b, "%.17g", x); return b;
}
inline std::string bars_json(const LBarcode& B) {
    std::string s = "[";
    for (size_t k = 0; k < B.size(); ++k) {
        s += k ? ", [" : "[";
        for (size_t i = 0; i < B[k].size(); ++i) s += (i ? ", [" : "[") + num_json(B[k][i].first) + ", " + num_json(B[k][i].second) + "]";
        s += "]";
    }
    return s + "]";
}
inline std::string homology_json(const LHomology& H) {
    std::string s = "[";
    for (size_t k = 0; k < H.free.size(); ++k) {
        s += k ? ", " : ""; s += "{\"degree\": " + std::to_string(k) + ", \"free\": " + std::to_string(H.free[k]) + ", \"torsion\": [";
        for (size_t i = 0; i < H.torsion[k].size(); ++i) s += (i ? ", \"" : "\"") + to_string(H.torsion[k][i]) + "\"";
        s += "]}";
    }
    return s + "]";
}
inline std::string to_json(const LaurentResult& L) {
    std::string s = "{\n  \"coefficients\": \"Z[t, t^-1]\",\n";
    s += "  \"barcode_over\": \"Q(t)\",\n  \"barcode\": " + bars_json(L.bars) + ",\n";
    s += "  \"bad_irreducibles\": [";
    for (size_t i = 0; i < L.bad_irreducibles.size(); ++i) {
        auto& b = L.bad_irreducibles[i];
        s += std::string(i ? ",\n    " : "\n    ") + "{\"polynomial\": \"" + to_string(b.g) + "\", \"irreducible\": " + (b.irreducible ? "true" : "false") +
             ", \"every_irreducible_factor_of\": " + (b.irreducible ? "false" : "true") + ", \"barcode\": " + bars_json(b.bars) + "}";
    }
    s += "],\n  \"bad_primes\": {";
    bool first = true;
    for (auto& kv : L.bad_primes) { s += std::string(first ? "" : ", ") + "\"" + std::to_string(kv.first) + "\": " + bars_json(kv.second); first = false; }
    s += "},\n  \"candidates\": {\"polynomials\": " + std::to_string(L.candidate_polys) + ", \"primes\": [";
    for (size_t i = 0; i < L.candidate_primes.size(); ++i) s += (i ? ", " : "") + std::to_string(L.candidate_primes[i]);
    s += "]},\n  \"homology_over\": \"Q[t, t^-1]\",\n  \"homology\": " + homology_json(L.homology) + ",\n";
    s += "  \"homology_bad_primes\": {";
    first = true;
    for (auto& kv : L.homology_bad_primes) { s += std::string(first ? "" : ", ") + "\"" + std::to_string(kv.first) + "\": " + homology_json(kv.second); first = false; }
    s += "},\n  \"homology_candidate_primes\": [";
    for (size_t i = 0; i < L.homology_candidate_primes.size(); ++i) s += (i ? ", " : "") + std::to_string(L.homology_candidate_primes[i]);
    s += "],\n  \"unresolved\": [";
    for (size_t i = 0; i < L.unresolved.size(); ++i) s += (i ? ", \"" : "\"") + L.unresolved[i] + "\"";
    s += std::string("],\n  \"complete\": ") + (L.complete ? "true" : "false") + "\n}\n";
    return s;
}

} // namespace laurent
} // namespace perseus
