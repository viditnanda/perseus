// pc_group.hpp -- integral homology and mod-p cohomology rings of a finite group given by a polycyclic
// generating sequence, through the Anick resolution of its collection rewriting system.
//
// Input: the multiplication table of G and a pcgs f_1, ..., f_k (element ids) with relative orders
// r_1, ..., r_k. Every element is uniquely f_1^a_1 ... f_k^a_k with 0 <= a_i < r_i, and with
// G_i = <f_i, ..., f_k>, both f_i^r_i and f_i^-1 f_j f_i (j > i) lie in G_{i+1}. The collected words are then
// the normal forms of a complete rewriting system whose minimal non-normal words ("tips") are f_i^r_i and
// f_j f_i (j > i).
//
// The Anick n-chains of that system are the critical cells of the algebraic Morse matching it induces on the
// bar resolution (Brown 1989; Anick 1986; Skoldberg 2006; Joellenbeck-Welker 2009). Here they are the
// non-increasing letter sequences j_1 >= ... >= j_n, so there are C(n+k-1, k-1) of them. A chain is stored as
// segments: the first is f_{j_1}; a smaller letter j is the segment f_j; repeating the previous letter j after a
// segment f_j^r gives the segment f_j^(r_j - r).
//
// The resolution P_n = Z[chains_n] (x) ZG of right ZG-modules is built by the contracting-homotopy recursion
//   d_1(f) = () (x) f - () (x) e
//   d_{n+1}(c t) = c (x) t - h(d_n(c) . t)
//   h(c (x) g) = (c s) (x) u - h(d(c s) . u - c (x) g)   when g = s u extends c to a chain, and 0 otherwise,
// memoised on the monomials c (x) g. That is at most |G| C(n+k-1, k-1) states in degree n: polynomial in the
// degree, where any enumeration of the bar construction is exponential.
//
// Integral homology: |G| kills H_n(G; Z) for n >= 1, so the only primes are those dividing |G| and, at a prime q
// with q^v || |G|, every nonzero elementary divisor of the augmented differentials has q-valuation <= v. The
// q-part is therefore read off a Smith form over Z/q^(v+1) (sparse elimination, least valuation first); this is
// exact, with no rational arithmetic and no candidate primes.
//
// Mod-p cohomology ring: F_p (x) P is split once across degrees by Gaussian elimination (fp_split), which gives a
// basis of cocycles and a dual basis of cycles at the size of H^*. Products by Yoneda composition: a cocycle alpha
// of degree a lifts through h to a chain map alpha~ : P_{a+*} -> P_*, and (alpha . beta)(c) = beta(eps alpha~_b(c)).
// By default only the minimal generators are lifted (their products span the decomposables, by associativity).
// This is the cup product: pulled back to the normalised bar resolution through the comparison map built from h,
// it agrees with the Alexander-Whitney cup product of the pulled-back cocycles up to a coboundary, with no sign
// (tests/t_pcgroup.cpp checks this at p = 2, 3, 5, in degrees where the sign (-1)^(ab) would be visible).
#pragma once
#include <vector>
#include <deque>
#include <unordered_map>
#include <map>
#include <tuple>
#include <string>
#include <stdexcept>
#include <cstdint>
#include <algorithm>
#include <chrono>
#include <numeric>
#include <queue>
#include <climits>
#include <cstdio>
#include <cstdlib>

namespace perseus {
namespace pcg {

typedef int64_t i64;
typedef uint64_t u64;
typedef std::vector<std::pair<u64, i64>> Vec;          // sorted keys, coefficients in [1, M)

inline i64 mulmod(i64 a, i64 b, i64 M) { return (i64)((__int128)a * b % M); }
inline i64 inv_mod(i64 a, i64 M) {
    i64 t = 0, nt = 1, r = M, nr = a % M; if (nr < 0) nr += M;
    while (nr) { i64 q = r / nr, x = t - q * nt; t = nt; nt = x; x = r - q * nr; r = nr; nr = x; }
    if (r != 1) throw std::logic_error("pc_group: inverse of a non-unit");
    return t < 0 ? t + M : t;
}
inline void normalize(std::vector<std::pair<u64, i64>>& acc, i64 M, Vec& out) {
    std::sort(acc.begin(), acc.end(), [](const std::pair<u64, i64>& a, const std::pair<u64, i64>& b) { return a.first < b.first; });
    out.clear();
    for (size_t i = 0; i < acc.size();) {
        u64 key = acc[i].first; i64 s = 0;
        for (; i < acc.size() && acc[i].first == key; ++i) { s += acc[i].second; if (s >= M) s -= M; }
        if (s) out.push_back({key, s});
    }
}

// ---- the group ----------------------------------------------------------------------------------------
struct PcGroup {
    int n = 0, k = 0, e = 0;
    std::vector<int32_t> T, inv;                 // multiplication table (n x n), inverses
    std::vector<int> gen, rel, stride;
    std::vector<int32_t> of_code, code;           // mixed-radix exponent code <-> element
    std::vector<int32_t> exps;                    // n x k exponent vectors
    std::vector<int> first_letter, first_run, length;

    int mul(int a, int b) const { return T[(size_t)a * n + b]; }
    int ex(int g, int i) const { return exps[(size_t)g * k + i]; }
    int power(int i, int a) const { return of_code[(size_t)a * stride[i]]; }

    PcGroup() {}
    PcGroup(const std::vector<std::vector<int>>& mult, const std::vector<int>& pcgs, std::vector<int> relative_orders, int identity) {
        n = (int)mult.size(); e = identity;
        if (n < 1) throw std::invalid_argument("pc_group: empty multiplication table");
        if (e < 0 || e >= n) throw std::invalid_argument("pc_group: identity out of range");
        T.assign((size_t)n * n, 0);
        for (int a = 0; a < n; ++a) {
            if ((int)mult[a].size() != n) throw std::invalid_argument("pc_group: the multiplication table must be square");
            for (int b = 0; b < n; ++b) { int x = mult[a][b]; if (x < 0 || x >= n) throw std::invalid_argument("pc_group: table entry out of range"); T[(size_t)a * n + b] = x; }
        }
        for (int a = 0; a < n; ++a) if (mul(e, a) != a || mul(a, e) != a) throw std::invalid_argument("pc_group: the identity is not an identity");
        std::vector<char> seen(n);
        for (int a = 0; a < n; ++a) {             // Latin square: rows and columns are permutations
            std::fill(seen.begin(), seen.end(), 0); for (int b = 0; b < n; ++b) { if (seen[mul(a, b)]) throw std::invalid_argument("pc_group: not a group (a row of the table repeats an element)"); seen[mul(a, b)] = 1; }
            std::fill(seen.begin(), seen.end(), 0); for (int b = 0; b < n; ++b) { if (seen[mul(b, a)]) throw std::invalid_argument("pc_group: not a group (a column of the table repeats an element)"); seen[mul(b, a)] = 1; }
        }
        inv.assign(n, -1);
        for (int a = 0; a < n; ++a) for (int b = 0; b < n; ++b) if (mul(a, b) == e) { inv[a] = b; break; }
        gen = pcgs; k = (int)gen.size();
        for (int g : gen) if (g < 0 || g >= n || g == e) throw std::invalid_argument("pc_group: a pcgs element is out of range or the identity");
        if (relative_orders.empty()) {            // infer: r_i is the least r with f_i^r in G_{i+1}
            relative_orders.assign(k, 0);
            std::vector<char> in(n, 0); in[e] = 1; std::vector<int> elems{e};
            for (int i = k - 1; i >= 0; --i) {
                int r = 1, x = gen[i];
                while (!in[x]) { x = mul(x, gen[i]); ++r; if (r > n) throw std::invalid_argument("pc_group: could not infer the relative orders"); }
                relative_orders[i] = r;
                std::vector<int> next; int y = e;
                for (int a = 0; a < r; ++a) { for (int z : elems) next.push_back(mul(y, z)); y = mul(y, gen[i]); }
                for (int z : next) in[z] = 1;
                elems.swap(next);
            }
        }
        rel = relative_orders;
        if ((int)rel.size() != k) throw std::invalid_argument("pc_group: relative_orders must have one entry per pcgs element");
        long long prod = 1; for (int r : rel) { if (r < 2) throw std::invalid_argument("pc_group: relative orders must be at least 2"); prod *= r; if (prod > n) break; }
        if (prod != n) throw std::invalid_argument("pc_group: the relative orders do not multiply to |G|");
        stride.assign(k, 1); for (int i = k - 2; i >= 0; --i) stride[i] = stride[i + 1] * rel[i + 1];
        of_code.assign(n, -1); code.assign(n, -1); exps.assign((size_t)n * k, 0);
        std::vector<std::vector<int>> pw(k);
        for (int i = 0; i < k; ++i) { pw[i].push_back(e); for (int a = 1; a < rel[i]; ++a) pw[i].push_back(mul(pw[i].back(), gen[i])); }
        for (int c = 0; c < n; ++c) {
            int g = e, rem = c; std::vector<int> dig(k);
            for (int i = 0; i < k; ++i) { dig[i] = rem / stride[i]; rem %= stride[i]; g = mul(g, pw[i][dig[i]]); }
            for (int i = 0; i < k; ++i) exps[(size_t)g * k + i] = dig[i];
            if (code[g] != -1) throw std::invalid_argument("pc_group: not a pcgs (two collected words give the same element)");
            code[g] = c; of_code[c] = g;
        }
        // Light's associativity test: the s with (x s) y = x (s y) for all x, y are closed under products, and every
        // element is a product of pcgs elements (the collected words above), so checking s in the pcgs is exact
        for (int s2 : gen) for (int x = 0; x < n; ++x) { int xs = mul(x, s2); for (int y = 0; y < n; ++y) if (mul(xs, y) != mul(x, mul(s2, y))) throw std::invalid_argument("pc_group: the table is not associative"); }
        // G_{i+1} = the elements whose exponents vanish at 1..i: check f_i^r_i and f_i^-1 f_j f_i (j > i) lie there
        auto in_tail = [&](int g, int i) { for (int j = 0; j <= i; ++j) if (ex(g, j)) return false; return true; };
        for (int i = 0; i < k; ++i) {
            if (!in_tail(mul(pw[i][rel[i] - 1], gen[i]), i)) throw std::invalid_argument("pc_group: not a pcgs (f_i^r_i is not in G_{i+1})");
            for (int j = i + 1; j < k; ++j)
                if (!in_tail(mul(mul(inv[gen[i]], gen[j]), gen[i]), i)) throw std::invalid_argument("pc_group: not a pcgs (G_{i+1} is not normalised by f_i)");
        }
        first_letter.assign(n, -1); first_run.assign(n, 0); length.assign(n, 0);
        for (int g = 0; g < n; ++g) for (int i = 0; i < k; ++i) { int a = ex(g, i); length[g] += a; if (a && first_letter[g] < 0) { first_letter[g] = i; first_run[g] = a; } }
    }
    int suffix(int g, int t) const {              // the element after the first t letters of the collected word
        int c = 0;
        for (int i = 0; i < k; ++i) { int a = ex(g, i), take = std::min(a, t); t -= take; c += (a - take) * stride[i]; }
        return of_code[c];
    }
};

// ---- Anick chains --------------------------------------------------------------------------------------
struct Chains {
    int k = 0;
    std::vector<std::vector<int32_t>> parent, seg, lastL, lastR, child;   // child[n][c * k + j]
    size_t size(int n) const { return parent[n].size(); }
    void build(const PcGroup& G, int D) {
        k = G.k;
        parent.assign(D + 1, {}); seg = lastL = lastR = child = parent;
        parent[0] = {-1}; seg[0] = {G.e}; lastL[0] = {k}; lastR[0] = {0};
        for (int n = 0; n < D; ++n) {
            child[n].assign(size(n) * k, -1);
            for (size_t c = 0; c < size(n); ++c) {
                int L = lastL[n][c], R = lastR[n][c];
                for (int j = 0; j < k && j <= L; ++j) {
                    int run = (j < L) ? 1 : G.rel[j] - R;
                    child[n][c * k + j] = (int32_t)parent[n + 1].size();
                    parent[n + 1].push_back((int32_t)c); seg[n + 1].push_back(G.power(j, run));
                    lastL[n + 1].push_back(j); lastR[n + 1].push_back(run);
                }
            }
        }
        child[D].assign(size(D) * k, -1);
    }
};

// ---- the resolution ------------------------------------------------------------------------------------
struct Anick {
    const PcGroup& G; Chains C; i64 M; int D;
    std::vector<std::vector<Vec>> d;                       // d[n][c] in P_{n-1}, n = 1..D
    std::vector<std::unordered_map<u64, int32_t>> memo;    // h_n: monomial -> index into pool[n] (-1: in progress)
    std::vector<std::deque<Vec>> pool;

    u64 key(size_t c, int g) const { return (u64)c * (u64)G.n + (u64)g; }

    Anick(const PcGroup& G_, int D_, i64 M_) : G(G_), M(M_), D(D_) {
        C.build(G, D); d.assign(D + 1, {}); memo.assign(D + 1, {}); pool.assign(D + 1, {});
        for (int n = 1; n <= D; ++n) {
            d[n].resize(C.size(n));
            for (size_t c = 0; c < C.size(n); ++c) {
                int p = C.parent[n][c], t = C.seg[n][c];
                std::vector<std::pair<u64, i64>> acc;
                if (n == 1) { acc.push_back({key(0, t), 1}); acc.push_back({key(0, G.e), M - 1}); }
                else {
                    Vec X; rmul(d[n - 1][p], t, X);
                    Vec H; hvec(n - 2, X, H);
                    acc.push_back({key(p, t), 1});
                    for (auto& kv : H) acc.push_back({kv.first, M - kv.second});
                }
                normalize(acc, M, d[n][c]);
            }
        }
    }
    void rmul(const Vec& v, int u, Vec& out) const {
        out.clear(); out.reserve(v.size());
        for (auto& kv : v) { u64 c = kv.first / (u64)G.n; int g = (int)(kv.first % (u64)G.n); out.push_back({key(c, G.mul(g, u)), kv.second}); }
        std::sort(out.begin(), out.end(), [](const std::pair<u64, i64>& a, const std::pair<u64, i64>& b) { return a.first < b.first; });
    }
    // does c (x) g (c a chain of degree n) extend? -> (child chain of degree n+1, u)
    bool ext(int n, size_t c, int g, int& ch, int& u) const {
        if (g == G.e) return false;
        int j = G.first_letter[g], s = G.first_run[g], i = C.lastL[n][c], r = C.lastR[n][c], t, letter;
        if (j < i) { t = 1; letter = j; }
        else if (j == i && r + s >= G.rel[i]) { t = G.rel[i] - r; letter = i; }
        else return false;
        ch = C.child[n][c * C.k + letter]; u = G.suffix(g, t);
        if (ch < 0) throw std::logic_error("pc_group: missing chain extension");
        return true;
    }
    int32_t store(int n, Vec&& v) { pool[n].push_back(std::move(v)); return (int32_t)pool[n].size() - 1; }
    // h_n on a monomial: index into pool[n]; iterative, with a guard against non-termination
    int32_t h(int n, u64 k0) {
        auto& mm = memo[n];
        auto it = mm.find(k0);
        if (it != mm.end()) { if (it->second < 0) throw std::logic_error("pc_group: the contracting homotopy does not terminate (is the input a pcgs?)"); return it->second; }
        if (n + 1 > D) throw std::logic_error("pc_group: h needs a higher degree than was built");
        struct Frame { u64 key, lead; Vec D; size_t pos; };
        std::vector<Frame> st;
        auto open = [&](u64 k) -> bool {
            size_t c = k / (u64)G.n; int g = (int)(k % (u64)G.n), ch, u;
            if (!ext(n, c, g, ch, u)) { mm[k] = store(n, Vec()); return false; }
            Frame f; f.key = k; f.lead = key(ch, u) ; f.pos = 0;
            rmul(d[n + 1][ch], u, f.D);
            auto p = std::lower_bound(f.D.begin(), f.D.end(), std::make_pair(k, (i64)0), [](const std::pair<u64, i64>& a, const std::pair<u64, i64>& b) { return a.first < b.first; });
            if (p == f.D.end() || p->first != k || p->second != 1) throw std::logic_error("pc_group: leading coefficient is not 1 (is the input a pcgs?)");
            f.D.erase(p);
            mm[k] = -1;
            st.push_back(std::move(f));
            return true;
        };
        if (!open(k0)) return mm[k0];
        while (!st.empty()) {
            bool pushed = false;
            while (st.back().pos < st.back().D.size()) {
                u64 k2 = st.back().D[st.back().pos].first;
                auto it2 = mm.find(k2);
                if (it2 == mm.end()) { if (open(k2)) { pushed = true; break; } ++st.back().pos; continue; }
                if (it2->second < 0) throw std::logic_error("pc_group: the contracting homotopy does not terminate (is the input a pcgs?)");
                ++st.back().pos;
            }
            if (pushed) continue;
            Frame& f = st.back();
            std::vector<std::pair<u64, i64>> acc; acc.push_back({f.lead, 1});
            for (auto& kv : f.D) { const Vec& hv = pool[n][mm[kv.first]]; for (auto& kw : hv) acc.push_back({kw.first, M - mulmod(kv.second, kw.second, M)}); }
            Vec res; normalize(acc, M, res);
            mm[f.key] = store(n, std::move(res));
            st.pop_back();
        }
        return mm[k0];
    }
    void hvec(int n, const Vec& X, Vec& out) {
        std::vector<std::pair<u64, i64>> acc;
        for (auto& kv : X) { int32_t ix = h(n, kv.first); for (auto& kw : pool[n][ix]) acc.push_back({kw.first, mulmod(kv.second, kw.second, M)}); }
        normalize(acc, M, out);
    }
    // the augmented differential eps(d_n) as sparse columns over the chains of degree n-1
    std::vector<std::vector<std::pair<int32_t, i64>>> eps(int n) const {
        std::vector<std::vector<std::pair<int32_t, i64>>> cols(C.size(n));
        for (size_t c = 0; c < C.size(n); ++c) {
            std::map<int32_t, i64> agg;
            for (auto& kv : d[n][c]) { i64& s = agg[(int32_t)(kv.first / (u64)G.n)]; s += kv.second; if (s >= M) s -= M; }
            for (auto& kv : agg) if (kv.second) cols[c].push_back(kv);
        }
        return cols;
    }
    size_t states() const { size_t s = 0; for (auto& m : memo) s += m.size(); return s; }
};

// ---- Smith form over Z/q^m: valuations of the nonzero elementary divisors -------------------------------
inline std::vector<int> local_elementary_valuations(std::vector<std::vector<std::pair<int32_t, i64>>> cols, size_t nrows, i64 q, int m) {
    i64 Q = 1; for (int i = 0; i < m; ++i) Q *= q;
    auto val = [&](i64 x) { int v = 0; while (v < m && x % q == 0) { x /= q; ++v; } return v; };
    for (auto& col : cols) {                      // reduce mod Q, drop zeros
        std::vector<std::pair<int32_t, i64>> c2;
        for (auto& kv : col) { i64 x = kv.second % Q; if (x < 0) x += Q; if (x) c2.push_back({kv.first, x}); }
        col.swap(c2);
    }
    std::vector<std::vector<int32_t>> rowcols(nrows);
    for (size_t j = 0; j < cols.size(); ++j) for (auto& kv : cols[j]) rowcols[kv.first].push_back((int32_t)j);
    std::vector<char> dead(cols.size(), 0), rowdead(nrows, 0);
    std::vector<int> version(cols.size(), 0);
    typedef std::tuple<int, size_t, int32_t, int> HE;   // (min valuation, size, column, version)
    std::priority_queue<HE, std::vector<HE>, std::greater<HE>> heap;
    auto minval = [&](size_t j) { int v = m; for (auto& kv : cols[j]) v = std::min(v, val(kv.second)); return v; };
    for (size_t j = 0; j < cols.size(); ++j) if (!cols[j].empty()) heap.push(HE(minval(j), cols[j].size(), (int32_t)j, 0));
    std::vector<int> out;
    std::vector<std::pair<int32_t, i64>> tmp;
    while (!heap.empty()) {
        HE top = heap.top(); heap.pop();
        int32_t j = std::get<2>(top);
        if (dead[j] || std::get<3>(top) != version[j] || cols[j].empty()) continue;
        int v = std::get<0>(top);
        int32_t r = -1; size_t best = SIZE_MAX; i64 a = 0;
        for (auto& kv : cols[j]) if (val(kv.second) == v && rowcols[kv.first].size() < best) { best = rowcols[kv.first].size(); r = kv.first; a = kv.second; }
        i64 qv = 1; for (int i = 0; i < v; ++i) qv *= q;
        i64 ainv = inv_mod((a / qv) % Q, Q);
        out.push_back(v);
        std::vector<int32_t> others = rowcols[r];
        std::sort(others.begin(), others.end()); others.erase(std::unique(others.begin(), others.end()), others.end());
        for (int32_t j2 : others) {
            if (j2 == j || dead[j2]) continue;
            auto& c2 = cols[j2];
            auto p = std::lower_bound(c2.begin(), c2.end(), std::make_pair(r, (i64)0), [](const std::pair<int32_t, i64>& x, const std::pair<int32_t, i64>& y) { return x.first < y.first; });
            if (p == c2.end() || p->first != r) continue;
            i64 f = mulmod(p->second / qv, ainv, Q);
            tmp.clear(); size_t x = 0, y = 0; const auto& cj = cols[j];
            while (x < c2.size() || y < cj.size()) {
                if (y == cj.size() || (x < c2.size() && c2[x].first < cj[y].first)) { tmp.push_back(c2[x]); ++x; }
                else {
                    i64 s = (Q - mulmod(f, cj[y].second, Q)) % Q; int32_t row = cj[y].first;
                    if (x < c2.size() && c2[x].first == row) { s += c2[x].second; if (s >= Q) s -= Q; ++x; }
                    else rowcols[row].push_back(j2);
                    if (s) tmp.push_back({row, s});
                    ++y;
                }
            }
            c2.swap(tmp);
            ++version[j2];
            if (!c2.empty()) heap.push(HE(minval(j2), c2.size(), j2, version[j2]));
        }
        dead[j] = 1; rowdead[r] = 1; rowcols[r].clear();
    }
    return out;
}

// ---- unit cancellation across a chain complex over Z/q^m ------------------------------------------------
// Gaussian elimination on a unit entry [tau : sigma] of d_n is a chain homotopy equivalence over the local ring:
// d_n becomes its Schur complement, the row tau of d_{n+1} and the column sigma of d_{n-1} are deleted
// (Gaussian elimination lemma; the matrices keep their original bases on the surviving cells). Pivots are
// taken by least Markowitz cost (|column| - 1)(|row| - 1), so collapses and coreductions, which cause no fill,
// come first. What is left has no unit entries.
struct LocalMat {
    size_t nrows = 0;
    std::vector<std::vector<std::pair<int32_t, i64>>> cols;
    std::vector<int32_t> rowcnt; std::vector<std::vector<int32_t>> rowcols;
    std::vector<char> coldead, rowdead; std::vector<int> version;
    void init(size_t nr, std::vector<std::vector<std::pair<int32_t, i64>>>&& c) {
        nrows = nr; cols = std::move(c); rowcnt.assign(nr, 0); rowcols.assign(nr, {}); coldead.assign(cols.size(), 0); rowdead.assign(nr, 0); version.assign(cols.size(), 0);
        for (size_t j = 0; j < cols.size(); ++j) for (auto& kv : cols[j]) { ++rowcnt[kv.first]; rowcols[kv.first].push_back((int32_t)j); }
    }
    void kill_row(int32_t r) {                    // delete a row from every column
        if (rowdead[r]) return;
        rowdead[r] = 1;
        for (int32_t j : rowcols[r]) {
            if (coldead[j]) continue;
            auto& c = cols[j];
            auto p = std::lower_bound(c.begin(), c.end(), std::make_pair(r, (i64)0), [](const std::pair<int32_t, i64>& x, const std::pair<int32_t, i64>& y) { return x.first < y.first; });
            if (p != c.end() && p->first == r) { c.erase(p); ++version[j]; }
        }
        rowcols[r].clear(); rowcnt[r] = 0;
    }
    void kill_col(int32_t j) {
        if (coldead[j]) return;
        coldead[j] = 1; for (auto& kv : cols[j]) --rowcnt[kv.first]; cols[j].clear();
    }
};
// cancel unit pivots in M; returns the pairs (row sigma, column tau)
inline std::vector<std::pair<int32_t, int32_t>> unit_cancel(LocalMat& M, i64 q, i64 Q) {
    typedef std::tuple<long long, int32_t, int> HE;
    std::priority_queue<HE, std::vector<HE>, std::greater<HE>> heap;
    auto cost = [&](int32_t j, int32_t* row) {
        long long best = LLONG_MAX; long long cs = (long long)M.cols[j].size() - 1;
        for (auto& kv : M.cols[j]) if (kv.second % q) { long long c = cs * (long long)(M.rowcnt[kv.first] - 1); if (c < best) { best = c; if (row) *row = kv.first; } }
        return best;
    };
    for (size_t j = 0; j < M.cols.size(); ++j) if (!M.coldead[j]) { long long c = cost((int32_t)j, nullptr); if (c != LLONG_MAX) heap.push(HE(c, (int32_t)j, M.version[j])); }
    std::vector<std::pair<int32_t, int32_t>> pairs;
    std::vector<std::pair<int32_t, i64>> tmp;
    while (!heap.empty()) {
        HE top = heap.top(); heap.pop();
        int32_t j = std::get<1>(top);
        if (M.coldead[j] || std::get<2>(top) != M.version[j]) continue;
        int32_t r = -1; long long c = cost(j, &r);
        if (c == LLONG_MAX) continue;
        if (!heap.empty() && c > std::get<0>(heap.top())) { heap.push(HE(c, j, M.version[j])); continue; }
        const auto& cj = M.cols[j];
        i64 a = 0; for (auto& kv : cj) if (kv.first == r) a = kv.second;
        i64 ainv = inv_mod(a, Q);
        std::vector<int32_t> others = M.rowcols[r];
        std::sort(others.begin(), others.end()); others.erase(std::unique(others.begin(), others.end()), others.end());
        for (int32_t j2 : others) {
            if (j2 == j || M.coldead[j2]) continue;
            auto& c2 = M.cols[j2];
            auto p = std::lower_bound(c2.begin(), c2.end(), std::make_pair(r, (i64)0), [](const std::pair<int32_t, i64>& x, const std::pair<int32_t, i64>& y) { return x.first < y.first; });
            if (p == c2.end() || p->first != r) continue;
            i64 f = mulmod(p->second, ainv, Q);
            tmp.clear(); size_t x = 0, y = 0;
            while (x < c2.size() || y < cj.size()) {
                if (y == cj.size() || (x < c2.size() && c2[x].first < cj[y].first)) { tmp.push_back(c2[x]); ++x; }
                else {
                    i64 s = (Q - mulmod(f, cj[y].second, Q)) % Q; int32_t row = cj[y].first; bool had = false;
                    if (x < c2.size() && c2[x].first == row) { s += c2[x].second; if (s >= Q) s -= Q; ++x; had = true; }
                    if (s) { tmp.push_back({row, s}); if (!had) { ++M.rowcnt[row]; M.rowcols[row].push_back(j2); } }
                    else if (had) --M.rowcnt[row];
                    ++y;
                }
            }
            c2.swap(tmp); ++M.version[j2];
            long long c3 = cost(j2, nullptr); if (c3 != LLONG_MAX) heap.push(HE(c3, j2, M.version[j2]));
        }
        pairs.push_back({r, j});
        M.kill_col(j); M.rowdead[r] = 1; M.rowcols[r].clear(); M.rowcnt[r] = 0;
    }
    return pairs;
}

// ---- results -------------------------------------------------------------------------------------------
struct GroupHomology {
    std::vector<size_t> free;
    std::vector<std::vector<u64>> torsion;            // invariant factors, each dividing the next
};
struct CohomologyRing {
    int p = 0;
    std::vector<size_t> dims;                         // dim H^n(G; F_p)
    std::vector<size_t> indecomposables;              // minimal ring generators in degree n (index n; n >= 1)
    std::vector<std::pair<int, int>> generators;      // (degree, basis index) of a choice of minimal generators
    bool all_products = false;                        // products: all (a, i, b, j) with a <= b; else (a, i) a generator
    std::map<std::tuple<int, int, int, int>, std::vector<uint32_t>> products;   // (a, i, b, j): coordinates of alpha_ai . beta_bj
    std::vector<std::vector<std::vector<std::pair<int32_t, uint32_t>>>> basis;   // basis[n][i]: a cocycle on the Anick n-chains
};
struct PcResult {
    bool has_homology = false, has_ring = false;
    GroupHomology H; CohomologyRing R;
    std::vector<size_t> chains;                       // Anick chains per degree
    size_t states = 0; double seconds = 0;
};

inline std::vector<std::pair<i64, int>> factor_small(i64 x) {
    std::vector<std::pair<i64, int>> f;
    for (i64 p = 2; p * p <= x; ++p) if (x % p == 0) { int v = 0; while (x % p == 0) { x /= p; ++v; } f.push_back({p, v}); }
    if (x > 1) f.push_back({x, 1});
    return f;
}

// H_n(G; Z) for n = 0..maxdim
inline GroupHomology group_homology(const PcGroup& G, int maxdim, size_t* chains_out = nullptr, size_t* states_out = nullptr) {
    auto primes = factor_small(G.n);
    i64 M = 1; for (auto& pv : primes) for (int i = 0; i <= pv.second; ++i) M *= pv.first;
    if (G.n == 1) M = 2;
    bool stats = std::getenv("PERSEUS_PC_STATS") != nullptr;
    auto t0 = std::chrono::steady_clock::now();
    Anick A(G, maxdim + 1, M);
    if (stats) fprintf(stderr, "  [anick resolution to degree %d: %.3f s, %zu states]\n", maxdim + 1, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(), A.states());
    std::vector<std::vector<std::vector<std::pair<int32_t, i64>>>> E(maxdim + 2);
    for (int n = 1; n <= maxdim + 1; ++n) E[n] = A.eps(n);
    std::vector<size_t> rank(maxdim + 3, 0);
    std::vector<std::vector<std::vector<int>>> vals(maxdim + 3);  // vals[n][prime index]: residual valuations of d_n
    for (int n = 1; n <= maxdim + 1; ++n) vals[n].resize(primes.size());
    for (size_t pi = 0; pi < primes.size(); ++pi) {
        i64 q = primes[pi].first; int m = primes[pi].second + 1; i64 Q = 1; for (int i = 0; i < m; ++i) Q *= q;
        std::vector<LocalMat> L(maxdim + 2);
        std::vector<size_t> units(maxdim + 2, 0);
        for (int n = 1; n <= maxdim + 1; ++n) {
            std::vector<std::vector<std::pair<int32_t, i64>>> cols(E[n].size());
            for (size_t j = 0; j < E[n].size(); ++j) for (auto& kv : E[n][j]) { i64 x = kv.second % Q; if (x) cols[j].push_back({kv.first, x}); }
            L[n].init(A.C.size(n - 1), std::move(cols));
            if (n >= 2) for (size_t t = 0; t < L[n - 1].coldead.size(); ++t) if (L[n - 1].coldead[t] == 2) L[n].kill_row((int32_t)t);   // cancelled tau cells of C_{n-1}
            auto pairs = unit_cancel(L[n], q, Q);
            units[n] = pairs.size();
            for (auto& pr : pairs) { L[n].coldead[pr.second] = 2; if (n >= 2) L[n - 1].kill_col(pr.first); }
            if (stats) fprintf(stderr, "  [d_%d mod %lld^%d: %zu x %zu, %zu unit pivots, %.3f s cumulative]\n", n, (long long)q, m, A.C.size(n - 1), A.C.size(n), pairs.size(), std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        }
        for (int n = 1; n <= maxdim + 1; ++n) {
            std::vector<std::vector<std::pair<int32_t, i64>>> rest;
            for (size_t j = 0; j < L[n].cols.size(); ++j) if (!L[n].coldead[j] && !L[n].cols[j].empty()) rest.push_back(L[n].cols[j]);
            vals[n][pi] = local_elementary_valuations(rest, L[n].nrows, q, m);
            size_t rk = units[n] + vals[n][pi].size();
            if (pi == 0) rank[n] = rk; else if (rk != rank[n]) throw std::logic_error("pc_group: ranks disagree between primes");
        }
    }
    GroupHomology H;
    for (int n = 0; n <= maxdim; ++n) {
        long long fr = (long long)A.C.size(n) - (long long)rank[n] - (long long)rank[n + 1];
        if (fr < 0) throw std::logic_error("pc_group: negative rank");
        H.free.push_back((size_t)fr);
        std::vector<std::vector<u64>> pp;               // per prime: prime powers, largest first
        size_t len = 0;
        for (size_t pi = 0; pi < primes.size(); ++pi) {
            std::vector<u64> l;
            if (vals[n + 1].size() > pi) for (int v : vals[n + 1][pi]) if (v > 0 && v < primes[pi].second + 1) { u64 x = 1; for (int i = 0; i < v; ++i) x *= (u64)primes[pi].first; l.push_back(x); }
            std::sort(l.rbegin(), l.rend()); len = std::max(len, l.size()); pp.push_back(l);
        }
        std::vector<u64> inv(len, 1);
        for (auto& l : pp) for (size_t i = 0; i < l.size(); ++i) inv[i] *= l[i];
        std::reverse(inv.begin(), inv.end());
        H.torsion.push_back(inv);
    }
    if (chains_out) for (int n = 0; n <= maxdim + 1; ++n) chains_out[n] = A.C.size(n);
    if (states_out) *states_out = A.states();
    return H;
}

// ---- linear algebra over F_p for the ring --------------------------------------------------------------
typedef std::vector<std::pair<int32_t, uint32_t>> SV;   // sparse vector over F_p, sorted indices
struct Echelon {                                        // rows keyed by pivot (least index), pivot coefficient 1
    i64 p; std::map<int32_t, SV> rows;
    explicit Echelon(i64 p_) : p(p_) {}
    // reduce v by the rows; if coef is given, add the multiple of each row used (keyed by pivot)
    void reduce(SV& v, std::map<int32_t, i64>* coef = nullptr) const {
        SV tmp;
        while (!v.empty()) {
            auto it = rows.find(v.front().first);
            if (it == rows.end()) return;
            i64 f = v.front().second;
            if (coef) { i64& c = (*coef)[it->first]; c = (c + f) % p; }
            const SV& r = it->second; tmp.clear(); size_t x = 0, y = 0;
            while (x < v.size() || y < r.size()) {
                if (y == r.size() || (x < v.size() && v[x].first < r[y].first)) { tmp.push_back(v[x]); ++x; }
                else {
                    i64 s = (p - (i64)r[y].second * f % p) % p;
                    if (x < v.size() && v[x].first == r[y].first) { s = (s + v[x].second) % p; ++x; }
                    if (s) tmp.push_back({r[y].first, (uint32_t)s});
                    ++y;
                }
            }
            v.swap(tmp);
        }
    }
    // insert (after reduction); returns the new pivot or -1 when v was in the span
    int32_t insert(SV v) {
        reduce(v);
        if (v.empty()) return -1;
        i64 inv = inv_mod(v.front().second, p);
        for (auto& kv : v) kv.second = (uint32_t)((i64)kv.second * inv % p);
        int32_t piv = v.front().first; rows[piv] = std::move(v); return piv;
    }
};

// ---- the F_p splitting of F_p (x) P: Gaussian elimination across degrees --------------------------------
// Over F_p every nonzero entry of eps(d_n) is a pivot, so eliminating to the end splits the chains of degree n
// into U_n (pivot columns of d_n), L_n (pivot rows of d_{n+1}) and the survivors S_n, with d_n[L_{n-1}, U_n]
// invertible and dim H_n = |S_n|. Eliminations in different degrees commute (one deletes the rows U_n of d_{n+1},
// which miss L_n), so the strong deformation retraction of F_p (x) P onto S has the block form
//   f0(s) = s - U_n-part   (the column transform of s: d f0(s) = 0),
//   g0(v) = v_S - d_{n+1}[S, U] A^{-1} v_L,   A = d_{n+1}[L_n, U_{n+1}],
//   K0(v) = A^{-1} v_L in U_{n+1},
// with g0 f0 = 1 and dK0 + K0 d = 1 - f0 g0. The elimination is recorded as column operations: W_j = d(t_j) is the
// pivot column j at the time it was pivoted (zero on the pivot rows of earlier pivots), and t_j its transform.
struct FpDegree {
    std::vector<int32_t> prow, pcol;                  // pivots in order: (row in degree n-1, column in degree n)
    std::vector<SV> W, T;                             // W[k] over chains_{n-1}, T[k] over chains_n
    std::vector<int32_t> rowpiv;                      // chains_{n-1} -> pivot index, or -1
    std::vector<SV> trans;                            // transform of every column (pivot columns: T[k])
    std::vector<char> is_pivot_col;
};
inline FpDegree fp_eliminate(std::vector<std::vector<std::pair<int32_t, i64>>> cols, size_t nrows, const std::vector<char>& rowkilled, i64 p) {
    FpDegree R; size_t m = cols.size();
    std::vector<SV> C(m), T(m);
    std::vector<int32_t> rowcnt(nrows, 0); std::vector<std::vector<int32_t>> rowcols(nrows);
    std::vector<char> dead(m, 0), rdead(nrows, 0); std::vector<int> version(m, 0);
    for (size_t j = 0; j < m; ++j) {
        for (auto& kv : cols[j]) { i64 x = kv.second % p; if (x && !rowkilled[kv.first]) C[j].push_back({kv.first, (uint32_t)x}); }
        T[j].push_back({(int32_t)j, 1u});
        for (auto& kv : C[j]) { ++rowcnt[kv.first]; rowcols[kv.first].push_back((int32_t)j); }
    }
    typedef std::tuple<long long, int32_t, int> HE;
    std::priority_queue<HE, std::vector<HE>, std::greater<HE>> heap;
    auto cost = [&](int32_t j, int32_t* row) {
        long long best = LLONG_MAX, cs = (long long)C[j].size() - 1 + (long long)T[j].size() / 4;
        for (auto& kv : C[j]) { long long c = cs * (long long)(rowcnt[kv.first] - 1); if (c < best) { best = c; if (row) *row = kv.first; } }
        return best;
    };
    for (size_t j = 0; j < m; ++j) if (!C[j].empty()) heap.push(HE(cost((int32_t)j, nullptr), (int32_t)j, 0));
    R.rowpiv.assign(nrows, -1); R.is_pivot_col.assign(m, 0);
    auto axpy = [&](SV& y, const SV& x, i64 f, SV& tmp, bool count, int32_t col) {   // y -= f x
        tmp.clear(); size_t a = 0, b = 0;
        while (a < y.size() || b < x.size()) {
            if (b == x.size() || (a < y.size() && y[a].first < x[b].first)) { tmp.push_back(y[a]); ++a; }
            else {
                i64 s = (p - (i64)x[b].second * f % p) % p; int32_t r = x[b].first; bool had = false;
                if (a < y.size() && y[a].first == r) { s = (s + y[a].second) % p; ++a; had = true; }
                if (s) { tmp.push_back({r, (uint32_t)s}); if (count && !had) { ++rowcnt[r]; rowcols[r].push_back(col); } }
                else if (count && had) --rowcnt[r];
                ++b;
            }
        }
        y.swap(tmp);
    };
    SV tmp;
    while (!heap.empty()) {
        HE top = heap.top(); heap.pop();
        int32_t j = std::get<1>(top);
        if (dead[j] || std::get<2>(top) != version[j] || C[j].empty()) continue;
        int32_t r = -1; long long c = cost(j, &r);
        if (!heap.empty() && c > std::get<0>(heap.top())) { heap.push(HE(c, j, version[j])); continue; }
        i64 a = 0; for (auto& kv : C[j]) if (kv.first == r) a = kv.second;
        i64 ainv = inv_mod(a, p);
        std::vector<int32_t> others = rowcols[r];
        std::sort(others.begin(), others.end()); others.erase(std::unique(others.begin(), others.end()), others.end());
        for (int32_t j2 : others) {
            if (j2 == j || dead[j2]) continue;
            auto pp = std::lower_bound(C[j2].begin(), C[j2].end(), std::make_pair(r, 0u), [](const std::pair<int32_t, uint32_t>& x, const std::pair<int32_t, uint32_t>& y) { return x.first < y.first; });
            if (pp == C[j2].end() || pp->first != r) continue;
            i64 f = (i64)pp->second * ainv % p;
            axpy(C[j2], C[j], f, tmp, true, j2);
            axpy(T[j2], T[j], f, tmp, false, j2);
            ++version[j2];
            if (!C[j2].empty()) heap.push(HE(cost(j2, nullptr), j2, version[j2]));
        }
        R.rowpiv[r] = (int32_t)R.prow.size();
        R.prow.push_back(r); R.pcol.push_back(j); R.W.push_back(C[j]); R.T.push_back(T[j]); R.is_pivot_col[j] = 1;
        dead[j] = 1; for (auto& kv : C[j]) --rowcnt[kv.first];
        rdead[r] = 1; rowcols[r].clear(); rowcnt[r] = 0;
    }
    for (size_t j = 0; j < m; ++j) if (!R.is_pivot_col[j] && !C[j].empty()) throw std::logic_error("pc_group: F_p elimination left a nonzero column");
    R.trans = std::move(T);
    return R;
}
struct FpSplit {
    i64 p; int D;
    std::vector<FpDegree> E;                          // E[n] eliminates eps(d_n), n = 1..D+1
    std::vector<std::vector<int32_t>> S;              // survivors per degree, n = 0..D
    std::vector<std::vector<int32_t>> sindex;         // chains_n -> index in S[n], or -1
};
inline FpSplit fp_split(const Anick& A, int D, i64 p) {
    FpSplit F; F.p = p; F.D = D; F.E.resize(D + 2);
    std::vector<char> killed(1, 0);                  // rows of eps(d_n) deleted: U_{n-1}
    for (int n = 1; n <= D + 1; ++n) {
        F.E[n] = fp_eliminate(A.eps(n), A.C.size(n - 1), killed, p);
        killed.assign(A.C.size(n), 0);
        for (int32_t j : F.E[n].pcol) killed[j] = 1;
    }
    F.S.resize(D + 1); F.sindex.resize(D + 1);
    for (int n = 0; n <= D; ++n) {
        size_t m = A.C.size(n); F.sindex[n].assign(m, -1);
        for (size_t c = 0; c < m; ++c) {
            bool up = n >= 1 && F.E[n].is_pivot_col[c], low = F.E[n + 1].rowpiv[c] >= 0;
            if (up && low) throw std::logic_error("pc_group: a cell is paired twice");
            if (!up && !low) { F.sindex[n][c] = (int32_t)F.S[n].size(); F.S[n].push_back((int32_t)c); }
        }
    }
    return F;
}
// the cocycles g0^*(s) = e_s - w^T pi_L on the chains of degree n, for every survivor s, where w^T W_L = W[s, :]:
// back substitution over the pivots of d_{n+1} in reverse order. Only the pivots reached from s are visited:
// w at the pivot row of k is nonzero only if W[k] meets s or a row where w is already nonzero.
inline std::vector<SV> fp_cocycles(const FpSplit& F, int n, size_t nchains) {
    const FpDegree& E = F.E[n + 1]; i64 p = F.p;
    std::vector<std::vector<int32_t>> in(nchains);          // row -> pivots k whose column W[k] has it (not as pivot)
    std::vector<i64> pivval(E.prow.size(), 0);
    for (size_t k = 0; k < E.prow.size(); ++k)
        for (auto& kv : E.W[k]) { if (kv.first == E.prow[k]) pivval[k] = inv_mod(kv.second, p); else in[kv.first].push_back((int32_t)k); }
    std::vector<i64> w(nchains, 0); std::vector<char> queued(E.prow.size(), 0);
    std::vector<SV> out;
    for (int32_t s : F.S[n]) {
        std::priority_queue<int32_t> heap; std::vector<int32_t> touched, qk;
        for (int32_t k : in[s]) if (!queued[k]) { queued[k] = 1; qk.push_back(k); heap.push(k); }
        while (!heap.empty()) {
            int32_t k = heap.top(); heap.pop();
            i64 b = 0;
            for (auto& kv : E.W[k]) {
                if (kv.first == s) b = (b + kv.second) % p;
                else if (kv.first != E.prow[k] && w[kv.first]) b = (b + p - w[kv.first] * kv.second % p) % p;
            }
            if (!b) continue;
            int32_t r = E.prow[k]; w[r] = b * pivval[k] % p; touched.push_back(r);
            for (int32_t k2 : in[r]) if (k2 < k && !queued[k2]) { queued[k2] = 1; qk.push_back(k2); heap.push(k2); }
        }
        SV v; v.push_back({s, 1u});
        for (int32_t r : touched) { v.push_back({r, (uint32_t)((p - w[r]) % p)}); w[r] = 0; }
        for (int32_t k : qk) queued[k] = 0;
        std::sort(v.begin(), v.end());
        out.push_back(std::move(v));
    }
    return out;
}

inline CohomologyRing cohomology_ring(const PcGroup& G, int maxdim, int p, size_t* states_out = nullptr, bool all_products = false) {
    const bool generators_only = !all_products;
    if (p < 2) throw std::invalid_argument("pc_group: ring needs a prime");
    for (i64 d = 2; d * d <= p; ++d) if (p % d == 0) throw std::invalid_argument("pc_group: ring needs a prime");
    bool stats = std::getenv("PERSEUS_PC_STATS") != nullptr;
    auto t0 = std::chrono::steady_clock::now();
    auto lap = [&](const char* what) { if (stats) fprintf(stderr, "  [ring: %s at %.3f s]\n", what, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()); };
    Anick A(G, maxdim + 1, p);
    lap("anick resolution");
    const u64 N = (u64)G.n;
    CohomologyRing R; R.p = p; R.all_products = all_products;
    FpSplit F = fp_split(A, maxdim, p);
    lap("F_p splitting");
    std::vector<std::vector<SV>> basis(maxdim + 1), reps(maxdim + 1);   // cocycles g0^*(s); cycles f0(s)
    for (int n = 0; n <= maxdim; ++n) {
        basis[n] = fp_cocycles(F, n, A.C.size(n));
        for (int32_t c : F.S[n]) reps[n].push_back(n == 0 ? SV{{c, 1u}} : F.E[n].trans[c]);
        R.dims.push_back(F.S[n].size());
    }
    lap("cohomology bases");
    std::vector<i64> dense;
    auto coords = [&](int n, const SV& v) {             // a cocycle is determined by its values on the f0(s)
        dense.assign(A.C.size(n), 0); for (auto& kv : v) dense[kv.first] = kv.second;
        std::vector<uint32_t> x(reps[n].size(), 0);
        for (size_t i = 0; i < reps[n].size(); ++i) { i64 t = 0; for (auto& kv : reps[n][i]) t = (t + dense[kv.first] * kv.second) % p; x[i] = (uint32_t)t; }
        return x;
    };
    // Products by lifting. The decomposables in degree n are spanned by the products gamma . x with gamma a
    // minimal generator (associativity: alpha . beta = sum gamma . (x . beta) when alpha = sum gamma . x), so
    // in generators-only mode only the generators are lifted, degree by degree: the generators of degree a are
    // the basis classes outside the span of the products already found in degree a. Otherwise every basis
    // class of degree a <= maxdim/2 is lifted and all products alpha . beta, a <= b, are kept.
    std::vector<Echelon> decomp; for (int n = 0; n <= maxdim; ++n) decomp.emplace_back(p);
    // chains on which a lift is needed: those of some f0(s), and those met by the differential of a needed chain
    std::vector<std::vector<char>> needed(maxdim + 1);
    for (int n = maxdim; n >= 1; --n) {
        needed[n].assign(A.C.size(n), 0);
        for (auto& r : reps[n]) for (auto& kv : r) needed[n][kv.first] = 1;
        if (n < maxdim) for (size_t c = 0; c < A.C.size(n + 1); ++c) if (needed[n + 1][c]) for (auto& kv : A.d[n + 1][c]) needed[n][kv.first / N] = 1;
    }
    R.indecomposables.assign(maxdim + 1, 0);
    auto lift = [&](int a, size_t ai, bool all_b) {
        const SV& alpha = basis[a][ai];
        std::vector<Vec> L(A.C.size(a));                // L[c] in P_j for c of degree a + j
        for (auto& kv : alpha) L[kv.first] = Vec{{(u64)G.e, (i64)kv.second}};
        for (int j = 1; a + j <= maxdim; ++j) {
            std::vector<Vec> L2(A.C.size(a + j));
            for (size_t c = 0; c < A.C.size(a + j); ++c) {
                if (!needed[a + j][c]) continue;
                std::vector<std::pair<u64, i64>> acc;
                for (auto& kv : A.d[a + j][c]) {
                    size_t cp = kv.first / N; int g = (int)(kv.first % N);
                    for (auto& kw : L[cp]) acc.push_back({(kw.first / N) * N + (u64)G.mul((int)(kw.first % N), g), mulmod(kv.second, kw.second, p)});
                }
                Vec X; normalize(acc, p, X);
                A.hvec(j - 1, X, L2[c]);
            }
            L.swap(L2);
            int b = j;
            if (!all_b && b < a) continue;
            std::vector<char> need(A.C.size(a + b), 0);   // the product is read only on the f0(s)
            for (auto& r : reps[a + b]) for (auto& kv : r) need[kv.first] = 1;
            for (size_t bj = 0; bj < basis[b].size(); ++bj) {
                std::vector<i64> beta(A.C.size(b), 0);
                for (auto& kv : basis[b][bj]) beta[kv.first] = kv.second;
                SV prod;
                for (size_t c = 0; c < A.C.size(a + b); ++c) {
                    if (!need[c]) continue;
                    i64 s = 0; for (auto& kw : L[c]) s = (s + kw.second * beta[kw.first / N]) % p;
                    if (s) prod.push_back({(int32_t)c, (uint32_t)s});
                }
                auto x = coords(a + b, prod);
                R.products[std::make_tuple(a, (int)ai, b, (int)bj)] = x;
                SV v; for (size_t i = 0; i < x.size(); ++i) if (x[i]) v.push_back({(int32_t)i, x[i]});
                decomp[a + b].insert(v);
            }
        }
    };
    for (int a = 1; a <= maxdim; ++a) {
        R.indecomposables[a] = R.dims[a] - decomp[a].rows.size();
        for (size_t i = 0; i < R.dims[a]; ++i) {
            bool gen = !decomp[a].rows.count((int32_t)i);
            if (gen) R.generators.push_back({a, (int)i});
            if (generators_only ? (gen && a < maxdim) : (2 * a <= maxdim)) lift(a, i, generators_only);
        }
    }
    lap("products");
    R.basis = basis;
    if (states_out) *states_out = A.states();
    return R;
}

inline PcResult pc_group_compute(const PcGroup& G, int maxdim, bool homology, int ring_prime, bool all_products = false) {
    auto t0 = std::chrono::steady_clock::now();
    PcResult res;
    if (maxdim < 0) throw std::invalid_argument("pc_group: maxdim must be >= 0");
    if (homology) {
        res.chains.assign(maxdim + 2, 0);
        res.H = group_homology(G, maxdim, res.chains.data(), &res.states); res.has_homology = true;
    }
    if (ring_prime) {
        if (G.n % ring_prime) throw std::invalid_argument("pc_group: the ring prime does not divide |G| (the ring is then trivial)");
        size_t s = 0; res.R = cohomology_ring(G, maxdim, ring_prime, &s, all_products); res.has_ring = true; res.states += s;
        if (res.chains.empty()) { Chains C; C.build(G, maxdim + 1); for (int n = 0; n <= maxdim + 1; ++n) res.chains.push_back(C.size(n)); }
    }
    res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return res;
}

}  // namespace pcg
}  // namespace perseus
