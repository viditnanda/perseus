// rips_engine.hpp -- Vietoris-Rips persistence without building the complex.
//
// Two engines over any Ring:
//  * cone-level cohomology (degrees >= 1): Ripser-style lazy cohomology with clearing, where the
//    filtration is refined by levels (rank of the longest edge). A level whose link L_e has an
//    apex v adjacent to all of L_e - v is a cone level (Boissonnat-Pritam edge domination); in the
//    apex-first order (level, flag, colex) its cells pair up as apparent pairs (beta, beta u v),
//    are never enumerated as columns, and their rows never become pivots of base cells
//    (projection lemma). Apparent pairs of non-cone levels are recognised implicitly.
//  * active H1 enumeration (homology direction): used when only H0, H1 are asked for.
#pragma once
#include "ring.hpp"
#include "types.hpp"
#include "io.hpp"
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace perseus {

inline double enclosing_radius(const RipsInput& R) {
    double enc = std::numeric_limits<double>::infinity();
    for (int i = 0; i < R.n; ++i) {
        double mx = 0;
        for (int j = 0; j < R.n; ++j) if (i != j) mx = std::max(mx, R.dist(i, j));
        enc = std::min(enc, mx);
    }
    return R.n <= 1 ? 0.0 : enc;
}

class RipsBackend {
public:
    typedef std::vector<uint64_t> Clear;     // sorted colex indices of cleared cells
    typedef uint64_t W;
    typedef uint64_t idx_t;
    bool keep_zero_length = false;           // only bars computed explicitly; apparent/cone pairs are never seen
    bool use_h1_engine = true;

    RipsBackend(const RipsInput& R, int maxdim) : dmax(maxdim) {
        N = R.n; WORDS = std::max(1, (N + 63) / 64);
        if (maxdim + 2 > 31) throw std::runtime_error("Rips: maxdim too large");
        double thr = R.threshold;
        BIN.assign((size_t)(N + 1) * (dmax + 4), 0);
        for (int i = 0; i <= N; ++i) { bin(i, 0) = 1; for (int j = 1; j < dmax + 4; ++j) bin(i, j) = i ? bin(i - 1, j - 1) + bin(i - 1, j) : 0; }
        // overflow guard for colex indices of (dmax+2)-subsets
        long double est = 1; for (int j = 0; j < dmax + 2; ++j) est = est * (N - j) / (j + 1);
        if (est > 1.8e19L) throw std::runtime_error("Rips: too many simplices to index in 64 bits; lower maxdim or n");
        struct E { double d; idx_t ci; int a, b; };
        std::vector<E> edges;
        for (int a = 1; a < N; ++a) for (int b = 0; b < a; ++b) {
            double d = R.dist(a, b);
            if (std::isfinite(d) && d <= thr) edges.push_back({d, bin(b, 1) + bin(a, 2), a, b});
        }
        std::sort(edges.begin(), edges.end(), [](const E& x, const E& y) { return x.d != y.d ? x.d < y.d : x.ci < y.ci; });
        NE = edges.size();
        RANK.assign((size_t)N * (N - 1) / 2 + 1, -1);
        EDIAM.resize(NE); EA.resize(NE); EB.resize(NE); APEX.assign(NE, -1); NB.assign((size_t)N * WORDS, 0);
        for (size_t r = 0; r < NE; ++r) {
            auto& e = edges[r];
            RANK[(size_t)e.a * (e.a - 1) / 2 + e.b] = (int32_t)r; EDIAM[r] = e.d; EA[r] = e.a; EB[r] = e.b;
            NB[(size_t)e.a * WORDS + (e.b >> 6)] |= W(1) << (e.b & 63); NB[(size_t)e.b * WORDS + (e.a >> 6)] |= W(1) << (e.a & 63);
        }
        // cone detection
        std::vector<W> G((size_t)N * WORDS, 0), L(WORDS);
        for (size_t r = 0; r < NE; ++r) {
            int a = EA[r], b = EB[r];
            bool any = false;
            for (int t = 0; t < WORDS; ++t) { L[t] = G[(size_t)a * WORDS + t] & G[(size_t)b * WORDS + t]; any |= L[t] != 0; }
            int apex = -1;
            if (any) for (int t = 0; t < WORDS && apex < 0; ++t) {
                W bb = L[t];
                while (bb && apex < 0) {
                    int bp = __builtin_ctzll(bb); bb &= bb - 1; int v = t * 64 + bp; bool ok = true;
                    for (int s2 = 0; s2 < WORDS && ok; ++s2) { W rest = L[s2]; if (s2 == t) rest &= ~(W(1) << bp); ok = (rest & ~G[(size_t)v * WORDS + s2]) == 0; }
                    if (ok) apex = v;
                }
            }
            APEX[r] = apex; if (apex >= 0) ++ncone;
            G[(size_t)a * WORDS + (b >> 6)] |= W(1) << (b & 63); G[(size_t)b * WORDS + (a >> 6)] |= W(1) << (a & 63);
        }
        // H0 by union-find (ring independent); cone-level edges never merge components
        std::vector<int> par(N); for (int i = 0; i < N; ++i) par[i] = i;
        auto find = [&](int x) { while (par[x] != x) x = par[x] = par[par[x]]; return x; };
        for (size_t r = 0; r < NE; ++r) {
            if (APEX[r] >= 0) continue;
            int u = find(EA[r]), w = find(EB[r]);
            if (u != w) {
                par[std::max(u, w)] = std::min(u, w);
                int32_t vv[2] = {std::min(EA[r], EB[r]), std::max(EA[r], EB[r])};
                tree.push_back(comb(vv, 2));
                if (EDIAM[r] > 0 || keep_zero_length) bars0.push_back({0.0, EDIAM[r]});
            }
        }
        for (int i = 0; i < N; ++i) if (find(i) == i) bars0.push_back({0.0, std::numeric_limits<double>::infinity()});
        std::sort(bars0.begin(), bars0.end()); std::sort(tree.begin(), tree.end());
    }
    int maxdim() const { return dmax; }
    size_t num_edges() const { return NE; }
    long long cone_levels() const { return ncone; }

    template <class R>
    BRun<Clear> run(const R& ring, int kmin, int kmax, const Clear* cleared_in) const {
        BRun<Clear> out; out.kmin = kmin;
        kmax = std::min(kmax, dmax);
        out.bars.assign(kmax + 1, {}); out.clear_out.assign(kmax + 1, {});
        for (int l = 0; l < R::lanes; ++l) out.moduli.push_back(ring.modulus(l));
        Clear cleared;
        if (kmin == 0) { out.bars[0] = bars0; out.clear_out[0] = tree; cleared = tree; }
        else if (cleared_in) cleared = *cleared_in;
        for (int k = std::max(1, kmin); k <= kmax; ++k) {
            if (k == 1 && dmax == 1 && use_h1_engine) { h1_degree(ring, out); continue; }
            cone_degree(ring, k, k == dmax, cleared, out);
            cleared = out.clear_out[k];
        }
        return out;
    }

private:
    int N = 0, WORDS = 1, dmax = 1;
    size_t NE = 0;
    long long ncone = 0;
    std::vector<idx_t> BIN;
    std::vector<int32_t> RANK, EA, EB, APEX;
    std::vector<double> EDIAM;
    std::vector<W> NB;
    std::vector<Bar> bars0;
    Clear tree;

    inline idx_t& bin(int i, int j) { return BIN[(size_t)i * (dmax + 4) + j]; }
    inline idx_t binc(int i, int j) const { return BIN[(size_t)i * (dmax + 4) + j]; }
    inline int rk(int a, int b) const { if (a < b) std::swap(a, b); return RANK[(size_t)a * (a - 1) / 2 + b]; }
    inline idx_t comb(const int32_t* v, int m) const { idx_t s = 0; for (int i = 0; i < m; ++i) s += binc(v[i], i + 1); return s; }
    void decode(idx_t idx, int m, int32_t* v) const {
        int top = N;
        for (int i = m; i >= 1; --i) {
            int lo = i - 1, hi = top - 1;
            while (lo < hi) { int mid = (lo + hi + 1) / 2; if (binc(mid, i) <= idx) lo = mid; else hi = mid - 1; }
            v[i - 1] = lo; idx -= binc(lo, i); top = lo;
        }
    }
    inline uint32_t level_of(const int32_t* v, int m) const {
        uint32_t L = 0; for (int i = 0; i < m; ++i) for (int j = i + 1; j < m; ++j) { uint32_t r = (uint32_t)rk(v[i], v[j]); if (r > L) L = r; }
        return L;
    }
    // common neighbourhood of the vertices of x into cand
    inline void common_nb(const int32_t* xv, int m, W* cand) const {
        for (int t = 0; t < WORDS; ++t) cand[t] = NB[(size_t)xv[0] * WORDS + t];
        for (int i = 1; i < m; ++i) for (int t = 0; t < WORDS; ++t) cand[t] &= NB[(size_t)xv[i] * WORDS + t];
    }

    template <class T> struct Ent { uint64_t lk; idx_t idx; T c; };
    template <class T> struct HeapCmp {
        bool operator()(const Ent<T>& a, const Ent<T>& b) const { return a.lk != b.lk ? a.lk > b.lk : a.idx > b.idx; }
    };

    template <class R, class Heap>
    void push_cobd(const R& ring, const int32_t* xv, int m, uint32_t xl, typename R::T mult, Heap& h, W* cand) const {
        common_nb(xv, m, cand);
        int32_t tv[32];
        typename R::T nm = ring.neg(mult);
        for (int t = 0; t < WORDS; ++t) {
            W b = cand[t];
            while (b) {
                int w = t * 64 + __builtin_ctzll(b); b &= b - 1;
                uint32_t L = xl;
                for (int i = 0; i < m; ++i) { uint32_t r = (uint32_t)rk(xv[i], w); if (r > L) L = r; }
                int ap = APEX[L], flag = 1;
                if (ap >= 0) {
                    bool has = (ap == w); for (int i = 0; i < m && !has; ++i) has = xv[i] == ap;
                    if (!has) continue;                       // base cell: row deleted (projection lemma)
                    flag = 0;
                }
                int pos = 0; while (pos < m && xv[pos] < w) ++pos;
                for (int i = 0; i < pos; ++i) tv[i] = xv[i];
                tv[pos] = w; for (int i = pos; i < m; ++i) tv[i + 1] = xv[i];
                h.push({(uint64_t)L * 2 + flag, comb(tv, m + 1), (pos & 1) ? nm : mult});
            }
        }
    }
    // first same-level cofacet of a cell in a non-cone level (increasing w = increasing colex index)
    bool first_same_level(const int32_t* xv, int m, uint32_t xl, idx_t& oidx, int& opos, W* cand) const {
        common_nb(xv, m, cand);
        int32_t tv[32];
        for (int t = 0; t < WORDS; ++t) {
            W b = cand[t];
            while (b) {
                int w = t * 64 + __builtin_ctzll(b); b &= b - 1;
                bool same = true;
                for (int i = 0; i < m && same; ++i) if ((uint32_t)rk(xv[i], w) > xl) same = false;
                if (!same) continue;
                int pos = 0; while (pos < m && xv[pos] < w) ++pos;
                for (int i = 0; i < pos; ++i) tv[i] = xv[i];
                tv[pos] = w; for (int i = pos; i < m; ++i) tv[i + 1] = xv[i];
                oidx = comb(tv, m + 1); opos = pos;
                return true;
            }
        }
        return false;
    }
    // youngest same-level facet: drop the smallest vertex outside the level edge
    bool youngest_facet(const int32_t* tv, int m, uint32_t L, int32_t* fv, int& dropped) const {
        int a = EA[L], b = EB[L];
        for (int i = 0; i < m; ++i) if (tv[i] != a && tv[i] != b) {
            int j = 0; for (int q = 0; q < m; ++q) if (q != i) fv[j++] = tv[q];
            dropped = i; return true;
        }
        return false;
    }
    bool apparent_upper(const int32_t* tv, int m, uint32_t L, idx_t tidx, int32_t* fv, int& dpos, W* cand) const {
        if (!youngest_facet(tv, m, L, fv, dpos)) return false;
        idx_t e; int p; if (!first_same_level(fv, m - 1, L, e, p, cand)) return false;
        return e == tidx;
    }
    // stream the k-cells (k+1 vertices) of non-cone levels, level by level
    template <class F>
    void stream_cells(int k, F&& emit) const {
        std::vector<W> G((size_t)N * WORDS, 0), scratch((size_t)(k + 2) * WORDS);
        std::vector<int32_t> cur;
        for (size_t r = 0; r < NE; ++r) {
            int a = EA[r], b = EB[r];
            if (APEX[r] < 0) {
                if (k == 1) { int32_t vv[2] = {std::min(a, b), std::max(a, b)}; emit((uint32_t)r, vv); }
                else {
                    for (int t = 0; t < WORDS; ++t) scratch[t] = G[(size_t)a * WORDS + t] & G[(size_t)b * WORDS + t];
                    auto rec = [&](auto&& self, const W* cand, int sz) -> void {
                        if (sz == k - 1) {
                            int32_t vv[32]; int mm = 0; vv[mm++] = a; vv[mm++] = b; for (int x : cur) vv[mm++] = x;
                            std::sort(vv, vv + mm); emit((uint32_t)r, vv); return;
                        }
                        W* nxt = &scratch[(size_t)(sz + 1) * WORDS];
                        for (int t = 0; t < WORDS; ++t) {
                            W bb = cand[t];
                            while (bb) {
                                int bp = __builtin_ctzll(bb); bb &= bb - 1; int x = t * 64 + bp;
                                const W* gx = &G[(size_t)x * WORDS];
                                for (int s2 = 0; s2 < WORDS; ++s2) nxt[s2] = cand[s2] & gx[s2];
                                int cw = x >> 6; for (int s2 = cw + 1; s2 < WORDS; ++s2) nxt[s2] = 0;
                                int rem = x & 63; nxt[cw] &= rem ? ((W(1) << rem) - 1) : 0;
                                cur.push_back(x); self(self, nxt, sz + 1); cur.pop_back();
                            }
                        }
                    };
                    rec(rec, &scratch[0], 0);
                }
            }
            G[(size_t)a * WORDS + (b >> 6)] |= W(1) << (b & 63); G[(size_t)b * WORDS + (a >> 6)] |= W(1) << (a & 63);
        }
    }

    template <class R>
    void cone_degree(const R& ring, int k, bool top, const Clear& cleared, BRun<Clear>& out) const {
        typedef typename R::T T;
        typedef Ent<T> E;
        typedef std::priority_queue<E, std::vector<E>, HeapCmp<T>> Heap;
        const T ONE = ring.one(), MONE = ring.neg(ring.one());
        auto pm = [&](int s) { return s ? MONE : ONE; };
        int mk = k + 1;
        std::vector<W> cand(WORDS);
        struct Col { uint32_t lvl; idx_t idx; };
        std::vector<Col> cols;
        double streamed = 0;
        stream_cells(k, [&](uint32_t lvl, const int32_t* cv) {
            streamed += 1;
            idx_t ci = comb(cv, mk);
            if (std::binary_search(cleared.begin(), cleared.end(), ci)) return;
            int32_t fv[32], tv[32]; int dp;
            if (k >= 2 && apparent_upper(cv, mk, lvl, ci, fv, dp, cand.data())) return;
            idx_t e0; int p0;
            if (first_same_level(cv, mk, lvl, e0, p0, cand.data())) {
                decode(e0, mk + 1, tv); youngest_facet(tv, mk + 1, lvl, fv, dp);
                if (comb(fv, mk) == ci) return;                   // apparent lower: paired, never a column
            }
            cols.push_back({lvl, ci});
        });
        std::sort(cols.begin(), cols.end(), [](const Col& x, const Col& y) { return x.lvl != y.lvl ? x.lvl > y.lvl : x.idx > y.idx; });
        out.columns += (long long)cols.size();
        // the matrix reduced here has at most `streamed` columns, and every row (a (k+1)-cell) has at most k+2
        // entries +-1: every minor is at most (k+2)^(s/2) in absolute value, s <= streamed
        out.note_bound(k, streamed * 0.5 * std::log2((double)k + 2), streamed);
        std::unordered_map<idx_t, std::pair<int, T>> pivot_of;           // pivot cell -> (column, inverse pivot coeff)
        std::unordered_map<int, std::vector<std::pair<idx_t, T>>> V;
        Clear next_cleared;
        const double INF = std::numeric_limits<double>::infinity();
        for (size_t j = 0; j < cols.size(); ++j) {
            double bd = EDIAM[cols[j].lvl];
            int32_t cv[32]; decode(cols[j].idx, mk, cv);
            {   // emergent pair: the first same-level cofacet is not yet a pivot and not apparent
                idx_t e0; int p0;
                if (first_same_level(cv, mk, cols[j].lvl, e0, p0, cand.data())) {
                    int32_t tv[32], fv[32]; int dp; decode(e0, mk + 1, tv);
                    if (!pivot_of.count(e0) && !apparent_upper(tv, mk + 1, cols[j].lvl, e0, fv, dp, cand.data())) {
                        pivot_of[e0] = {(int)j, pm(p0 & 1)};          // coefficient +-1 is its own inverse
                        if (!top) next_cleared.push_back(e0);
                        continue;
                    }
                }
            }
            Heap h; std::vector<std::pair<idx_t, T>> wv;
            push_cobd(ring, cv, mk, cols[j].lvl, ONE, h, cand.data());
            auto add_cell = [&](idx_t cidx, T fac) {
                int32_t v[32]; decode(cidx, mk, v);
                push_cobd(ring, v, mk, level_of(v, mk), fac, h, cand.data()); wv.push_back({cidx, fac});
            };
            E piv; bool found = false;
            while (true) {
                found = false;
                while (!h.empty()) {
                    E t = h.top(); h.pop();
                    T c = t.c;
                    while (!h.empty() && h.top().idx == t.idx && h.top().lk == t.lk) { c = ring.add(c, h.top().c); h.pop(); }
                    if (!ring.is_zero(c)) { t.c = c; piv = t; found = true; break; }
                }
                if (!found) break;
                auto it = pivot_of.find(piv.idx);
                if (it != pivot_of.end()) {
                    int i = it->second.first; T factor = ring.neg(ring.mul(piv.c, it->second.second));
                    h.push(piv); ++out.additions;
                    add_cell(cols[i].idx, factor);
                    auto vt = V.find(i);
                    if (vt != V.end()) for (auto& pr : vt->second) add_cell(pr.first, ring.mul(factor, pr.second));
                    continue;
                }
                uint32_t pl = (uint32_t)(piv.lk >> 1); int flag = (int)(piv.lk & 1);
                int32_t tv[32]; decode(piv.idx, mk + 1, tv);
                if (flag == 0) {                                       // cone cell: add the coboundary of its partner
                    int ap = APEX[pl], apos = 0; while (tv[apos] != ap) ++apos;
                    int32_t bv[32]; int jj = 0; for (int i = 0; i <= mk; ++i) if (i != apos) bv[jj++] = tv[i];
                    T factor = ring.neg(ring.mul(piv.c, pm(apos & 1)));
                    h.push(piv);
                    push_cobd(ring, bv, mk, pl, factor, h, cand.data()); wv.push_back({comb(bv, mk), factor});
                    continue;
                }
                int32_t fv[32]; int dp;
                if (apparent_upper(tv, mk + 1, pl, piv.idx, fv, dp, cand.data())) {
                    T factor = ring.neg(ring.mul(piv.c, pm(dp & 1)));
                    h.push(piv);
                    push_cobd(ring, fv, mk, pl, factor, h, cand.data()); wv.push_back({comb(fv, mk), factor});
                    continue;
                }
                break;
            }
            if (found) {
                double dd = EDIAM[piv.lk >> 1];
                if (ring.partly_zero(piv.c)) out.unlucky = true;
                if (!ring.is_pm1(piv.c)) {
                    Pivot pv; pv.k = k; pv.key = cols[j].idx; pv.birth = bd; pv.death = dd; ring.store(piv.c, pv); out.nonunit.push_back(std::move(pv));
                }
                pivot_of[piv.idx] = {(int)j, ring.inv(piv.c)};
                if (!top) next_cleared.push_back(piv.idx);
                std::sort(wv.begin(), wv.end(), [](const std::pair<idx_t, T>& x, const std::pair<idx_t, T>& y) { return x.first < y.first; });
                std::vector<std::pair<idx_t, T>> mv;
                for (size_t a = 0; a < wv.size();) {
                    size_t b = a; T c = ring.zero();
                    while (b < wv.size() && wv[b].first == wv[a].first) c = ring.add(c, wv[b++].second);
                    if (!ring.is_zero(c)) mv.push_back({wv[a].first, c});
                    a = b;
                }
                if (!mv.empty()) V[(int)j] = std::move(mv);
                if (dd > bd) out.bars[k].push_back({bd, dd});
            } else out.bars[k].push_back({bd, INF});
        }
        std::sort(out.bars[k].begin(), out.bars[k].end());
        std::sort(next_cleared.begin(), next_cleared.end());
        out.clear_out[k] = std::move(next_cleared);
    }

    // ---- H1, homology direction, active enumeration ------------------------------------
    template <class R>
    void h1_degree(const R& ring, BRun<Clear>& out) const {
        typedef typename R::T T;
        typedef std::vector<std::pair<int, T>> SVec;
        const T ONE = ring.one(), MONE = ring.neg(ring.one());
        auto axpy = [&](const SVec& x, T f, const SVec& y, SVec& o) {
            o.clear(); size_t i = 0, j = 0;
            while (i < x.size() || j < y.size()) {
                if (j == y.size() || (i < x.size() && x[i].first < y[j].first)) { o.push_back(x[i++]); continue; }
                if (i == x.size() || y[j].first < x[i].first) { T c = ring.mul(f, y[j].second); if (!ring.is_zero(c)) o.push_back({y[j].first, c}); ++j; continue; }
                T c = ring.add(x[i].second, ring.mul(f, y[j].second));
                if (!ring.is_zero(c)) o.push_back({x[i].first, c}); ++i; ++j;
            }
        };
        auto lb = [](SVec& v, int key) { return std::lower_bound(v.begin(), v.end(), key, [](const std::pair<int, T>& e, int kk) { return e.first < kk; }); };
        std::vector<W> G((size_t)N * WORDS, 0), ACT((size_t)N * WORDS, 0);
        std::vector<int> par(N); for (int i = 0; i < N; ++i) par[i] = i;
        auto find = [&](int x) { while (par[x] != x) x = par[x] = par[par[x]]; return x; };
        std::unordered_map<int, SVec> Rm; std::unordered_map<int, std::vector<int>> Ri;
        auto setact = [&](int a, int b, bool on) {
            W ma = W(1) << (b & 63), mb = W(1) << (a & 63);
            if (on) { ACT[(size_t)a * WORDS + (b >> 6)] |= ma; ACT[(size_t)b * WORDS + (a >> 6)] |= mb; }
            else { ACT[(size_t)a * WORDS + (b >> 6)] &= ~ma; ACT[(size_t)b * WORDS + (a >> 6)] &= ~mb; }
        };
        SVec bar, tmp;
        auto rof = [&](int e) -> const SVec* { auto it = Rm.find(e); return it == Rm.end() ? nullptr : &it->second; };
        auto process = [&](int x0, int x1, int x2, int yr) {
            int e12 = rk(x1, x2), e02 = rk(x0, x2), e01 = rk(x0, x1);
            bar.clear();
            const SVec* r12 = rof(e12); const SVec* r02 = rof(e02); const SVec* r01 = rof(e01);
            if (r12) { axpy(bar, ONE, *r12, tmp); bar.swap(tmp); }
            if (r02) { axpy(bar, MONE, *r02, tmp); bar.swap(tmp); }
            if (r01) { axpy(bar, ONE, *r01, tmp); bar.swap(tmp); }
            if (bar.empty()) return;
            int g = bar.back().first; T cg = bar.back().second;
            ++out.additions;
            if (ring.partly_zero(cg)) out.unlucky = true;
            if (!ring.is_pm1(cg)) {
                Pivot pv; pv.k = 1; pv.key = (uint64_t)g; pv.birth = EDIAM[g]; pv.death = EDIAM[yr]; ring.store(cg, pv); out.nonunit.push_back(std::move(pv));
            }
            if (EDIAM[g] < EDIAM[yr]) out.bars[1].push_back({EDIAM[g], EDIAM[yr]});
            T fac0 = ring.neg(ring.inv(cg));
            std::vector<int> zs; zs.swap(Ri[g]); Ri.erase(g);
            for (int z : zs) {
                auto it = Rm.find(z); if (it == Rm.end()) continue;
                SVec& rz = it->second;
                auto pos = lb(rz, g);
                if (pos == rz.end() || pos->first != g) continue;
                T fac = ring.mul(fac0, pos->second);
                axpy(rz, fac, bar, tmp);
                for (auto& pr : tmp) { auto q = lb(rz, pr.first); if (q == rz.end() || q->first != pr.first) Ri[pr.first].push_back(z); }
                rz.swap(tmp);
                if (rz.empty()) { Rm.erase(it); setact(EA[z], EB[z], false); }
            }
        };
        std::vector<W> L(WORDS), C(WORDS);
        for (size_t r = 0; r < NE; ++r) {
            int a = EA[r], b = EB[r];
            bool any = false;
            for (int t = 0; t < WORDS; ++t) { L[t] = G[(size_t)a * WORDS + t] & G[(size_t)b * WORDS + t]; any |= L[t] != 0; }
            int u = find(a), v = find(b);
            bool born = (u == v);
            if (!born) par[std::max(u, v)] = std::min(u, v);
            int first_w = -1;
            if (born && any) for (int t = 0; t < WORDS && first_w < 0; ++t) { W x = L[t] & ~ACT[(size_t)a * WORDS + t] & ~ACT[(size_t)b * WORDS + t]; if (x) first_w = t * 64 + __builtin_ctzll(x); }
            if (born && first_w >= 0) { /* zero-length class, pivot +-1, never stored */ }
            else if (born) { Rm[(int)r] = SVec{{(int)r, ONE}}; Ri[(int)r].push_back((int)r); setact(a, b, true); }
            if (any) {
                bool yact = Rm.count((int)r) > 0;
                const W* src = yact ? L.data() : C.data();
                if (!yact) for (int t = 0; t < WORDS; ++t) C[t] = L[t] & (ACT[(size_t)a * WORDS + t] | ACT[(size_t)b * WORDS + t]);
                for (int t = 0; t < WORDS; ++t) {
                    W x = src[t];
                    while (x) {
                        int w = t * 64 + __builtin_ctzll(x); x &= x - 1;
                        if (w == first_w) continue;
                        if (!Rm.count((int)r)) {
                            bool aw = (ACT[(size_t)a * WORDS + (w >> 6)] >> (w & 63)) & 1, bw = (ACT[(size_t)b * WORDS + (w >> 6)] >> (w & 63)) & 1;
                            if (!aw && !bw) continue;
                        }
                        int xs[3] = {a, b, w}; std::sort(xs, xs + 3); process(xs[0], xs[1], xs[2], (int)r);
                    }
                }
            }
            G[(size_t)a * WORDS + (b >> 6)] |= W(1) << (b & 63); G[(size_t)b * WORDS + (a >> 6)] |= W(1) << (a & 63);
        }
        const double INF = std::numeric_limits<double>::infinity();
        // d_2 has at most NE rows and columns with 3 entries +-1
        out.note_bound(1, (double)NE * 0.5 * std::log2(3.0), (double)NE);
        for (auto& kv : Ri) out.bars[1].push_back({EDIAM[kv.first], INF});
        std::sort(out.bars[1].begin(), out.bars[1].end());
    }
};

} // namespace perseus
