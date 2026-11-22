// collapse.hpp -- iterated edge collapse of a Vietoris-Rips (flag) filtration.
//
// Boissonnat-Pritam (SoCG 2020) edge collapse in the backward "shift and trim" form of Glisse-Pritam
// (SoCG 2022). Edges are processed from the last to the first; every step is a valid operation on the
// current filtration, so the result has persistence modules isomorphic (over Z) to the input's:
//   * an edge e = uv is dominated at time tau if some w in the common neighbourhood C_tau(u, v) is
//     adjacent at time tau to every other vertex of C_tau(u, v);
//   * a dominated edge may be shifted to the next time at which it stops being dominated, and an edge
//     dominated until the end of the filtration may be removed.
//
// Implementation. During the backward sweep the graph at the current time is kept as bitsets B(v), and
// the final graph as bitsets F(v). For the edge e = uv at time t:
//   C = B(u) & B(v);  a dominator is w in C with C - w contained in B(w)            (word-parallel)
//   future common neighbours X = F(u) & F(v) - C;  if X is contained in B(w) as well, w dominates e at
//   every later time (its neighbourhood only grows), so e is removed at once           (word-parallel)
// Only otherwise are the future common neighbours visited in time order: domination can be lost only
// when a new common neighbour x arrives, and the current dominator survives iff wx is already present.
#pragma once
#include "io.hpp"
#include <vector>
#include <algorithm>
#include <cstdint>
#include <chrono>
#include <limits>
#include <cstdio>
#include <cstdlib>

namespace perseus {

struct CollapseStats { size_t edges_before = 0, edges_after = 0, shifted = 0; int passes = 0; double seconds = 0; bool used = true; };

// Returns a RipsInput whose lower triangle has collapsed edges removed (+inf) and shifted edges moved to
// later filtration values. Only edges <= R.threshold are considered (others are absent). Passes stop when
// one removes fewer than min_gain * (remaining edges).
// If worth_it > 0: give up (st.used = false, nothing returned) when the first pass removes less than
// worth_it of the edges -- shifting without removing destroys the cone levels the Rips engine exploits.
// The first pass is abandoned early if, after its first third, the removal rate is below worth_it / 4.
inline RipsInput edge_collapse(const RipsInput& R, int max_passes, CollapseStats& st, double min_gain = 0.0, double worth_it = 0.0, std::string* trace = nullptr) {
    typedef uint64_t W;
    auto t0 = std::chrono::steady_clock::now();
    const int n = R.n, WD = std::max(1, (n + 63) / 64);
    const int32_t NONE = std::numeric_limits<int32_t>::max();
    struct E { double d; int a, b; };
    std::vector<E> edges;
    for (int a = 1; a < n; ++a) for (int b = 0; b < a; ++b) {
        double d = R.dist(a, b);
        if (std::isfinite(d) && d <= R.threshold) edges.push_back({d, a, b});
    }
    std::sort(edges.begin(), edges.end(), [](const E& x, const E& y) { return x.d != y.d ? x.d < y.d : (x.a != y.a ? x.a < y.a : x.b < y.b); });
    const size_t m = edges.size();
    st.edges_before = m;
    std::vector<double> val(m); for (size_t r = 0; r < m; ++r) val[r] = edges[r].d;
    std::vector<int32_t> T((size_t)n * n, NONE);                 // current insertion rank of each edge
    auto Tm = [&](int x, int y) -> int32_t& { return T[(size_t)x * n + y]; };
    for (size_t r = 0; r < m; ++r) Tm(edges[r].a, edges[r].b) = Tm(edges[r].b, edges[r].a) = (int32_t)r;
    std::vector<int32_t> cur(m); for (size_t r = 0; r < m; ++r) cur[r] = (int32_t)r;
    std::vector<W> B((size_t)n * WD), F((size_t)n * WD), Cb(WD), X(WD);
    auto setb = [&](std::vector<W>& S, int x, int y, bool on) {
        W mx = W(1) << (y & 63), my = W(1) << (x & 63);
        if (on) { S[(size_t)x * WD + (y >> 6)] |= mx; S[(size_t)y * WD + (x >> 6)] |= my; }
        else { S[(size_t)x * WD + (y >> 6)] &= ~mx; S[(size_t)y * WD + (x >> 6)] &= ~my; }
    };
    // does w dominate the set S at the current time (S - w contained in B(w))?
    auto dominates = [&](int w, const std::vector<W>& S) {
        const W* bw = &B[(size_t)w * WD];
        for (int k = 0; k < WD; ++k) {
            W rest = S[k] & ~bw[k];
            if (k == (w >> 6)) rest &= ~(W(1) << (w & 63));
            if (rest) return false;
        }
        return true;
    };
    std::vector<std::pair<int32_t, int>> ev;
    std::vector<int> Clist;
    size_t alive_before = m;
    int hint = -1;
    for (int pass = 0; pass < max_passes; ++pass) {
        size_t changed = 0;
        std::vector<int32_t> order;
        for (size_t r = 0; r < m; ++r) if (cur[r] != NONE) order.push_back((int32_t)r);
        std::sort(order.begin(), order.end(), [&](int32_t x, int32_t y) { return cur[x] != cur[y] ? cur[x] > cur[y] : x > y; });
        std::fill(B.begin(), B.end(), 0); std::fill(F.begin(), F.end(), 0);
        for (int32_t r : order) { setb(B, edges[r].a, edges[r].b, true); setb(F, edges[r].a, edges[r].b, true); }
        for (size_t oi = 0; oi < order.size(); ++oi) {
            int32_t r = order[oi];
            int u = edges[r].a, v = edges[r].b;
            int32_t t = cur[r];
            const W* bu = &B[(size_t)u * WD]; const W* bv = &B[(size_t)v * WD];
            bool any = false;
            for (int k = 0; k < WD; ++k) { Cb[k] = bu[k] & bv[k]; any |= Cb[k] != 0; }
            int w = -1;
            if (any) {
                if (hint >= 0 && ((Cb[hint >> 6] >> (hint & 63)) & 1) && dominates(hint, Cb)) w = hint;
                for (int k = WD; k-- > 0 && w < 0;) {
                    W bits = Cb[k];
                    while (bits && w < 0) { int x = k * 64 + 63 - __builtin_clzll(bits); bits &= ~(W(1) << (x & 63)); if (dominates(x, Cb)) w = x; }
                }
            }
            if (w >= 0) {
                hint = w;
                const W* fu = &F[(size_t)u * WD]; const W* fv = &F[(size_t)v * WD]; const W* bw = &B[(size_t)w * WD];
                bool covered = true;
                for (int k = 0; k < WD; ++k) {
                    X[k] = fu[k] & fv[k] & ~Cb[k];
                    if (X[k] & ~bw[k]) covered = false;
                }
                int32_t newt = NONE;
                if (!covered) {
                    ev.clear();
                    for (int k = 0; k < WD; ++k) { W bits = X[k]; while (bits) { int x = k * 64 + __builtin_ctzll(bits); bits &= bits - 1; ev.push_back({std::max(Tm(u, x), Tm(v, x)), x}); } }
                    std::sort(ev.begin(), ev.end());
                    Clist.clear();
                    for (int k = 0; k < WD; ++k) { W bits = Cb[k]; while (bits) { Clist.push_back(k * 64 + __builtin_ctzll(bits)); bits &= bits - 1; } }
                    auto ok = [&](int c, int32_t tau) { for (int y : Clist) if (y != c && Tm(c, y) > tau) return false; return true; };
                    for (size_t i = 0; i < ev.size();) {
                        int32_t tau = ev[i].first; bool broke = false; size_t j = i;
                        for (; j < ev.size() && ev[j].first == tau; ++j) { Clist.push_back(ev[j].second); if (Tm(w, ev[j].second) > tau) broke = true; }
                        if (broke) {
                            int w2 = -1;
                            if (ok(ev[j - 1].second, tau)) w2 = ev[j - 1].second;
                            for (size_t q = Clist.size(); q-- > 0 && w2 < 0;) if (ok(Clist[q], tau)) w2 = Clist[q];
                            if (w2 < 0) { newt = tau; break; }
                            w = w2;
                        }
                        i = j;
                    }
                }
                if (newt == NONE) { Tm(u, v) = Tm(v, u) = NONE; cur[r] = NONE; setb(F, u, v, false); ++changed; }
                else if (newt != t) { Tm(u, v) = Tm(v, u) = newt; cur[r] = newt; ++changed; }
            }
            setb(B, u, v, false);           // the sweep moves below time t
            if (pass == 0 && worth_it > 0 && (oi + 1 == order.size() / 10 || oi + 1 == order.size() / 3)) {
                size_t removed = 0; for (size_t q = 0; q <= oi; ++q) removed += cur[order[q]] == NONE;
                if ((double)removed < 0.25 * worth_it * (double)(oi + 1)) { st.used = false; break; }
            }
        }
        if (!st.used) { st.passes = 1; break; }
        st.passes = pass + 1;
        size_t alive_now = 0; for (size_t r = 0; r < m; ++r) alive_now += cur[r] != NONE;
        if (getenv("PERSEUS_COLLAPSE_TRACE")) {
            char line[192];
            snprintf(line, sizeof(line), "    pass %d: %zu changed, %zu edges left, %.3f s\n", pass + 1, changed, alive_now,
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
            if (trace) trace->append(line);
            else fputs(line, stderr);
        }
        size_t removed_pass = alive_before - alive_now; alive_before = alive_now;
        if (pass == 0 && worth_it > 0 && (double)removed_pass < worth_it * (double)m) { st.used = false; break; }
        if (!changed || (double)removed_pass < min_gain * (double)alive_now) break;
    }
    if (!st.used) { st.edges_after = m; st.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); return RipsInput(); }
    RipsInput out; out.n = n; out.maxdim = R.maxdim; out.threshold = R.threshold;
    out.lower.assign((size_t)n * (n - 1) / 2, std::numeric_limits<double>::infinity());
    size_t kept = 0;
    for (size_t r = 0; r < m; ++r) if (cur[r] != NONE) {
        int a = edges[r].a, b = edges[r].b;
        out.lower[(size_t)a * (a - 1) / 2 + b] = val[cur[r]];
        ++kept; if (cur[r] != (int32_t)r) ++st.shifted;
    }
    st.edges_after = kept;
    st.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return out;
}

} // namespace perseus
