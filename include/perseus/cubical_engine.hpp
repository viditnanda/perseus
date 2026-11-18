// cubical_engine.hpp -- persistent (co)homology of a cubical complex given by top-cell values, without
// building the complex (T-construction: a cell is born with the earliest top cube containing it).
//
// Cells are positions in the doubled ("Khalimsky") grid: coordinate c_i odd along the axes where the cell
// has extent. Births are stored as indices into the sorted distinct input values (32 bits per position),
// cells are ranked inside their dimension by (birth, position) with a stable radix sort, and coboundaries
// are generated on the fly. Degree 0 is union-find; degree d-1 is union-find on the dual graph; the
// degrees between are cohomology with clearing over any Ring, with apparent pairs recognised without
// building their columns: if tau is the oldest coface of sigma and sigma the youngest face of tau, no
// column processed before sigma can have pivot tau, so (sigma, tau) is a pair with pivot +-1.
#pragma once
#include "io.hpp"
#include "ring.hpp"
#include "types.hpp"
#include "cubical_morse.hpp"
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <cstring>
#include <chrono>
#include <cstdio>
#include <cstdlib>

namespace perseus {

class CubicalClassic {
public:
    typedef std::vector<uint32_t> Clear;     // positions of cleared cells in the next dimension
    bool keep_zero_length = false;
    static constexpr uint32_t ABSENT = 0xFFFFFFFFu;

    CubicalClassic(const CubicalGrid& g, int maxdim) {
        d = (int)g.ext.size();
        if (d == 0 || d > 16) throw std::runtime_error("cubical: dimension must be 1..16");
        dmax = std::min(maxdim < 0 ? d : maxdim, d);
        uint64_t total = 1;
        for (int i = 0; i < d; ++i) {
            per[i] = !g.periodic.empty() && g.periodic[i];
            L[i] = per[i] ? 2 * g.ext[i] : 2 * g.ext[i] + 1;
            stride[i] = (int64_t)total; total *= (uint64_t)L[i];
            recip[i] = (uint64_t)(((unsigned __int128)1 << 64) / (uint64_t)L[i]) + 1;   // ceil(2^64 / L): exact for pos < 2^32
            nvx[i] = per[i] ? g.ext[i] : g.ext[i] + 1; ntx[i] = g.ext[i];               // vertices / top cells along axis i
        }
        if (total >= (1ULL << 32)) throw std::runtime_error("cubical: grid too large for the implicit engine");
        N = (int64_t)total;
        auto T0 = std::chrono::steady_clock::now(); bool dbg = getenv("PERSEUS_DEBUG") != nullptr;
        auto lapt = [&](const char* w) { if (!dbg) return; auto t = std::chrono::steady_clock::now(); fprintf(stderr, "    [%s %.3f s]\n", w, std::chrono::duration<double>(t - T0).count()); T0 = t; };
        // the top cubes as indices into the sorted distinct finite values. Integer data over a modest range
        // (images) index directly by value - min; otherwise (value, cube) pairs are sorted once.
        bi.assign(N, ABSENT);
        {
            size_t nt = g.top.size();
            std::vector<int64_t> tpos(nt);
            {
                int64_t a[16] = {0};
                for (size_t t = 0; t < nt; ++t) {
                    int64_t pos = 0; for (int i = 0; i < d; ++i) pos += (2 * a[i] + 1) * stride[i];
                    tpos[t] = pos;
                    for (int i = 0; i < d; ++i) { if (++a[i] < g.ext[i]) break; a[i] = 0; }
                }
            }
            double lo = std::numeric_limits<double>::infinity(), hi = -lo; bool ints = true; size_t nf = 0;
            for (double x : g.top) if (std::isfinite(x)) { lo = std::min(lo, x); hi = std::max(hi, x); ints = ints && x == std::floor(x); ++nf; }
            if (nf && ints && hi - lo < std::max<double>(1 << 20, 4.0 * (double)nf)) {
                size_t range = (size_t)(hi - lo) + 1;
                std::vector<uint32_t> used(range, 0);
                for (double x : g.top) if (std::isfinite(x)) used[(size_t)(x - lo)] = 1;
                std::vector<uint32_t> idx(range); uint32_t m = 0;
                for (size_t r = 0; r < range; ++r) { idx[r] = m; if (used[r]) { vals.push_back(lo + (double)r); ++m; } }
                for (size_t t = 0; t < nt; ++t) { double x = g.top[t]; if (std::isfinite(x)) bi[tpos[t]] = idx[(size_t)(x - lo)]; }
            } else if (nf) {
                std::vector<std::pair<double, uint32_t>> pr; pr.reserve(nf);
                for (size_t t = 0; t < nt; ++t) if (std::isfinite(g.top[t])) pr.push_back({g.top[t], (uint32_t)t});
                std::sort(pr.begin(), pr.end());
                for (size_t i = 0; i < pr.size(); ++i) {
                    if (i == 0 || pr[i].first != pr[i - 1].first) vals.push_back(pr[i].first);
                    bi[tpos[pr[i].second]] = (uint32_t)(vals.size() - 1);
                }
            }
        }
        {
            lapt("values");
            for (int ax = 0; ax < d; ++ax) {                    // even coordinates take the min of odd neighbours
                const int64_t st = stride[ax], La = L[ax], block = st * La;
                for (int64_t outer = 0; outer < N; outer += block)
                    for (int64_t c = 0; c < La; c += 2) {
                        uint32_t* row = &bi[outer + c * st];
                        int64_t lo = c - 1, hi = c + 1;
                        if (per[ax]) { lo = (lo + La) % La; hi %= La; }
                        const uint32_t* pl = lo >= 0 ? &bi[outer + lo * st] : nullptr;
                        const uint32_t* ph = hi < La ? &bi[outer + hi * st] : nullptr;
                        if (pl && ph) for (int64_t in = 0; in < st; ++in) row[in] = std::min(pl[in], ph[in]);
                        else if (pl) for (int64_t in = 0; in < st; ++in) row[in] = pl[in];
                        else for (int64_t in = 0; in < st; ++in) row[in] = ph[in];
                    }
            }
        }
        lapt("propagate");
        // cells by dimension: counts for all; positions sorted by (birth, position) and ranks only for the
        // dimensions 1..d-1 (columns and rows of the reductions). Vertices and top cells are handled by
        // union-finds indexed by grid position, comparing (birth, position) directly.
        cells.assign(d + 1, {}); cnt_.assign(d + 1, 0);
        {
            int64_t c[16] = {0}; int k = 0;
            for (int64_t pos = 0; pos < N; ++pos) {
                if (bi[pos] != ABSENT) ++cnt_[k];
                for (int i = 0; i < d; ++i) { k -= (int)(c[i] & 1); if (++c[i] < L[i]) { k += (int)(c[i] & 1); break; } c[i] = 0; }
            }
            lapt("count");
            std::vector<std::vector<uint32_t>> key(d + 1);
            std::vector<size_t> fill(d + 1, 0);
            for (int kk = 0; kk <= d; ++kk) if (sorted_dim(kk)) { cells[kk].resize(cnt_[kk]); key[kk].resize(cnt_[kk]); }
            for (int i = 0; i < d; ++i) c[i] = 0;
            k = 0;
            for (int64_t pos = 0; pos < N; ++pos) {
                if (bi[pos] != ABSENT && sorted_dim(k)) { size_t t = fill[k]++; cells[k][t] = (uint32_t)pos; key[k][t] = bi[pos]; }
                for (int i = 0; i < d; ++i) { k -= (int)(c[i] & 1); if (++c[i] < L[i]) { k += (int)(c[i] & 1); break; } c[i] = 0; }
            }
            lapt("bucket");
            if (d >= 3) rk.assign(N, -1);                  // ranks are only needed for rows: dimensions 2..d-1
            std::vector<uint32_t> k2, v2;
            for (int kk = 0; kk <= d; ++kk) {
                if (!sorted_dim(kk)) continue;
                radix_sort(key[kk], cells[kk], k2, v2);
                if (kk >= 2) for (size_t i = 0; i < cells[kk].size(); ++i) rk[cells[kk][i]] = (int32_t)i;
                std::vector<uint32_t>().swap(key[kk]);
            }
            lapt("sort and rank");
        }
    }
    // (birth, position) keys of the vertices / top cells by compact grid index (no divisions)
    void compact_keys(bool topcells, std::vector<uint64_t>& out) const {
        const int64_t* n = topcells ? ntx : nvx;
        size_t m = 1; for (int i = 0; i < d; ++i) m *= (size_t)n[i];
        out.resize(m);
        int64_t c[16] = {0}; int64_t pos = topcells ? 0 : 0;
        if (topcells) for (int i = 0; i < d; ++i) pos += stride[i];
        for (size_t t = 0; t < m; ++t) {
            out[t] = bi[pos] == ABSENT ? ~0ULL : ((uint64_t)bi[pos] << 32 | (uint64_t)pos);
            for (int i = 0; i < d; ++i) { pos += 2 * stride[i]; if (++c[i] < n[i]) break; pos -= 2 * stride[i] * n[i]; c[i] = 0; }
        }
    }
    int maxdim() const { return dmax; }
    size_t count(int k) const { return k >= 0 && k <= d ? cnt_[k] : 0; }

    template <class R>
    BRun<Clear> run(const R& ring, int kmin, int kmax, const Clear* cleared_in) const { return run(ring, kmin, kmax, cleared_in, keep_zero_length); }
    template <class R>
    BRun<Clear> run(const R& ring, int kmin, int kmax, const Clear* cleared_in, bool keep) const {
        typedef typename R::T T;
        struct Ent { int32_t row; T c; };
        BRun<Clear> out; out.kmin = kmin;
        kmax = std::min(kmax, dmax);
        int nd = std::max(0, kmax + 1);
        out.bars.assign(nd, {}); out.clear_out.assign(nd, {});
        for (int l = 0; l < R::lanes; ++l) out.moduli.push_back(ring.modulus(l));
        const T ONE = ring.one(), MONE = ring.neg(ring.one());
        const double INF = std::numeric_limits<double>::infinity();
        Clear cleared; if (cleared_in) cleared = *cleared_in;
        std::vector<uint64_t> cbits;                       // cleared positions of the current degree
        auto setbits = [&](const Clear& cl) { cbits.assign((size_t)(N + 63) / 64, 0); for (uint32_t q : cl) cbits[q >> 6] |= 1ULL << (q & 63); };
        auto isclr = [&](uint64_t q) { return (cbits[q >> 6] >> (q & 63)) & 1; };
        int kstart = kmin;
        if (kmin == 0 && kmax >= 0) {                     // H0 by union-find (elder rule), ring independent
            // vertices by compact grid index (spatially local); a component is represented by its oldest
            // vertex, compared by (birth, position)
            size_t nv = 1; for (int i = 0; i < d; ++i) nv *= (size_t)nvx[i];
            std::vector<int32_t> par(nv); for (size_t i = 0; i < nv; ++i) par[i] = (int32_t)i;
            auto find = [&](int32_t x) { while (par[x] != x) x = par[x] = par[par[x]]; return x; };
            std::vector<uint64_t> VK; compact_keys(false, VK);
            auto vkey = [&](int32_t v) { return VK[v]; };
            int64_t c[16];
            for (size_t e = 0; e < cells[1].size(); ++e) {
                int64_t pos = cells[1][e]; decode(pos, c);
                int ax = 0; while (!(c[ax] & 1)) ++ax;
                int64_t lo = c[ax] - 1, hi = c[ax] + 1; if (per[ax] && hi == L[ax]) hi = 0;
                if (lo == hi) continue;                   // boundary zero (periodic axis of extent 1)
                int64_t vl = 0, vh = 0, m = 1;
                for (int i = 0; i < d; ++i) { int64_t ci = i == ax ? lo : c[i]; int64_t ch2 = i == ax ? hi : c[i]; vl += (ci / 2) * m; vh += (ch2 / 2) * m; m *= nvx[i]; }
                int32_t a = find((int32_t)vl), bb = find((int32_t)vh);
                if (a == bb) continue;
                uint64_t ka = vkey(a), kb = vkey(bb);
                int32_t young = ka > kb ? a : bb, old = ka > kb ? bb : a;
                par[young] = old;
                double bv = vals[(uint32_t)(std::max(ka, kb) >> 32)], dv = val(pos);
                if (dv > bv || keep) out.bars[0].push_back({bv, dv});
                if (dmax > 0) out.clear_out[0].push_back((uint32_t)pos);
                ++out.columns;
            }
            for (size_t i = 0; i < nv; ++i) {
                if (find((int32_t)i) != (int32_t)i) continue;
                uint64_t kv = vkey((int32_t)i); if ((uint32_t)(kv >> 32) == ABSENT) continue;   // not in the complex
                out.bars[0].push_back({vals[(uint32_t)(kv >> 32)], INF});
            }
            std::sort(out.bars[0].begin(), out.bars[0].end());
            cleared = out.clear_out[0];
            kstart = 1;
        }
        std::vector<Ent> w, tmp, cob, heap;
        auto cmp = [](const Ent& x, const Ent& y) { return x.row > y.row; };   // min-heap on the row
        size_t dbg_maxw = 0; long long dbg_merged = 0;
        for (int k = kstart; k <= kmax; ++k) {
            size_t nk = count(k), nk1 = count(k + 1);
            setbits(cleared);
            Clear next_cleared;
            bool top = (k == dmax);
            if (k == d) {                                 // top cells: essential unless cleared
                size_t nt = 1; for (int i = 0; i < d; ++i) nt *= (size_t)ntx[i];
                for (size_t t = 0; t < nt; ++t) {
                    int64_t pos = top_pos((int64_t)t);
                    if (bi[pos] == ABSENT || isclr((uint64_t)pos)) continue;
                    ++out.columns; out.bars[k].push_back({val(pos), INF});
                }
                std::sort(out.bars[k].begin(), out.bars[k].end());
                continue;
            }
            if (k == d - 1 && k >= 1) {
                // dual union-find: every (d-1)-cell has at most two cofaces with opposite signs, so delta^(d-1)
                // is the incidence matrix of the dual graph (d-cells, plus a virtual node V for cofaces outside
                // the complex). Columns in decreasing order with pivot = smallest row is H0 in the reversed
                // filtration: a component is represented by its largest rank (V largest of all); joining two
                // components pairs the column with the smaller representative. All pivots are +-1.
                size_t nt = 1; for (int i = 0; i < d; ++i) nt *= (size_t)ntx[i];
                const int32_t V = (int32_t)nt;            // top cells by compact grid index, V virtual
                std::vector<int32_t> par(nt + 1); for (size_t i = 0; i <= nt; ++i) par[i] = (int32_t)i;
                auto find = [&](int32_t x) { while (par[x] != x) x = par[x] = par[par[x]]; return x; };
                std::vector<uint64_t> TK; compact_keys(true, TK);
                auto tkey = [&](int32_t t) -> uint64_t { return t == V ? ~0ULL : TK[t]; };
                int64_t c[16];
                for (size_t jj = nk; jj-- > 0;) {
                    int32_t j = (int32_t)jj;
                    int64_t pos = cells[k][j];
                    if (isclr((uint64_t)pos)) continue;
                    ++out.columns;
                    decode(pos, c);
                    int i = 0; while (c[i] & 1) ++i;               // the one even axis
                    int64_t q[2]; int nq = 0;
                    int64_t lo = c[i] - 1, hi = c[i] + 1;
                    if (per[i] && lo < 0) lo = L[i] - 1;
                    if (lo >= 0) q[nq++] = lo;
                    if (hi < L[i]) q[nq++] = hi;
                    if (nq == 2 && q[0] == q[1]) nq = 0;           // the same coface from both sides: zero column
                    int32_t ends[2]; int ne = 0;
                    for (int t = 0; t < nq; ++t) {
                        int64_t qp = pos + (q[t] - c[i]) * stride[i];
                        if (bi[qp] == ABSENT) { ends[ne++] = V; continue; }
                        int64_t ti = 0, m = 1;
                        for (int a2 = 0; a2 < d; ++a2) { int64_t ca = a2 == i ? q[t] : c[a2]; ti += ((ca - 1) / 2) * m; m *= ntx[a2]; }
                        ends[ne++] = (int32_t)ti;
                    }
                    if (ne == 1) ends[ne++] = V;                   // a face on the boundary of the grid
                    double bj = val(pos);
                    int32_t ra = ne ? find(ends[0]) : 0, rb = ne ? find(ends[1]) : 0;
                    if (ne == 0 || ra == rb) { out.bars[k].push_back({bj, INF}); continue; }
                    // representatives are the youngest of their component (largest key, V largest of all)
                    uint64_t ka = tkey(ra), kb = tkey(rb);
                    int32_t lo2 = ka < kb ? ra : rb, hi2 = ka < kb ? rb : ra;
                    par[lo2] = hi2;
                    int64_t lp = (int64_t)(uint32_t)TK[lo2];
                    double dv = val(lp);
                    if (dv > bj || keep) out.bars[k].push_back({bj, dv});
                    if (!top) next_cleared.push_back((uint32_t)lp);
                }
                std::sort(out.bars[k].begin(), out.bars[k].end());
                out.clear_out[k] = next_cleared;
                cleared.swap(next_cleared);
                continue;
            }
            out.note_bound(k, (double)std::min(nk, nk1) * 0.5 * std::log2(2.0 * (k + 1)), (double)std::min(nk, nk1));
            std::vector<int32_t> piv_col(nk1, -1);
            std::vector<T> piv_inv(nk1);
            std::unordered_map<int32_t, std::vector<Ent>> stored;
            // coboundary of the k-cell of rank j: entries (row rank, coefficient), sorted, duplicates merged
            auto coboundary = [&](int32_t j, std::vector<Ent>& o) {
                int64_t pos = cells[k][j], c[16]; decode(pos, c);
                Ent buf[32]; int n = 0, odd = 0;
                for (int i = 0; i < d; ++i) {
                    if (c[i] & 1) { ++odd; continue; }
                    bool neg = odd & 1;                   // sign (-1)^j, j = odd axes before i
                    int64_t lo = c[i] - 1, hi = c[i] + 1;
                    if (per[i] && lo < 0) lo = L[i] - 1;
                    if (lo >= 0) { int64_t q = pos + (lo - c[i]) * stride[i]; if (bi[q] != ABSENT) buf[n++] = {rk[q], neg ? MONE : ONE}; }
                    if (hi < L[i]) { int64_t q = pos + (hi - c[i]) * stride[i]; if (bi[q] != ABSENT) buf[n++] = {rk[q], neg ? ONE : MONE}; }
                }
                for (int a = 1; a < n; ++a) { Ent x = buf[a]; int b2 = a; while (b2 > 0 && buf[b2 - 1].row > x.row) { buf[b2] = buf[b2 - 1]; --b2; } buf[b2] = x; }
                o.clear();
                for (int a = 0; a < n;) {
                    int e = a; T s = buf[a].c; while (++e < n && buf[e].row == buf[a].row) s = ring.add(s, buf[e].c);
                    if (!ring.is_zero(s)) o.push_back({buf[a].row, s});
                    a = e;
                }
            };
            // apparent pair: the oldest coface tau of sigma has sigma as its youngest face; returns tau's
            // position and the coefficient of tau in delta(sigma), or -1
            auto apparent = [&](int64_t pos, const int64_t* c, T& coef) -> int64_t {
                uint64_t best = ~0ULL; int64_t tau = -1; int tax = -1; bool tneg = false; int odd = 0;
                for (int i = 0; i < d; ++i) {
                    if (c[i] & 1) { ++odd; continue; }
                    int64_t lo = c[i] - 1, hi = c[i] + 1;
                    if (per[i] && lo < 0) lo = L[i] - 1;
                    for (int s2 = 0; s2 < 2; ++s2) {
                        int64_t cc = s2 ? hi : lo;
                        if (cc < 0 || cc >= L[i]) continue;
                        int64_t q = pos + (cc - c[i]) * stride[i];
                        if (bi[q] == ABSENT) continue;
                        uint64_t key = (uint64_t)bi[q] << 32 | (uint64_t)q;
                        if (key < best) { best = key; tau = q; tax = i; tneg = ((odd & 1) != 0) != (s2 == 1); }
                    }
                }
                if (tau < 0) return -1;
                if (per[tax] && L[tax] == 2) return -1;   // tau meets sigma from both sides: use the general path
                // youngest face of tau
                int64_t ct[16]; std::memcpy(ct, c, sizeof(int64_t) * d); ct[tax] = (tau - pos) / stride[tax] + c[tax];
                uint64_t worst = 0; int64_t yf = -1;
                for (int i = 0; i < d; ++i) {
                    if (!(ct[i] & 1)) continue;
                    int64_t lo = ct[i] - 1, hi = ct[i] + 1;
                    if (per[i] && hi == L[i]) hi = 0;
                    for (int64_t cc : {lo, hi}) {
                        int64_t f = tau + (cc - ct[i]) * stride[i];
                        uint64_t key = (uint64_t)bi[f] << 32 | (uint64_t)f;
                        if (key >= worst) { worst = key; yf = f; }
                    }
                }
                if (yf != pos) return -1;
                coef = tneg ? MONE : ONE;
                return tau;
            };
            int64_t cbuf[16];
            for (size_t jj = nk; jj-- > 0;) {
                int32_t j = (int32_t)jj;
                int64_t pos = cells[k][j];
                if (isclr((uint64_t)pos)) continue;
                ++out.columns;
                double bj = val(pos);
                {
                    decode(pos, cbuf);
                    T coef; int64_t tau = apparent(pos, cbuf, coef);
                    if (tau >= 0) {
                        int32_t r = rk[tau];
                        piv_col[r] = j; piv_inv[r] = coef;     // +-1 is its own inverse
                        double dv = val(tau);
                        if (dv > bj || keep) out.bars[k].push_back({bj, dv});
                        if (!top) next_cleared.push_back((uint32_t)tau);
                        continue;
                    }
                }
                // working column: a sorted vector while it is short (an addition merges), a min-heap on the row
                // with lazy cancellation once it is long (an addition pushes the few entries of the added column
                // and never copies the working column)
                coboundary(j, w);
                bool ch = false, found = false, heapmode = false; Ent piv{0, ring.zero()};
                heap.clear();
                while (true) {
                    found = false;
                    if (!heapmode) { if (!w.empty()) { piv = w[0]; found = true; } }
                    else while (!heap.empty()) {                    // pop the smallest row, summing duplicates
                        std::pop_heap(heap.begin(), heap.end(), cmp); Ent e = heap.back(); heap.pop_back();
                        while (!heap.empty() && heap.front().row == e.row) { std::pop_heap(heap.begin(), heap.end(), cmp); e.c = ring.add(e.c, heap.back().c); heap.pop_back(); }
                        if (!ring.is_zero(e.c)) { piv = e; found = true; break; }
                    }
                    if (!found) break;
                    int32_t i = piv_col[piv.row];
                    if (i < 0) break;
                    T f = ring.neg(ring.mul(piv.c, piv_inv[piv.row]));
                    const std::vector<Ent>* ci;
                    auto it = stored.find(i);
                    if (it != stored.end()) ci = &it->second; else { coboundary(i, cob); ci = &cob; }
                    if (!heapmode && w.size() > 2048) {               // switch: w[0] cancels against (*ci)[0]
                        heapmode = true; heap.assign(w.begin() + 1, w.end()); std::make_heap(heap.begin(), heap.end(), cmp);
                    }
                    if (heapmode) {
                        for (size_t q = 1; q < ci->size(); ++q) {   // entry 0 is the pivot: it cancels piv
                            T sc = ring.mul(f, (*ci)[q].c);
                            if (!ring.is_zero(sc)) { heap.push_back({(*ci)[q].row, sc}); std::push_heap(heap.begin(), heap.end(), cmp); }
                        }
                    } else {
                        tmp.clear(); size_t a2 = 1, bq = 1, be = ci->size();
                        while (a2 < w.size() && bq < be) {
                            if (w[a2].row < (*ci)[bq].row) tmp.push_back(w[a2++]);
                            else if ((*ci)[bq].row < w[a2].row) { T sc = ring.mul(f, (*ci)[bq].c); if (!ring.is_zero(sc)) tmp.push_back({(*ci)[bq].row, sc}); ++bq; }
                            else { T sc = ring.add(w[a2].c, ring.mul(f, (*ci)[bq].c)); if (!ring.is_zero(sc)) tmp.push_back({w[a2].row, sc}); ++a2; ++bq; }
                        }
                        for (; a2 < w.size(); ++a2) tmp.push_back(w[a2]);
                        for (; bq < be; ++bq) { T sc = ring.mul(f, (*ci)[bq].c); if (!ring.is_zero(sc)) tmp.push_back({(*ci)[bq].row, sc}); }
                        w.swap(tmp);
                    }
                    ch = true; ++out.additions; dbg_merged += (long long)ci->size(); dbg_maxw = std::max(dbg_maxw, heapmode ? heap.size() : w.size());
                }
                if (!found) { out.bars[k].push_back({bj, INF}); continue; }
                T pc = piv.c; int32_t r = piv.row;
                if (ring.partly_zero(pc)) out.unlucky = true;
                piv_col[r] = j; piv_inv[r] = ring.inv(pc);
                double dv = val(cells[k + 1][r]);
                if (dv > bj || keep) out.bars[k].push_back({bj, dv});
                if (!top) next_cleared.push_back(cells[k + 1][r]);
                if (!ring.is_pm1(pc)) { Pivot pv; pv.k = k; pv.key = (uint64_t)j; pv.birth = bj; pv.death = dv; ring.store(pc, pv); out.nonunit.push_back(std::move(pv)); }
                if (ch && !heapmode) stored[j] = w;
                else if (ch) {                                     // materialise the reduced column: pivot, then the rest
                    std::vector<Ent>& RC = stored[j]; RC.clear(); RC.push_back(piv);
                    std::sort(heap.begin(), heap.end(), [](const Ent& x, const Ent& y) { return x.row < y.row; });
                    for (size_t q = 0; q < heap.size();) {
                        size_t e = q; T sum = heap[q].c; while (++e < heap.size() && heap[e].row == heap[q].row) sum = ring.add(sum, heap[e].c);
                        if (!ring.is_zero(sum)) RC.push_back({heap[q].row, sum});
                        q = e;
                    }
                }
            }
            if (getenv("PERSEUS_DEBUG")) fprintf(stderr, "    [degree %d: %lld columns, %lld additions so far, max working length %zu, total merged entries %lld]\n", k, out.columns, out.additions, dbg_maxw, dbg_merged);
            std::sort(out.bars[k].begin(), out.bars[k].end());
            out.clear_out[k] = next_cleared;
            cleared.swap(next_cleared);
        }
        return out;
    }

private:
    int d = 0, dmax = 0;
    int64_t N = 0;
    int64_t L[16], stride[16];
    uint64_t recip[16];
    uint8_t per[16];
    int64_t nvx[16], ntx[16];
    std::vector<size_t> cnt_;
    bool sorted_dim(int k) const { return (k >= 1 && k <= d - 1) || (k == 1 && d == 1); }
    // position of a vertex / top cell given by its compact grid index
    int64_t vertex_pos(int64_t v) const { int64_t p = 0; for (int i = 0; i < d; ++i) { p += 2 * (v % nvx[i]) * stride[i]; v /= nvx[i]; } return p; }
    int64_t top_pos(int64_t t) const { int64_t p = 0; for (int i = 0; i < d; ++i) { p += (2 * (t % ntx[i]) + 1) * stride[i]; t /= ntx[i]; } return p; }

    std::vector<double> vals;           // sorted distinct input values
    std::vector<uint32_t> bi;           // per position: index into vals, or ABSENT
    std::vector<int32_t> rk;
    std::vector<std::vector<uint32_t>> cells;
    double val(int64_t pos) const { return vals[bi[pos]]; }
    // coordinates of a position: one chain of exact reciprocal divisions
    void decode(int64_t pos, int64_t* c) const {
        uint64_t q = (uint64_t)pos;
        for (int i = 0; i < d; ++i) {
            uint64_t qn = (uint64_t)(((unsigned __int128)q * recip[i]) >> 64);
            c[i] = (int64_t)(q - qn * (uint64_t)L[i]); q = qn;
        }
    }
    static void radix_sort(std::vector<uint32_t>& key, std::vector<uint32_t>& val, std::vector<uint32_t>& k2, std::vector<uint32_t>& v2) {
        size_t n = key.size(); if (n < 2) return;
        k2.resize(n); v2.resize(n);
        uint32_t allor = 0, alland = ~0u; for (size_t i = 0; i < n; ++i) { allor |= key[i]; alland &= key[i]; }
        for (int shift = 0; shift < 32; shift += 8) {
            if ((((allor ^ alland) >> shift) & 0xFF) == 0) continue;   // this digit is the same everywhere
            size_t cnt[257] = {0};
            for (size_t i = 0; i < n; ++i) cnt[((key[i] >> shift) & 0xFF) + 1]++;
            for (int i = 0; i < 256; ++i) cnt[i + 1] += cnt[i];
            for (size_t i = 0; i < n; ++i) { size_t t = cnt[(key[i] >> shift) & 0xFF]++; k2[t] = key[i]; v2[t] = val[i]; }
            key.swap(k2); val.swap(v2);
        }
    }
};

// The cubical backend: the discrete-gradient engine (cubical_morse.hpp) where it applies (dimensions 2..4,
// no periodic axes), the classic engine otherwise. engine = "auto", "morse" or "classic"; the environment
// variable PERSEUS_CUBICAL overrides "auto".
class CubicalBackend {
public:
    typedef std::vector<uint32_t> Clear;
    bool keep_zero_length = false;
    CubicalBackend(const CubicalGrid& g, int maxdim, std::string engine = "auto") {
        if (engine == "auto" || engine == "implicit") { const char* e = getenv("PERSEUS_CUBICAL"); engine = e && *e ? e : "auto"; }
        if (engine != "auto" && engine != "morse" && engine != "classic") throw std::runtime_error("cubical engine must be auto, morse or classic");
        if (engine == "morse" && !CubicalMorse::supports(g)) throw std::runtime_error("cubical: the morse engine needs dimension 2..4 and no periodic axes");
        if (engine != "classic" && CubicalMorse::supports(g)) {
            try { m.reset(new CubicalMorse(g, maxdim)); return; }
            catch (std::overflow_error&) { if (engine == "morse") throw; }    // a path sum beyond int64: classic engine
        }
        c.reset(new CubicalClassic(g, maxdim));
    }
    int maxdim() const { return m ? m->maxdim() : c->maxdim(); }
    const char* name() const { return m ? "morse" : "classic"; }
    const CubicalMorse* morse() const { return m.get(); }
    template <class R>
    BRun<Clear> run(const R& ring, int kmin, int kmax, const Clear* in) const {
        return m ? m->run(ring, kmin, kmax, in, keep_zero_length) : c->run(ring, kmin, kmax, in, keep_zero_length);
    }
private:
    std::unique_ptr<CubicalClassic> c;
    std::unique_ptr<CubicalMorse> m;
};

} // namespace perseus
