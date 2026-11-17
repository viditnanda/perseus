// explicit_engine.hpp -- persistent cohomology with clearing for an explicit filtered cell
// complex, over any Ring from ring.hpp (F_2, F_p, or Q carried through several large primes).
//
// Degree k: columns are the k-cells (not cleared) in decreasing filtration order, each holding its
// coboundary; the pivot of a column is its earliest (k+1)-cell. A column reducing to zero is an
// essential class; the pivot rows of degree k are cleared from the columns of degree k+1.
// Reduced columns are stored only when the reduction actually changed them.
#pragma once
#include "complex.hpp"
#include "ring.hpp"
#include "types.hpp"
#include <vector>
#include <cstdint>
#include <algorithm>
#include <numeric>
#include <map>
#include <climits>
#include <cmath>
#include <cstdint>

namespace perseus {

// complex sorted into degrees, with coboundaries by local index
struct Prepared {
    int D = -1;                                   // top dimension
    std::vector<std::vector<int32_t>> cells;      // per degree: global ids in filtration order
    std::vector<std::vector<double>> birth;       // per degree: births in that order
    // coboundary of local k-cell s: rows cob_row[k][cob_off[k][s] .. cob_off[k][s+1]) (local (k+1)-indices, ascending)
    std::vector<std::vector<int64_t>> cob_off;
    std::vector<std::vector<int32_t>> cob_row;
    std::vector<std::vector<int64_t>> cob_coef;

    explicit Prepared(const FilteredComplex& C) {
        size_t n = C.size();
        D = C.max_dim();
        struct Key { double b; int32_t d, id; };
        std::vector<Key> keys(n);
        for (size_t i = 0; i < n; ++i) keys[i] = {C.birth[i], (int32_t)C.dim[i], (int32_t)i};
        std::sort(keys.begin(), keys.end(), [](const Key& x, const Key& y) {
            if (x.b != y.b) return x.b < y.b;
            if (x.d != y.d) return x.d < y.d;
            return x.id < y.id;
        });
        std::vector<int32_t> ord(n);
        for (size_t i = 0; i < n; ++i) ord[i] = keys[i].id;
        std::vector<Key>().swap(keys);
        cells.assign(D + 1, {}); birth.assign(D + 1, {});
        std::vector<int32_t> local(n);
        for (int32_t g : ord) { int k = C.dim[g]; local[g] = (int32_t)cells[k].size(); cells[k].push_back(g); birth[k].push_back(C.birth[g]); }
        cob_off.assign(D + 1, {}); cob_row.assign(D + 1, {}); cob_coef.assign(D + 1, {});
        for (int k = 0; k <= D; ++k) {
            size_t nk = cells[k].size();
            std::vector<int64_t> cnt(nk + 1, 0);
            if (k < D) for (int32_t t : cells[k + 1]) for (int64_t e = C.off[t]; e < C.off[t + 1]; ++e) cnt[local[C.face[e]] + 1]++;
            for (size_t s = 0; s < nk; ++s) cnt[s + 1] += cnt[s];
            cob_off[k] = cnt;
            cob_row[k].resize(cnt[nk]); cob_coef[k].resize(cnt[nk]);
            if (k < D) {
                std::vector<int64_t> fill(cnt.begin(), cnt.end() - 1);
                for (size_t tl = 0; tl < cells[k + 1].size(); ++tl) {           // ascending rows -> sorted coboundaries
                    int32_t t = cells[k + 1][tl];
                    for (int64_t e = C.off[t]; e < C.off[t + 1]; ++e) {
                        int64_t p = fill[local[C.face[e]]]++;
                        cob_row[k][p] = (int32_t)tl; cob_coef[k][p] = C.coef[e];
                    }
                }
                // a face listed twice in one boundary gives adjacent equal rows: merge them
                bool dup = false;
                for (size_t s = 0; s < nk && !dup; ++s) for (int64_t e = cnt[s] + 1; e < cnt[s + 1]; ++e) if (cob_row[k][e] == cob_row[k][e - 1]) { dup = true; break; }
                if (dup) {
                    std::vector<int64_t> noff(nk + 1, 0); std::vector<int32_t> nr; std::vector<int64_t> nc;
                    for (size_t s = 0; s < nk; ++s) {
                        for (int64_t e = cnt[s]; e < cnt[s + 1];) {
                            int64_t f = e, v = 0; while (f < cnt[s + 1] && cob_row[k][f] == cob_row[k][e]) v += cob_coef[k][f++];
                            if (v != 0) { nr.push_back(cob_row[k][e]); nc.push_back(v); }
                            e = f;
                        }
                        noff[s + 1] = (int64_t)nr.size();
                    }
                    cob_off[k] = noff; cob_row[k] = nr; cob_coef[k] = nc;
                }
            }
        }
        // a priori bound data per degree k (pairs between k- and (k+1)-cells): every minor of d_{k+1} has
        // log2|det| <= hb[k] (Hadamard, by columns or by rows), and there are at most pb[k] pairs
        hb.assign(D + 1, 0.0); pb.assign(D + 1, 0.0);
        for (int k = 0; k < D; ++k) {
            std::vector<double> col2(cells[k + 1].size(), 0.0);
            double rows = 0;
            for (size_t s2 = 0; s2 < cells[k].size(); ++s2) {
                double r2 = 0;
                for (int64_t e = cob_off[k][s2]; e < cob_off[k][s2 + 1]; ++e) { double c = (double)cob_coef[k][e]; r2 += c * c; col2[cob_row[k][e]] += c * c; }
                if (r2 > 1) rows += 0.5 * std::log2(r2);
            }
            double cols = 0; for (double c2 : col2) if (c2 > 1) cols += 0.5 * std::log2(c2);
            hb[k] = std::min(rows, cols);
            pb[k] = (double)std::min(cells[k].size(), cells[k + 1].size());
        }
        // degree 0 by union-find is valid when every 1-cell has boundary 0 or +-(v - w)
        graph_like = true;
        if (D >= 1) {
            ends.assign(cells[1].size(), {-1, -1});
            for (size_t el = 0; el < cells[1].size() && graph_like; ++el) {
                int32_t t = cells[1][el];
                int64_t nb = C.off[t + 1] - C.off[t];
                if (nb == 0) continue;
                if (nb == 2) {                       // the common case, without a map
                    int64_t e0 = C.off[t];
                    int32_t f0 = local[C.face[e0]], f1 = local[C.face[e0 + 1]];
                    int64_t c0 = C.coef[e0], c1 = C.coef[e0 + 1];
                    if (f0 != f1 && c0 + c1 == 0 && (c0 == 1 || c0 == -1)) { ends[el] = {f0, f1}; continue; }
                }
                std::map<int32_t, int64_t> acc;
                for (int64_t e = C.off[t]; e < C.off[t + 1]; ++e) acc[local[C.face[e]]] += C.coef[e];
                std::vector<std::pair<int32_t, int64_t>> nz; for (auto& kv : acc) if (kv.second) nz.push_back(kv);
                if (nz.empty()) continue;
                if (nz.size() == 2 && nz[0].second + nz[1].second == 0 && (nz[0].second == 1 || nz[0].second == -1)) ends[el] = {nz[0].first, nz[1].first};
                else graph_like = false;
            }
        }
    }
    bool graph_like = false;
    std::vector<double> hb, pb;
    std::vector<std::pair<int32_t, int32_t>> ends;   // local endpoints of 1-cells (-1 if boundary 0)
    size_t count(int k) const { return (k >= 0 && k <= D) ? cells[k].size() : 0; }

    // boundaries by local index (k-cell -> (k-1)-cells, ascending), built on demand by transposing
    mutable std::vector<std::vector<int64_t>> bd_off;
    mutable std::vector<std::vector<int32_t>> bd_row;
    mutable std::vector<std::vector<int64_t>> bd_coef;
    void build_boundaries() const {
        if (!bd_off.empty()) return;
        bd_off.assign(D + 1, {}); bd_row.assign(D + 1, {}); bd_coef.assign(D + 1, {});
        for (int k = 1; k <= D; ++k) {
            size_t nk = cells[k].size(), nk1 = cells[k - 1].size();
            std::vector<int64_t> cnt(nk + 1, 0);
            for (int64_t e = 0; e < (int64_t)cob_row[k - 1].size(); ++e) cnt[cob_row[k - 1][e] + 1]++;
            for (size_t t = 0; t < nk; ++t) cnt[t + 1] += cnt[t];
            bd_off[k] = cnt; bd_row[k].resize(cnt[nk]); bd_coef[k].resize(cnt[nk]);
            std::vector<int64_t> fill(cnt.begin(), cnt.end() - 1);
            for (size_t s = 0; s < nk1; ++s)
                for (int64_t e = cob_off[k - 1][s]; e < cob_off[k - 1][s + 1]; ++e) {
                    int64_t q = fill[cob_row[k - 1][e]]++;
                    bd_row[k][q] = (int32_t)s; bd_coef[k][q] = cob_coef[k - 1][e];
                }
        }
    }
};

// Backend for certify(): degrees kmin..kmax; Clear = the k-cells to skip (pivot rows of degree k-1).
struct ExplicitBackend {
    typedef std::vector<int32_t> Clear;
    const Prepared& P;
    int dmax;
    bool keep_zero_length = false;
    enum Direction { COHOMOLOGY = 0, HOMOLOGY = 1, FORWARD = 2 };
    int direction = COHOMOLOGY;

    ExplicitBackend(const Prepared& p, int maxdim, int dir = COHOMOLOGY) : P(p), dmax(std::min(maxdim, p.D)), direction(dir) {
        if (direction == HOMOLOGY || direction == FORWARD) P.build_boundaries();
    }
    int maxdim() const { return dmax; }

    template <class R>
    BRun<Clear> run(const R& ring, int kmin, int kmax, const Clear* cleared_in) const {
        BRun<Clear> out;
        if (direction == HOMOLOGY) out = run_homology(ring, kmin, kmax);
        else if (direction == FORWARD) out = run_forward(ring, kmin, kmax, cleared_in);
        else out = run_cohomology(ring, kmin, kmax, cleared_in);
        // degrees whose pairs depend on arithmetic (degree 0 is exact union-find when the complex is graph-like;
        // the homology direction also reduces d_kmin to find the essential kmin-classes)
        int lo = direction == HOMOLOGY ? std::max(0, kmin - 1) : kmin;
        for (int k = lo; k <= std::min(kmax, dmax); ++k) {
            if (k == 0 && P.graph_like) continue;
            if (k < P.D) out.note_bound(k, P.hb[k], P.pb[k]);
        }
        return out;
    }

    // homology direction with clearing from the top: reduce d_{j} for j = kmax+1 down to max(kmin, 1);
    // pairs (low, column) of d_{k+1} give the degree-k bars; essential k-cells are those that are neither
    // lows of d_{k+1} nor carry a nonzero reduced column in d_k.
    template <class R>
    BRun<Clear> run_homology(const R& ring, int kmin, int kmax) const {
        typedef typename R::T T;
        struct Ent { int32_t row; T c; };
        BRun<Clear> out; out.kmin = kmin;
        kmax = std::min(kmax, dmax);
        int nd = std::max(0, kmax + 1);
        out.bars.assign(nd, {}); out.clear_out.assign(nd, {});
        for (int l = 0; l < R::lanes; ++l) out.moduli.push_back(ring.modulus(l));
        const T ONE = ring.one(), MONE = ring.neg(ring.one());
        auto conv = [&](int64_t c) -> T { return c == 1 ? ONE : (c == -1 ? MONE : ring.from_int(c)); };
        const double INF = std::numeric_limits<double>::infinity();
        int jtop = std::min(kmax + 1, P.D), jbot = std::max(kmin, 1);
        // is_low[k][s]: k-cell s is the low of some column of d_{k+1}; killer[k][s]: column of d_k nonzero
        std::vector<std::vector<uint8_t>> is_low(P.D + 2), nonzero(P.D + 2);
        for (int k = 0; k <= P.D; ++k) { is_low[k].assign(P.count(k), 0); nonzero[k].assign(P.count(k), 0); }
        bool uf = P.graph_like && jbot <= 1;
        for (int j = jtop; j >= jbot; --j) {
            if (j == 1 && uf) break;                       // d_1 by union-find below
            size_t nj = P.count(j), nl = P.count(j - 1);
            const int64_t* off = P.bd_off[j].data();
            const int32_t* rows = P.bd_row[j].data();
            const int64_t* coefs = P.bd_coef[j].data();
            std::vector<int32_t> piv_col(nl, -1);
            std::vector<T> piv_inv(nl);
            std::vector<std::vector<Ent>> stored(nj);
            std::vector<uint8_t> changed(nj, 0);
            std::vector<Ent> w, tmp;
            int deg = j - 1;
            const std::vector<double>& lbirth = P.birth[j - 1];

            bool record = deg >= kmin && deg <= kmax;
            for (size_t jj = 0; jj < nj; ++jj) {
                int32_t col = (int32_t)jj;
                if (is_low[j][col]) continue;               // clearing: positive cell, column reduces to zero
                ++out.columns;
                w.clear();
                for (int64_t e = off[col]; e < off[col + 1]; ++e) { T x = conv(coefs[e]); if (!ring.is_zero(x)) w.push_back({rows[e], x}); }
                bool ch = false;
                size_t pick = 0;                          // index in w of the chosen pivot entry
                while (!w.empty()) {
                    // top level of w: entries with the birth of w.back(); rows sorted, so a suffix of w
                    double top = lbirth[w.back().row];
                    size_t s0 = w.size(); while (s0 > 0 && lbirth[w[s0 - 1].row] == top) --s0;
                    // an entry of the top level sitting on an existing pivot row must be eliminated
                    size_t hit = SIZE_MAX;
                    for (size_t q = w.size(); q-- > s0;) if (piv_col[w[q].row] >= 0) { hit = q; break; }
                    if (hit == SIZE_MAX) {                // free choice of pivot inside the level
                        pick = w.size() - 1;
                        break;
                    }
                    int32_t r = w[hit].row, i = piv_col[r];
                    T f = ring.neg(ring.mul(w[hit].c, piv_inv[r]));
                    tmp.clear(); size_t a = 0;
                    if (changed[i]) {
                        const std::vector<Ent>& ci = stored[i]; size_t b = 0;
                        while (a < w.size() || b < ci.size()) {
                            if (b == ci.size() || (a < w.size() && w[a].row < ci[b].row)) tmp.push_back(w[a++]);
                            else if (a == w.size() || ci[b].row < w[a].row) { T s = ring.mul(f, ci[b].c); if (!ring.is_zero(s)) tmp.push_back({ci[b].row, s}); ++b; }
                            else { T s = ring.add(w[a].c, ring.mul(f, ci[b].c)); if (!ring.is_zero(s)) tmp.push_back({w[a].row, s}); ++a; ++b; }
                        }
                    } else {
                        int64_t b = off[i], be = off[i + 1];
                        while (a < w.size() || b < be) {
                            if (b == be || (a < w.size() && w[a].row < rows[b])) tmp.push_back(w[a++]);
                            else if (a == w.size() || rows[b] < w[a].row) { T s = ring.mul(f, conv(coefs[b])); if (!ring.is_zero(s)) tmp.push_back({rows[b], s}); ++b; }
                            else { T s = ring.add(w[a].c, ring.mul(f, conv(coefs[b]))); if (!ring.is_zero(s)) tmp.push_back({w[a].row, s}); ++a; ++b; }
                        }
                    }
                    w.swap(tmp); ch = true; ++out.additions;
                }
                if (w.empty()) continue;                      // positive j-cell
                T pc = w[pick].c; int32_t r = w[pick].row;
                if (ring.partly_zero(pc)) out.unlucky = true;
                piv_col[r] = col; piv_inv[r] = ring.inv(pc);
                is_low[j - 1][r] = 1; nonzero[j][col] = 1;
                if (record) {
                    double b = P.birth[deg][r], d = P.birth[j][col];
                    if (d > b || keep_zero_length) out.bars[deg].push_back({b, d});
                    if (!ring.is_pm1(pc)) {
                        Pivot pv; pv.k = deg; pv.key = (uint64_t)col; pv.birth = b; pv.death = d; ring.store(pc, pv); out.nonunit.push_back(std::move(pv));
                    }
                }
                if (ch) { stored[col] = w; changed[col] = 1; }
            }
        }
        if (uf && P.D >= 1) {                                 // d_1 by union-find (elder rule)
            size_t n0 = P.count(0);
            std::vector<int32_t> par(n0); for (size_t i = 0; i < n0; ++i) par[i] = (int32_t)i;
            auto find = [&](int32_t x) { while (par[x] != x) x = par[x] = par[par[x]]; return x; };
            for (size_t el = 0; el < P.count(1); ++el) {
                if (is_low[1][el]) continue;
                auto ev = P.ends[el]; if (ev.first < 0) continue;
                int32_t a = find(ev.first), b = find(ev.second);
                if (a == b) continue;
                int32_t young = std::max(a, b); par[young] = std::min(a, b);
                is_low[0][young] = 1; nonzero[1][el] = 1;
                if (kmin == 0) {
                    double bb = P.birth[0][young], dd = P.birth[1][el];
                    if (dd > bb || keep_zero_length) out.bars[0].push_back({bb, dd});
                }
            }
        }
        for (int k = kmin; k <= kmax; ++k) {
            for (size_t s = 0; s < P.count(k); ++s)
                if (!is_low[k][s] && !nonzero[k][s]) out.bars[k].push_back({P.birth[k][s], INF});
            std::sort(out.bars[k].begin(), out.bars[k].end());
        }
        return out;
    }

    // forward cohomology ("annotations"): cells are inserted in filtration order; every k-cell carries its
    // coordinates in the basis of the currently alive k-classes (indexed by their generating k-cells).
    // A (k+1)-cell whose boundary has zero coordinates creates a class; otherwise it kills the youngest
    // class in its boundary's coordinates, and every k-cell carrying that coordinate is updated (inverted
    // index). Degrees kmin..kmax; the positive kmin-cells come from cleared_in (or union-find when kmin = 0).
    template <class R>
    BRun<Clear> run_forward(const R& ring, int kmin, int kmax, const Clear* cleared_in) const {
        typedef typename R::T T;
        typedef std::vector<std::pair<int32_t, T>> SVec;
        BRun<Clear> out; out.kmin = kmin;
        kmax = std::min(kmax, dmax);
        int nd = std::max(0, kmax + 1);
        out.bars.assign(nd, {}); out.clear_out.assign(nd, {});
        for (int l = 0; l < R::lanes; ++l) out.moduli.push_back(ring.modulus(l));
        const T ONE = ring.one(), MONE = ring.neg(ring.one());
        auto conv = [&](int64_t c) -> T { return c == 1 ? ONE : (c == -1 ? MONE : ring.from_int(c)); };
        const double INF = std::numeric_limits<double>::infinity();
        int ks = kmin;                                   // lowest degree tracked by annotations
        std::vector<uint8_t> negative_low;               // negative ks-cells
        if (kmin == 0 && P.graph_like && P.D >= 1) {     // degree 0 by union-find
            size_t n0 = P.count(0);
            std::vector<int32_t> par(n0); for (size_t i = 0; i < n0; ++i) par[i] = (int32_t)i;
            auto find = [&](int32_t x) { while (par[x] != x) x = par[x] = par[par[x]]; return x; };
            negative_low.assign(P.count(1), 0);
            for (size_t el = 0; el < P.count(1); ++el) {
                auto ev = P.ends[el]; if (ev.first < 0) continue;
                int32_t a = find(ev.first), b = find(ev.second);
                if (a == b) continue;
                int32_t young = std::max(a, b); par[young] = std::min(a, b);
                double bb = P.birth[0][young], dd = P.birth[1][el];
                if (dd > bb || keep_zero_length) out.bars[0].push_back({bb, dd});
                if (dmax > 0) out.clear_out[0].push_back((int32_t)el);
                negative_low[el] = 1;
            }
            for (size_t i = 0; i < n0; ++i) if (find((int32_t)i) == (int32_t)i) out.bars[0].push_back({P.birth[0][i], INF});
            std::sort(out.bars[0].begin(), out.bars[0].end());
            ks = 1;
            if (kmax < 1) return out;
        } else {
            negative_low.assign(P.count(ks), 0);
            if (cleared_in && ks > 0) for (int32_t s : *cleared_in) negative_low[s] = 1;
        }
        int top = std::min(kmax + 1, P.D);
        std::vector<std::vector<SVec>> ann(nd);
        std::vector<std::vector<std::vector<int32_t>>> inv(nd);
        std::vector<std::vector<uint8_t>> alive(nd);
        for (int k = ks; k <= kmax; ++k) { ann[k].resize(P.count(k)); inv[k].resize(P.count(k)); alive[k].assign(P.count(k), 0); }
        // merge cells of dimensions ks..top in filtration order (birth, then dimension)
        std::vector<size_t> pos(top + 2, 0);
        SVec a, tmp;
        auto axpy = [&](const SVec& x, T f, const SVec& y, SVec& o) {   // o = x + f*y
            o.clear(); size_t i = 0, j = 0;
            while (i < x.size() || j < y.size()) {
                if (j == y.size() || (i < x.size() && x[i].first < y[j].first)) { o.push_back(x[i++]); continue; }
                if (i == x.size() || y[j].first < x[i].first) { T v = ring.mul(f, y[j].second); if (!ring.is_zero(v)) o.push_back({y[j].first, v}); ++j; continue; }
                T v = ring.add(x[i].second, ring.mul(f, y[j].second));
                if (!ring.is_zero(v)) o.push_back({x[i].first, v}); ++i; ++j;
            }
        };
        auto lb = [](SVec& v, int32_t key) { return std::lower_bound(v.begin(), v.end(), key, [](const std::pair<int32_t, T>& e, int32_t kk) { return e.first < kk; }); };
        while (true) {
            int d = -1; double bb = INF;
            for (int k = ks; k <= top; ++k) if (pos[k] < P.count(k) && (d < 0 || P.birth[k][pos[k]] < bb)) { d = k; bb = P.birth[k][pos[k]]; }
            if (d < 0) break;
            int32_t t = (int32_t)pos[d]++;
            bool positive;
            if (d == ks) positive = !negative_low[t];
            else {
                // coordinates of the boundary in degree d-1
                a.clear();
                for (int64_t e = P.bd_off[d][t]; e < P.bd_off[d][t + 1]; ++e) {
                    const SVec& s = ann[d - 1][P.bd_row[d][e]];
                    if (s.empty()) continue;
                    T x = conv(P.bd_coef[d][e]); if (ring.is_zero(x)) continue;
                    axpy(a, x, s, tmp); a.swap(tmp);
                }
                positive = a.empty();
                if (!positive) {
                    ++out.columns;
                    int32_t g = a.back().first; T cg = a.back().second;
                    if (ring.partly_zero(cg)) out.unlucky = true;
                    int k = d - 1;
                    double b0 = P.birth[k][g], d0 = P.birth[d][t];
                    if (d0 > b0 || keep_zero_length) out.bars[k].push_back({b0, d0});
                    if (!ring.is_pm1(cg)) {
                        Pivot pv; pv.k = k; pv.key = (uint64_t)t; pv.birth = b0; pv.death = d0; ring.store(cg, pv); out.nonunit.push_back(std::move(pv));
                    }
                    if (k < dmax) out.clear_out[k].push_back(t);
                    alive[k][g] = 0;
                    T fac0 = ring.neg(ring.inv(cg));
                    std::vector<int32_t> zs; zs.swap(inv[k][g]);
                    for (int32_t z : zs) {
                        SVec& rz = ann[k][z];
                        auto it = lb(rz, g);
                        if (it == rz.end() || it->first != g) continue;
                        T fac = ring.mul(fac0, it->second);
                        axpy(rz, fac, a, tmp);
                        for (auto& pr : tmp) { auto q = lb(rz, pr.first); if (q == rz.end() || q->first != pr.first) inv[k][pr.first].push_back(z); }
                        rz.swap(tmp); ++out.additions;
                    }
                }
            }
            if (positive && d <= kmax) {                  // a new class in degree d, generated by t
                ann[d][t] = SVec{{t, ONE}}; inv[d][t].push_back(t); alive[d][t] = 1;
            }
        }
        for (int k = ks; k <= kmax; ++k) {
            for (size_t g = 0; g < P.count(k); ++g) if (alive[k][g]) out.bars[k].push_back({P.birth[k][g], INF});
            std::sort(out.bars[k].begin(), out.bars[k].end());
        }
        return out;
    }

    template <class R>
    BRun<Clear> run_cohomology(const R& ring, int kmin, int kmax, const Clear* cleared_in) const {
        typedef typename R::T T;
        struct Ent { int32_t row; T c; };
        BRun<Clear> out; out.kmin = kmin;
        kmax = std::min(kmax, dmax);
        int nd = std::max(0, kmax + 1);
        out.bars.assign(nd, {}); out.clear_out.assign(nd, {});
        for (int l = 0; l < R::lanes; ++l) out.moduli.push_back(ring.modulus(l));
        const T ONE = ring.one(), MONE = ring.neg(ring.one());
        auto conv = [&](int64_t c) -> T { return c == 1 ? ONE : (c == -1 ? MONE : ring.from_int(c)); };
        const double INF = std::numeric_limits<double>::infinity();

        Clear cleared;
        if (cleared_in) cleared = *cleared_in;
        int kstart = kmin;
        if (kmin == 0 && P.graph_like && kmax >= 0) {   // H0 by union-find (elder rule), ring independent
            size_t n0 = P.count(0);
            std::vector<int32_t> par(n0); for (size_t i = 0; i < n0; ++i) par[i] = (int32_t)i;
            auto find = [&](int32_t x) { while (par[x] != x) x = par[x] = par[par[x]]; return x; };
            for (size_t el = 0; el < P.count(1); ++el) {
                auto ev = P.ends[el]; if (ev.first < 0) continue;
                int32_t a = find(ev.first), b = find(ev.second);
                if (a == b) continue;
                int32_t young = std::max(a, b), old = std::min(a, b);
                par[young] = old;
                double bb = P.birth[0][young], dd = P.birth[1][el];
                if (dd > bb || keep_zero_length) out.bars[0].push_back({bb, dd});
                if (dmax > 0) out.clear_out[0].push_back((int32_t)el);
                ++out.columns;
            }
            for (size_t i = 0; i < n0; ++i) if (find((int32_t)i) == (int32_t)i) out.bars[0].push_back({P.birth[0][i], INF});
            std::sort(out.bars[0].begin(), out.bars[0].end());
            cleared = out.clear_out[0];
            kstart = 1;
        }
        for (int k = kstart; k <= kmax; ++k) {
            size_t nk = P.count(k), nk1 = P.count(k + 1);
            std::vector<uint8_t> skip(nk, 0);
            for (int32_t s : cleared) skip[s] = 1;
            std::vector<int32_t> piv_col(nk1, -1);
            std::vector<T> piv_inv(nk1);
            std::vector<std::vector<Ent>> stored(nk);      // reduced column, if the reduction changed it
            std::vector<uint8_t> changed(nk, 0);
            const int64_t* off = P.cob_off[k].data();
            const int32_t* rows = P.cob_row[k].data();
            const int64_t* coefs = P.cob_coef[k].data();
            const std::vector<double>& bk = P.birth[k];
            std::vector<Ent> w, tmp;
            Clear next_cleared;
            bool top = (k == dmax);                        // pivot rows of the last degree are not needed
            for (size_t jj = nk; jj-- > 0;) {
                int32_t j = (int32_t)jj;
                if (skip[j]) continue;
                ++out.columns;
                w.clear();
                for (int64_t e = off[j]; e < off[j + 1]; ++e) { T x = conv(coefs[e]); if (!ring.is_zero(x)) w.push_back({rows[e], x}); }
                bool ch = false;
                while (!w.empty()) {
                    int32_t r = w[0].row;
                    int32_t i = piv_col[r];
                    if (i < 0) break;
                    T f = ring.neg(ring.mul(w[0].c, piv_inv[r]));
                    tmp.clear();
                    size_t a = 1;          // w[0] cancels by construction
                    if (changed[i]) {
                        const std::vector<Ent>& ci = stored[i];
                        size_t b = 1, be = ci.size();
                        while (a < w.size() && b < be) {
                            if (w[a].row < ci[b].row) tmp.push_back(w[a++]);
                            else if (ci[b].row < w[a].row) { tmp.push_back({ci[b].row, ring.mul(f, ci[b].c)}); ++b; }
                            else { T s = ring.add(w[a].c, ring.mul(f, ci[b].c)); if (!ring.is_zero(s)) tmp.push_back({w[a].row, s}); ++a; ++b; }
                        }
                        for (; a < w.size(); ++a) tmp.push_back(w[a]);
                        for (; b < be; ++b) tmp.push_back({ci[b].row, ring.mul(f, ci[b].c)});
                    } else {
                        int64_t b = off[i], be = off[i + 1];
                        while (b < be && ring.is_zero(conv(coefs[b]))) ++b;
                        ++b;
                        while (a < w.size() && b < be) {
                            if (w[a].row < rows[b]) tmp.push_back(w[a++]);
                            else if (rows[b] < w[a].row) { T s = ring.mul(f, conv(coefs[b])); if (!ring.is_zero(s)) tmp.push_back({rows[b], s}); ++b; }
                            else { T s = ring.add(w[a].c, ring.mul(f, conv(coefs[b]))); if (!ring.is_zero(s)) tmp.push_back({w[a].row, s}); ++a; ++b; }
                        }
                        for (; a < w.size(); ++a) tmp.push_back(w[a]);
                        for (; b < be; ++b) { T s = ring.mul(f, conv(coefs[b])); if (!ring.is_zero(s)) tmp.push_back({rows[b], s}); }
                    }
                    w.swap(tmp); ch = true; ++out.additions;
                }
                if (w.empty()) { out.bars[k].push_back({bk[j], INF}); continue; }
                T pc = w[0].c;
                if (ring.partly_zero(pc)) out.unlucky = true;
                int32_t r = w[0].row;
                piv_col[r] = j; piv_inv[r] = ring.inv(pc);
                double d = P.birth[k + 1][r];
                if (d > bk[j] || keep_zero_length) out.bars[k].push_back({bk[j], d});
                if (!top) next_cleared.push_back(r);
                if (!ring.is_pm1(pc)) {
                    Pivot pv; pv.k = k; pv.key = (uint64_t)j; pv.birth = bk[j]; pv.death = d; ring.store(pc, pv); out.nonunit.push_back(std::move(pv));
                }
                if (ch) { stored[j] = w; changed[j] = 1; }
            }
            std::sort(out.bars[k].begin(), out.bars[k].end());
            out.clear_out[k] = next_cleared;
            cleared.swap(next_cleared);
        }
        return out;
    }
};

} // namespace perseus

namespace perseus {
// direction heuristic (calibrated by bench/direction.py): forward cohomology wins when some degree pair
// (k, k+1) has many more (k+1)-cells than k-cells and the k-cells arrive in few large batches (e.g.
// Linial-Meshulam complexes: all edges at once, triangles one by one). Otherwise backward cohomology
// with clearing (cubical, Rips, general complexes).
inline int auto_direction(const Prepared& P, int maxdim) {
    for (int k = 1; k < std::min(P.D, maxdim + 1); ++k) {
        size_t nk = P.count(k), nk1 = P.count(k + 1);
        if (nk == 0 || nk1 < 4 * nk) continue;
        size_t levels = 1;
        for (size_t s = 1; s < nk; ++s) if (P.birth[k][s] != P.birth[k][s - 1]) ++levels;
        if (levels * 10 <= nk) return 2;
    }
    return 0;
}
}
