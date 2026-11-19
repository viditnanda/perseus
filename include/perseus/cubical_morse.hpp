// cubical_morse.hpp -- persistent (co)homology of a top-cell cubical filtration through a filtered discrete
// gradient (Mischaikow-Nanda 2013), with the critical cells reduced instead of all cells.
//
// The filtration. Top cubes are totally ordered by (birth index, position); a cell enters with its first top
// coface tau. The cells entering with tau form the block L(tau) = cl(tau) - A(tau), where A(tau) is the union
// of the faces tau shares with the earlier top cubes among its 3^d - 1 neighbours. A(tau) depends only on
// which neighbours are earlier, so every block is one of finitely many "patterns" (subcomplexes of the
// boundary of the d-cube).
//
// The gradient. Each pattern gets an acyclic matching on L(tau): a deterministic greedy pass, then, when it
// leaves more critical cells than rank H_*(cl tau, A; F_2) (the minimum), an exact search over removal
// orders. Every pattern met in 2D and 3D is perfect, and so is every pattern we have met in 4D; a pattern
// that is not costs only extra critical cells, never correctness. Matchings live inside blocks, so the
// gradient is filtered, and the Morse complex is filtered chain-homotopy equivalent (over Z) to the cubical
// complex. Barcodes therefore agree over every coefficient ring, and so does the bad-prime list. Positive-
// length bars depend only on the filtration by values; zero-length bars (--zero-length) are recovered by
// counting, n_k(v) = births_k(v) + deaths_{k-1}(v).
//
// Order inside a block: a topological order of the block's modified Hasse diagram (matched pairs point up,
// other incidences down). With ord = reverse topological position, the key (birth, tau position, ord)
// strictly increases along every coboundary V-path step and strictly decreases along every boundary step,
// so it is a filtration order on the Morse complex and a processing order for path sums.
//
// The reductions. H0 is union-find on critical vertices: the Morse boundary of a critical edge joins the
// ends of the two descending V-paths from its vertices (memoised, so every vertex is walked once). Degree
// d-1 is the dual union-find on critical top cells, the same way. Degrees 1..d-2 use the explicit Morse
// coboundary (integer path sums, computed once) and cohomology with clearing, as in the classic engine.
#pragma once
#include "io.hpp"
#include "ring.hpp"
#include "bitcol.hpp"
#include "types.hpp"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <cstring>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <functional>
#include <map>

namespace perseus {

namespace cubdmt {

struct M128 {
    uint64_t lo = 0, hi = 0;
    bool test(int i) const { return i < 64 ? (lo >> i) & 1 : (hi >> (i - 64)) & 1; }
    void set(int i) { if (i < 64) lo |= 1ULL << i; else hi |= 1ULL << (i - 64); }
    void clr(int i) { if (i < 64) lo &= ~(1ULL << i); else hi &= ~(1ULL << (i - 64)); }
    M128& operator|=(const M128& o) { lo |= o.lo; hi |= o.hi; return *this; }
    bool operator==(const M128& o) const { return lo == o.lo && hi == o.hi; }
    bool any() const { return lo || hi; }
};

// one block's matching: for each cell of L, its local index, its mate code (0 critical, else 1 + 2*axis +
// side: the mate is the cell at local offset +-1 along axis) and its order inside the block
struct PCell { uint8_t loc, code, ord; };
struct Pattern { uint32_t start = 0; uint16_t len = 0; uint8_t critical = 0, minimum = 0; };

// the closed d-cube as a local complex: cells = offset vectors in {-1,0,1}^d, index sum (o_i + 1) 3^i
class LocalCube {
public:
    int d = 0, n = 0, center = 0;
    std::vector<int> dim;
    std::vector<std::vector<int>> faces, cofaces;
    std::vector<M128> up;               // up[e]: the cells of the face shared with the neighbour at offset e
    std::vector<M128> facm, cofm;       // faces / cofaces of each cell as masks
    std::vector<int> off;               // off[i * d + a] in {-1, 0, 1}
    std::vector<int> posd;              // index of each cell among the cells of its dimension
    M128 dmask[17];                     // the cells of each dimension
    mutable std::vector<int> m1, m2, ordv, sq, sq2;   // scratch
    long long exact_calls = 0, imperfect = 0;

    explicit LocalCube(int d_) : d(d_) {
        n = 1; for (int i = 0; i < d; ++i) n *= 3;
        if (n > 128) throw std::runtime_error("LocalCube: dimension too large");
        off.assign((size_t)n * d, 0); dim.assign(n, 0); faces.assign(n, {}); cofaces.assign(n, {});
        for (int i = 0; i < n; ++i) { int r = i; for (int a = 0; a < d; ++a) { off[i * d + a] = r % 3 - 1; r /= 3; if (off[i * d + a] == 0) ++dim[i]; } }
        center = (n - 1) / 2;
        { int cnt[17] = {0}; posd.assign(n, 0); for (int i = 0; i < n; ++i) { posd[i] = cnt[dim[i]]++; dmask[dim[i]].set(i); } }
        int p3[16]; p3[0] = 1; for (int a = 1; a < d; ++a) p3[a] = 3 * p3[a - 1];
        for (int i = 0; i < n; ++i)
            for (int a = 0; a < d; ++a) if (off[i * d + a] == 0)
                for (int s : {-1, 1}) { int j = i + s * p3[a]; faces[i].push_back(j); cofaces[j].push_back(i); }
        facm.assign(n, M128()); cofm.assign(n, M128());
        for (int i = 0; i < n; ++i) { for (int f : faces[i]) facm[i].set(f); for (int c : cofaces[i]) cofm[i].set(c); }
        up.assign(n, M128());
        for (int e = 0; e < n; ++e) {
            if (e == center) continue;
            for (int o = 0; o < n; ++o) {
                bool in = true;
                for (int a = 0; a < d && in; ++a) if (off[e * d + a] != 0 && off[e * d + a] != off[o * d + a]) in = false;
                if (in) up[e].set(o);
            }
        }
    }

    // sum of F_2 ranks of H_*(cl tau, A)
    int minimum(const M128& A) const {
        int nl = 0; for (int i = 0; i < n; ++i) if (!A.test(i)) ++nl;
        const std::vector<int>& pos = posd;
        int rk = 0;
        for (int p = 1; p <= d; ++p) {
            uint64_t rows[64]; bool have[64] = {false};
            for (int c = 0; c < n; ++c) {
                if (dim[c] != p || A.test(c)) continue;
                uint64_t v = 0;
                for (int f : faces[c]) if (!A.test(f)) v |= 1ULL << pos[f];
                while (v) {
                    int h = 63 - __builtin_clzll(v);
                    if (have[h]) v ^= rows[h]; else { have[h] = true; rows[h] = v; ++rk; break; }
                }
            }
        }
        return nl - 2 * rk;
    }

    // greedy removal from the top (and, if coreduce, from the bottom), on bit masks; returns the mate array
    // (-1 critical, -2 not in L) and the number of critical cells
    int greedy(const M128& A, bool coreduce, std::vector<int>& mate, std::vector<int>* seq = nullptr) const {
        M128 live; for (int i = 0; i < n; ++i) if (!A.test(i)) live.set(i);
        mate.assign(n, -2);
        if (seq) seq->clear();
        for (int i = 0; i < n; ++i) if (live.test(i)) mate[i] = -1;
        auto cnt = [&](const M128& m) { return __builtin_popcountll(m.lo & live.lo) + __builtin_popcountll(m.hi & live.hi); };
        auto none = [&](const M128& m) { return !(m.lo & live.lo) && !(m.hi & live.hi); };
        int crit = 0;
        while (live.any()) {
            int pf = -1, pc = -1;
            for (int w = 1; w >= 0 && pf < 0; --w) {                  // top: c maximal, f's only live coface is c
                uint64_t bits = w ? live.hi : live.lo;
                while (bits && pf < 0) {
                    int c = 63 - __builtin_clzll(bits) + 64 * w; bits &= ~(1ULL << (c - 64 * w));
                    if (!none(cofm[c])) continue;
                    for (int f : faces[c]) if (live.test(f) && cnt(cofm[f]) == 1) { pf = f; pc = c; break; }
                }
            }
            if (pf < 0 && coreduce)
                for (int f = 0; f < n && pf < 0; ++f) {                // bottom: f minimal, c's only live face is f
                    if (!live.test(f) || !none(facm[f])) continue;
                    for (int c : cofaces[f]) if (live.test(c) && cnt(facm[c]) == 1) { pf = f; pc = c; break; }
                }
            if (pf >= 0) { mate[pf] = pc; mate[pc] = pf; live.clr(pf); live.clr(pc); if (seq) { seq->push_back(pf); seq->push_back(pc); } continue; }
            int best = -1;                                            // stuck: a maximal cell of highest dimension
            for (int c = 0; c < n; ++c) if (live.test(c) && none(cofm[c]) && (best < 0 || dim[c] > dim[best])) best = c;
            ++crit; live.clr(best); if (seq) seq->push_back(best);
        }
        return crit;
    }

    // exact: can L be removed from the top with at most `budget` critical cells? fills mate on success
    bool exact(const M128& A, int budget, std::vector<int>& mate, long long cap) const {
        M128 full; for (int i = 0; i < n; ++i) if (!A.test(i)) full.set(i);
        struct H { size_t operator()(const std::pair<M128, int>& k) const { return std::hash<uint64_t>()(k.first.lo * 0x9E3779B97F4A7C15ULL ^ k.first.hi ^ ((uint64_t)k.second << 56)); } };
        struct E { bool operator()(const std::pair<M128, int>& a, const std::pair<M128, int>& b) const { return a.first == b.first && a.second == b.second; } };
        std::unordered_set<std::pair<M128, int>, H, E> dead;
        long long states = 0;
        std::vector<std::pair<int, int>> moves;            // (f, c) pair or (-1, c) critical
        std::function<bool(M128, int)> go = [&](M128 live, int b) -> bool {
            if (!live.any()) return true;
            if (dead.count({live, b})) return false;
            if (++states > cap) return false;
            std::vector<int> mx;
            for (int c = 0; c < n; ++c) {
                if (!live.test(c)) continue;
                bool m = true; for (int j : cofaces[c]) if (live.test(j)) { m = false; break; }
                if (m) mx.push_back(c);
            }
            for (int c : mx)
                for (int f : faces[c]) {
                    if (!live.test(f)) continue;
                    int k = 0; for (int j : cofaces[f]) k += live.test(j);
                    if (k != 1) continue;
                    M128 nl = live; nl.clr(c); nl.clr(f);
                    moves.push_back({f, c});
                    if (go(nl, b)) return true;
                    moves.pop_back();
                }
            if (b > 0)
                for (int c : mx) {
                    M128 nl = live; nl.clr(c);
                    moves.push_back({-1, c});
                    if (go(nl, b - 1)) return true;
                    moves.pop_back();
                }
            if (states <= cap) dead.insert({live, b});
            return false;
        };
        if (!go(full, budget)) return false;
        mate.assign(n, -2);
        for (int i = 0; i < n; ++i) if (!A.test(i)) mate[i] = -1;
        for (auto& m : moves) if (m.first >= 0) { mate[m.first] = m.second; mate[m.second] = m.first; }
        return true;
    }

    // reverse topological positions of the modified Hasse diagram on L (matched f -> c, otherwise c -> f);
    // false if it has a cycle
    bool order(const M128& A, const std::vector<int>& mate, std::vector<int>& ord) const {
        // successors of u: its matched coface, and its faces in L other than a face it is matched with
        int indeg[128] = {0}, q[128]; int nl = 0, qn = 0;
        for (int c = 0; c < n; ++c) {
            if (A.test(c)) continue;
            ++nl;
            for (int f : faces[c]) { if (A.test(f)) continue; if (mate[f] == c) ++indeg[c]; else ++indeg[f]; }
        }
        for (int i = 0; i < n; ++i) if (!A.test(i) && indeg[i] == 0) q[qn++] = i;
        ord.assign(n, -1); int t = 0;
        for (int h = 0; h < qn; ++h) {
            int u = q[h]; ord[u] = t++;
            if (mate[u] >= 0) for (int c : cofaces[u]) if (c == mate[u] && --indeg[c] == 0) q[qn++] = c;
            for (int f : faces[u]) if (!A.test(f) && mate[f] != u && --indeg[f] == 0) q[qn++] = f;
        }
        if (t != nl) return false;
        for (int i = 0; i < n; ++i) if (ord[i] >= 0) ord[i] = nl - 1 - ord[i];
        return true;
    }

    // the block matching for the pattern A, appended to pool. A removal from the top lists every matched
    // face just before its coface and every other incidence coface first, so its sequence is a topological
    // order of the modified Hasse diagram. The greedy count is optimal when it equals |chi(L)|, a lower
    // bound for sum rank H_*(cl tau, A); otherwise the rank decides and, if needed, an exact search.
    Pattern solve(const M128& A, std::vector<PCell>& pool) {
        Pattern P;
        std::vector<int>& mate = m1; std::vector<int>& ord = ordv;
        std::vector<int>& seq = sq;
        int crit = greedy(A, false, mate, &seq);
        int nl = (int)seq.size();
        ord.assign(n, -1);
        for (int t = 0; t < nl; ++t) ord[seq[t]] = nl - 1 - t;
        int chi = 0;
        for (int p = 0; p <= d; ++p) {
            int c = __builtin_popcountll(dmask[p].lo & ~A.lo) + __builtin_popcountll(dmask[p].hi & ~A.hi);
            chi += (p & 1) ? -c : c;
        }
        int mn = crit == std::abs(chi) ? crit : minimum(A);
        P.minimum = (uint8_t)mn;
        if (crit > mn) {
            int c2 = greedy(A, true, m2);
            if (c2 < crit && order(A, m2, sq2)) { crit = c2; mate.swap(m2); ord.swap(sq2); }
        }
        if (crit > mn) {
            ++exact_calls;
            for (int bgt = mn; bgt < crit; ++bgt) if (exact(A, bgt, m2, 200000)) { crit = bgt; mate.swap(m2); if (!order(A, mate, ord)) throw std::runtime_error("cubical morse: cyclic block matching"); break; }
        }
        if (crit > mn) ++imperfect;
        P.critical = (uint8_t)crit; P.start = (uint32_t)pool.size();
        for (int i = 0; i < n; ++i) {
            if (A.test(i)) continue;
            PCell c; c.loc = (uint8_t)i; c.ord = (uint8_t)ord[i]; c.code = 0;
            if (mate[i] >= 0)
                for (int a = 0; a < d; ++a) { int df = off[mate[i] * d + a] - off[i * d + a]; if (df) c.code = (uint8_t)(1 + 2 * a + (df > 0)); }
            pool.push_back(c);
        }
        P.len = (uint16_t)(pool.size() - P.start);
        return P;
    }
};

// open-addressing map from a 128-bit block mask to a pattern index
class PatternCache {
public:
    std::vector<M128> keys; std::vector<int32_t> vals; size_t mask = 0, used = 0;
    PatternCache() { rehash(1 << 12); }
    static uint64_t h(const M128& k) { uint64_t x = k.lo * 0x9E3779B97F4A7C15ULL ^ (k.hi + 0x632BE59BD9B4E019ULL) * 0xC2B2AE3D27D4EB4FULL; return x ^ (x >> 29); }
    int32_t* find_or_slot(const M128& k) {
        size_t i = h(k) & mask;
        while (vals[i] >= 0) { if (keys[i] == k) return &vals[i]; i = (i + 1) & mask; }
        return &vals[i];
    }
    int32_t get(const M128& k) const {
        size_t i = h(k) & mask;
        while (vals[i] >= 0) { if (keys[i] == k) return vals[i]; i = (i + 1) & mask; }
        return -1;
    }
    void put(const M128& k, int32_t v) {
        if (2 * (used + 1) > keys.size()) rehash(keys.size() * 2);
        size_t i = h(k) & mask;
        while (vals[i] >= 0) i = (i + 1) & mask;
        keys[i] = k; vals[i] = v; ++used;
    }
    void rehash(size_t sz) {
        std::vector<M128> ok; std::vector<int32_t> ov; ok.swap(keys); ov.swap(vals);
        keys.assign(sz, M128()); vals.assign(sz, -1); mask = sz - 1; used = 0;
        for (size_t i = 0; i < ok.size(); ++i) if (ov[i] >= 0) put(ok[i], ov[i]);
    }
};

// position -> int64 scratch map, cleared in O(1) by a stamp; keeps its insertion list
class PosMap {
public:
    PosMap() { grow(1 << 10); }
    void reset() { if (++cur == 0) { std::fill(stamp.begin(), stamp.end(), 0); cur = 1; } list.clear(); }
    // the slot of k, inserted with value 0 if absent (inserted = true)
    int64_t& at(uint32_t k, bool& inserted) {
        if (2 * (list.size() + 1) > key.size()) grow(key.size() * 2);
        size_t i = h(k) & mask;
        while (stamp[i] == cur) { if (key[i] == k) { inserted = false; return val[i]; } i = (i + 1) & mask; }
        stamp[i] = cur; key[i] = k; val[i] = 0; list.push_back((uint32_t)i); inserted = true; return val[i];
    }
    size_t size() const { return list.size(); }
    uint32_t key_of(size_t j) const { return key[list[j]]; }
    int64_t val_of(size_t j) const { return val[list[j]]; }
private:
    std::vector<uint32_t> key, stamp, list; std::vector<int64_t> val; uint32_t cur = 1; size_t mask = 0;
    static size_t h(uint32_t k) { return (size_t)((uint64_t)k * 0x9E3779B97F4A7C15ULL >> 20); }
    void grow(size_t sz) {
        std::vector<uint32_t> ok = key, os = stamp, ol = list; std::vector<int64_t> ov = val; uint32_t oc = cur;
        key.assign(sz, 0); stamp.assign(sz, 0); val.assign(sz, 0); list.clear(); mask = sz - 1; cur = 1;
        for (uint32_t j : ol) if (os[j] == oc) { bool ins; at(ok[j], ins) = ov[j]; }
    }
};

} // namespace cubdmt

class CubicalMorse {
public:
    typedef std::vector<uint32_t> Clear;     // critical indices (in the next degree) of the pivots
    bool keep_zero_length = false;
    static constexpr uint32_t ABSENT = 0xFFFFFFFFu;

    // whether this engine handles the grid (dimension 2..4, no periodic axes)
    static bool supports(const CubicalGrid& g) {
        int d = (int)g.ext.size();
        if (d < 2 || d > 4) return false;
        for (size_t i = 0; i < g.periodic.size(); ++i) if (g.periodic[i]) return false;
        return true;
    }

    CubicalMorse(const CubicalGrid& g, int maxdim) {
        d = (int)g.ext.size();
        if (!supports(g)) throw std::runtime_error("cubical morse: unsupported grid");
        dmax = std::min(maxdim < 0 ? d : maxdim, d);
        uint64_t total = 1;
        for (int i = 0; i < d; ++i) {
            L[i] = 2 * g.ext[i] + 1; stride[i] = (int64_t)total; total *= (uint64_t)L[i];
            recip[i] = (uint64_t)(((unsigned __int128)1 << 64) / (uint64_t)L[i]) + 1;
            nvx[i] = g.ext[i] + 1; ntx[i] = g.ext[i];
        }
        if (total >= (1ULL << 32)) throw std::runtime_error("cubical: grid too large for the implicit engine");
        N = (int64_t)total;
        auto T0 = std::chrono::steady_clock::now(); dbg = getenv("PERSEUS_DEBUG") != nullptr;
        auto lapt = [&](const char* w) { if (!dbg) return; auto t = std::chrono::steady_clock::now(); fprintf(stderr, "    [morse %s %.3f s]\n", w, std::chrono::duration<double>(t - T0).count()); T0 = t; };
        births(g);
        lapt("births");
        gradient();
        lapt("gradient");
        collect();
        lapt("critical cells");
        h0_ends();
        lapt("H0 paths");
        if (dmax >= d - 1) top_ends();
        lapt("dual paths");
        morse_columns();
        lapt("Morse complex");
    }

    int maxdim() const { return dmax; }
    size_t critical(int k) const { return k >= 0 && k <= d ? crit[k].size() : 0; }
    size_t count(int k) const { return k >= 0 && k <= d ? cnt_[k] : 0; }
    size_t patterns() const { return pats.size(); }
    long long imperfect_blocks() const { return imperfect_blocks_; }
    // --- test hooks -------------------------------------------------------------------------------------
    // nonzero entries of delta_M^(k+1) delta_M^(k) over Z, over all pairs of stored degrees
    long long check_dd() const {
        long long bad = 0;
        for (int k = 1; k + 1 < (int)cols.size(); ++k) {
            const MCols &A = cols[k], &B2 = cols[k + 1];
            if (A.start.empty() || B2.start.empty()) continue;
            for (size_t j = 0; j + 1 < A.start.size(); ++j) {
                std::unordered_map<int32_t, int64_t> acc;
                for (uint64_t q = A.start[j]; q < A.start[j + 1]; ++q)
                    for (uint64_t r = B2.start[A.row[q]]; r < B2.start[A.row[q] + 1]; ++r) acc[B2.row[r]] += A.coef[q] * B2.coef[r];
                for (auto& kv : acc) if (kv.second) ++bad;
            }
        }
        return bad;
    }
    // coboundary entries pointing to critical cells that are not later in the filtration order
    long long check_order() const {
        long long bad = 0;
        for (int k = 1; k < (int)cols.size(); ++k) {
            const MCols& M = cols[k]; if (M.start.empty()) continue;
            for (size_t j = 0; j + 1 < M.start.size(); ++j)
                for (uint64_t q = M.start[j]; q < M.start[j + 1]; ++q) {
                    int32_t r = M.row[q];
                    bool later = ckey[k + 1][r] > ckey[k][j] || (ckey[k + 1][r] == ckey[k][j] && cord[k + 1][r] > cord[k][j]);
                    if (!later) ++bad;
                }
        }
        return bad;
    }
    // stored columns that differ from a recomputation by memoised recursion over V-paths; apparent-pair
    // flags that are wrong; H0 / dual ends that differ from a direct walk
    long long check_columns() const {
        long long bad = 0;
        for (int k = 1; k < (int)cols.size(); ++k) {
            const MCols& M = cols[k];
            if (M.start.empty()) continue;
            std::unordered_map<int64_t, std::map<int64_t, int64_t>> memo;
            std::function<const std::map<int64_t, int64_t>&(int64_t)> R = [&](int64_t u) -> const std::map<int64_t, int64_t>& {
                auto it = memo.find(u); if (it != memo.end()) return it->second;
                std::map<int64_t, int64_t> out;
                int64_t cu[16]; decode(u, cu);
                int64_t skip = mate_pos(u);
                int odd = 0;
                for (int i = 0; i < d; ++i) {
                    if (cu[i] & 1) { ++odd; continue; }
                    for (int s2 = 0; s2 < 2; ++s2) {
                        int64_t x = cu[i] + (s2 ? 1 : -1);
                        if (x < 0 || x >= L[i]) continue;
                        int64_t rho = u + (s2 ? 1 : -1) * stride[i];
                        if (B(rho) == ABSENT || rho == skip) continue;
                        int64_t inc = (s2 == 0) != (odd & 1) ? 1 : -1;
                        if (MT(rho) == 0) { out[rho] += inc; continue; }
                        int ax = (MT(rho) - 1) >> 1; int64_t sd = (MT(rho) - 1) & 1 ? 1 : -1;
                        int64_t crho[16]; std::memcpy(crho, cu, sizeof(int64_t) * d); crho[i] = x;
                        if (!(crho[ax] & 1)) continue;
                        int oddb = 0; for (int a = 0; a < ax; ++a) oddb += (int)(crho[a] & 1);
                        int64_t incf = ((sd > 0) != (oddb & 1)) ? 1 : -1;
                        for (auto& kv : R(rho + sd * stride[ax])) out[kv.first] += -inc * incf * kv.second;
                    }
                }
                return memo[u] = out;
            };
            std::vector<int32_t> yf(crit[k + 1].size(), -1);
            for (size_t j = 0; j + 1 < M.start.size(); ++j) {
                std::map<int64_t, int64_t> want;
                for (auto& kv : R(crit[k][j])) if (kv.second) want[crit_index(k + 1, (uint32_t)kv.first)] = kv.second;
                std::map<int64_t, int64_t> got;
                for (uint64_t q = M.start[j]; q < M.start[j + 1]; ++q) { got[M.row[q]] = M.coef[q]; yf[M.row[q]] = (int32_t)j; }
                if (want != got) ++bad;
            }
            for (size_t j = 0; j + 1 < M.start.size(); ++j) {
                bool app = M.start[j + 1] > M.start[j] && (M.coef[M.start[j]] == 1 || M.coef[M.start[j]] == -1) && yf[M.row[M.start[j]]] == (int32_t)j;
                if (app != (bool)M.apparent[j]) ++bad;
            }
        }
        return bad;
    }

    template <class R>
    BRun<Clear> run(const R& ring, int kmin, int kmax, const Clear* cleared_in) const { return run(ring, kmin, kmax, cleared_in, keep_zero_length); }

    template <class R>
    BRun<Clear> run(const R& ring, int kmin, int kmax, const Clear* cleared_in, bool keep) const {
        typedef typename R::T T;
        kmax = std::min(kmax, dmax);
        if (keep && kmin > 0) {
            // zero-length bars of degree k are counted from degree k-1: run from degree 0, report from kmin
            BRun<Clear> full = run(ring, 0, kmax, nullptr, true);
            BRun<Clear> out; out.kmin = kmin; out.moduli = full.moduli; out.unlucky = full.unlucky;
            out.bars.assign(full.bars.size(), {}); out.clear_out = full.clear_out;
            for (size_t k = kmin; k < full.bars.size(); ++k) out.bars[k] = full.bars[k];
            for (auto& pv : full.nonunit) if (pv.k >= kmin) out.nonunit.push_back(pv);
            out.additions = full.additions; out.columns = full.columns; out.hbits = full.hbits; out.pairs = full.pairs;
            return out;
        }
        struct Ent { int32_t row; T c; };
        BRun<Clear> out; out.kmin = kmin;
        int nd = std::max(0, kmax + 1);
        out.bars.assign(nd, {}); out.clear_out.assign(nd, {});
        for (int l = 0; l < R::lanes; ++l) out.moduli.push_back(ring.modulus(l));
        const T ONE = ring.one(), MONE = ring.neg(ring.one());
        const double INF = std::numeric_limits<double>::infinity();
        std::vector<char> clr;                             // cleared critical cells of the current degree
        auto set_cleared = [&](int k, const Clear& cl) {         // clearing data: critical indices in degree k
            clr.assign(crit[k].size(), 0);
            for (uint32_t i : cl) clr[i] = 1;
        };
        Clear cleared; if (cleared_in) cleared = *cleared_in;
        for (int k = kmin; k <= kmax; ++k) {
            set_cleared(k, cleared);
            Clear next;
            bool top = (k == dmax);
            const std::vector<uint32_t>& C = crit[k];
            const std::vector<double>& CV = cval[k];
            if (k == d) {
                for (size_t i = 0; i < C.size(); ++i) { if (clr[i]) continue; ++out.columns; out.bars[k].push_back({CV[i], INF}); }
            } else if (k == 0) {
                // elder rule on critical vertices (index = filtration order: the representative is the oldest)
                size_t nv = C.size();
                std::vector<int32_t> par(nv); for (size_t i = 0; i < nv; ++i) par[i] = (int32_t)i;
                auto find = [&](int32_t x) { while (par[x] != x) x = par[x] = par[par[x]]; return x; };
                const std::vector<uint32_t>& E = crit[1]; const std::vector<double>& EV = cval[1];
                for (size_t e = 0; e < E.size(); ++e) {
                    int32_t a = find(ends0[e].first), b = find(ends0[e].second);
                    if (a == b) continue;
                    int32_t young = std::max(a, b), old = std::min(a, b);
                    par[young] = old;
                    double bv = CV[young], dv = EV[e];
                    if (dv > bv) out.bars[0].push_back({bv, dv});
                    if (!top) next.push_back((uint32_t)e);
                    ++out.columns;
                }
                for (size_t i = 0; i < nv; ++i) if (find((int32_t)i) == (int32_t)i) out.bars[0].push_back({CV[i], INF});
            } else if (k == d - 1) {
                // dual union-find on critical top cells plus the virtual node V (index ntop, youngest of all)
                const std::vector<uint32_t>& Tc = crit[d]; const std::vector<double>& TV = cval[d];
                int32_t V = (int32_t)Tc.size();
                std::vector<int32_t> par(Tc.size() + 1); for (size_t i = 0; i < par.size(); ++i) par[i] = (int32_t)i;
                auto find = [&](int32_t x) { while (par[x] != x) x = par[x] = par[par[x]]; return x; };
                for (size_t jj = C.size(); jj-- > 0;) {
                    if (clr[jj]) continue;
                    ++out.columns;
                    double bj = CV[jj];
                    int32_t ea = endsd[jj].first, eb = endsd[jj].second;
                    if (ea < 0) ea = V;
                    if (eb < 0) eb = V;
                    int32_t ra = find(ea), rb = find(eb);
                    if (ra == rb) { out.bars[k].push_back({bj, INF}); continue; }
                    int32_t lo2 = std::min(ra, rb), hi2 = std::max(ra, rb);   // representatives are the youngest
                    par[lo2] = hi2;
                    double dv = TV[lo2];
                    if (dv > bj) out.bars[k].push_back({bj, dv});
                    if (!top) next.push_back((uint32_t)lo2);
                }
            } else {
                // cohomology with clearing on the explicit Morse coboundary; apparent pairs need no column
                const MCols& M = cols[k];
                const std::vector<double>& CnV = cval[k + 1];
                size_t nk = C.size(), nk1 = CnV.size();
                out.note_bound(k, M.hbits, (double)std::min(nk, nk1));
                std::vector<int32_t> piv_col(nk1, -1);
                std::vector<T> piv_inv(nk1);
                std::unordered_map<int32_t, std::vector<Ent>> stored;
                auto conv = [&](int64_t c) { return c == 1 ? ONE : c == -1 ? MONE : ring.from_int(c); };
                auto column = [&](int32_t j, std::vector<Ent>& o) {
                    o.clear();
                    for (uint64_t q = M.start[j]; q < M.start[j + 1]; ++q) {       // path sums vanish mod p when p divides them
                        T c = conv(M.coef[q]);
                        if (!ring.is_zero(c)) o.push_back({M.row[q], c});
                    }
                };
                // Working columns are sorted vectors merged with each added column; a column longer than
                // LONG moves into a bit tree (bitcol.hpp), where duplicate rows cancel on arrival and the
                // pivot is found by descending the tree. The tree is allocated on first use in a degree.
                std::vector<Ent> w, tmp, cob;
                BitCol<R> tc; bool tc_ready = false;
                const size_t LONG = long_column_threshold();
                size_t dbg_maxw = 0, dbg_stored = 0; long long dbg_app = 0, add0 = out.additions, dbg_merged = 0; auto dbg_t0 = std::chrono::steady_clock::now();
                for (size_t jj = nk; jj-- > 0;) {
                    int32_t j = (int32_t)jj;
                    if (clr[j]) continue;
                    ++out.columns;
                    double bj = CV[j];
                    if (M.apparent[j]) {
                        ++dbg_app;
                        int32_t r = M.row[M.start[j]];
                        piv_col[r] = j; piv_inv[r] = conv(M.coef[M.start[j]]);     // +-1 is its own inverse
                        double dv = CnV[r];
                        if (dv > bj) out.bars[k].push_back({bj, dv});
                        if (!top) next.push_back((uint32_t)r);
                        continue;
                    }
                    column(j, w);
                    bool ch = false, found = false, treemode = false; Ent piv{0, ring.zero()};
                    while (true) {
                        found = false;
                        if (!treemode) { if (!w.empty()) { piv = w[0]; found = true; } }
                        else if (!tc.empty()) { uint32_t r0 = tc.min_row(); piv = {(int32_t)r0, tc.coef(ring, r0)}; found = true; }
                        if (!found) break;
                        int32_t i = piv_col[piv.row];
                        if (i < 0) break;
                        T f = ring.neg(ring.mul(piv.c, piv_inv[piv.row]));
                        const std::vector<Ent>* ci;
                        auto it = stored.find(i);
                        if (it != stored.end()) ci = &it->second; else { column(i, cob); ci = &cob; }
                        if (!treemode && w.size() > LONG) {
                            treemode = true;
                            if (!tc_ready) { tc.init(nk1); tc_ready = true; }
                            for (const Ent& e : w) tc.add(ring, (uint32_t)e.row, e.c);
                        }
                        if (treemode) {
                            tc.erase((uint32_t)piv.row);           // the pivot cancels exactly
                            for (size_t q = 1; q < ci->size(); ++q) {
                                T sc = ring.mul(f, (*ci)[q].c);
                                if (!ring.is_zero(sc)) tc.add(ring, (uint32_t)(*ci)[q].row, sc);
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
                        ch = true; ++out.additions; dbg_merged += (long long)ci->size(); dbg_maxw = std::max(dbg_maxw, treemode ? tc.count : w.size());
                    }
                    if (!found) { out.bars[k].push_back({bj, INF}); continue; }
                    T pc = piv.c; int32_t r = piv.row;
                    if (ring.partly_zero(pc)) out.unlucky = true;
                    piv_col[r] = j; piv_inv[r] = ring.inv(pc);
                    double dv = CnV[r];
                    if (dv > bj) out.bars[k].push_back({bj, dv});
                    if (!top) next.push_back((uint32_t)r);
                    if (!ring.is_pm1(pc)) { Pivot pv; pv.k = k; pv.key = (uint64_t)j; pv.birth = bj; pv.death = dv; ring.store(pc, pv); out.nonunit.push_back(std::move(pv)); }
                    if (ch && !treemode) { stored[j] = w; dbg_stored += w.size(); }
                    else if (ch) {                                // the tree still holds the pivot, its smallest row
                        std::vector<Ent>& RC = stored[j]; RC.clear(); RC.reserve(tc.count);
                        tc.drain(ring, [&](uint32_t r0, T c0) { RC.push_back({(int32_t)r0, c0}); });
                        dbg_stored += RC.size();
                    }
                }
                if (dbg) fprintf(stderr, "    [morse reduction degree %d: %.3f s, %lld columns (%lld apparent), %lld additions, %lld merged entries, longest working column %zu, %zu stored entries]\n",
                                 k, std::chrono::duration<double>(std::chrono::steady_clock::now() - dbg_t0).count(), out.columns, dbg_app, out.additions - add0, dbg_merged, dbg_maxw, dbg_stored);
            }
            std::sort(out.bars[k].begin(), out.bars[k].end());
            out.clear_out[k] = next;
            cleared.swap(next);
        }
        if (keep) add_zero_length(out, kmax);
        return out;
    }

private:
    int d = 0, dmax = 0; bool dbg = false;
    int64_t N = 0, L[16], stride[16], nvx[16], ntx[16];
    uint64_t recip[16];
    std::vector<double> vals;           // sorted distinct input values
    // per position, packed: birth index (bits 0-31, ABSENT if not in the complex), mate code (32-39: 0 critical,
    // else 1 + 2*axis + side), local index in its block (40-47), order inside the block (48-55). Once the
    // critical cells are sorted, each one holds instead its index in its dimension (bits 32-62) and bit 63
    // (when the counts allow; idx_in_cw), and its key moves to ckey / cord.
    std::vector<uint64_t> cw;
    std::vector<int64_t> offpos;        // position offset of each local cell of a block
    std::vector<cubdmt::Pattern> pats;
    std::vector<cubdmt::PCell> pool;
    long long imperfect_blocks_ = 0;
    std::vector<size_t> cnt_;
    std::vector<std::vector<uint32_t>> crit;                        // critical positions by dimension, filtration order
    std::vector<std::vector<double>> cval;                          // their values
    std::vector<std::vector<uint64_t>> ckey; std::vector<std::vector<uint8_t>> cord;   // and their keys
    bool idx_in_cw = false;
    std::vector<std::vector<std::pair<uint32_t, int32_t>>> bypos;   // (position, index), sorted by position
    std::vector<std::pair<int32_t, int32_t>> ends0, endsd;          // terminal critical cells of the two V-paths
    struct MCols { std::vector<uint64_t> start; std::vector<int32_t> row; std::vector<int64_t> coef; std::vector<char> apparent; double hbits = 0; };
    std::vector<MCols> cols;

    uint32_t B(int64_t p) const { return (uint32_t)cw[p]; }
    static constexpr uint64_t CRIT = 1ULL << 63;
    uint8_t MT(int64_t p) const { uint64_t w = cw[p]; return (w & CRIT) ? 0 : (uint8_t)(w >> 32); }
    uint8_t LOC(int64_t p) const { return (uint8_t)(cw[p] >> 40); }
    uint8_t ORD(int64_t p) const { return (uint8_t)(cw[p] >> 48); }
    double val(int64_t pos) const { return vals[B(pos)]; }
    // position of the cell's mate, or -1
    int64_t mate_pos(int64_t p) const { uint8_t m = MT(p); if (!m) return -1; return p + (((m - 1) & 1) ? 1 : -1) * stride[(m - 1) >> 1]; }
    // filtration key of a cell: (birth, position of its block's top cube), then its order inside the block
    uint64_t key_hi(int64_t p) const { uint64_t w = cw[p]; return (w & 0xFFFFFFFFULL) << 32 | (uint64_t)(p - offpos[(w >> 40) & 0xFF]); }
    void decode(int64_t pos, int64_t* c) const {
        uint64_t q = (uint64_t)pos;
        for (int i = 0; i < d; ++i) { uint64_t qn = (uint64_t)(((unsigned __int128)q * recip[i]) >> 64); c[i] = (int64_t)(q - qn * (uint64_t)L[i]); q = qn; }
    }
    int32_t crit_index(int k, uint32_t pos) const {
        const auto& Bp = bypos[k];
        auto it = std::lower_bound(Bp.begin(), Bp.end(), std::make_pair(pos, (int32_t)-1));
        return it != Bp.end() && it->first == pos ? it->second : -1;
    }

    // LSD radix sort of (64-bit key, payload) pairs by key, stable, 11-bit digits, skipping constant digits
    static void radix64(std::vector<uint64_t>& key, std::vector<uint32_t>& pay, int from = 0) {
        size_t n = key.size(); if (n < 2) return;
        std::vector<uint64_t> k2(n); std::vector<uint32_t> p2(n);
        uint64_t allor = 0, alland = ~0ULL; for (uint64_t x : key) { allor |= x; alland &= x; }
        size_t cnt[2049];
        for (int sh = from; sh < 64; sh += 11) {
            if ((((allor ^ alland) >> sh) & 0x7FF) == 0) continue;
            std::fill(cnt, cnt + 2049, 0);
            for (size_t i = 0; i < n; ++i) cnt[((key[i] >> sh) & 0x7FF) + 1]++;
            for (int i = 0; i < 2048; ++i) cnt[i + 1] += cnt[i];
            for (size_t i = 0; i < n; ++i) { size_t t = cnt[(key[i] >> sh) & 0x7FF]++; k2[t] = key[i]; p2[t] = pay[i]; }
            key.swap(k2); pay.swap(p2);
        }
    }

    // births as indices into the sorted distinct values, in the low half of cw; faces take the minimum
    // over their top cofaces
    void births(const CubicalGrid& g) {
        auto T0 = std::chrono::steady_clock::now();
        auto lap = [&](const char* w) { if (!dbg) return; auto t = std::chrono::steady_clock::now(); fprintf(stderr, "      [births: %s %.3f s]\n", w, std::chrono::duration<double>(t - T0).count()); T0 = t; };
        const size_t nt = g.top.size();
        std::vector<uint32_t> rank(nt, ABSENT);            // value index of each top cube (axis 0 fastest)
        double lo = std::numeric_limits<double>::infinity(), hi = -lo; bool ints = true; size_t nf = 0;
        for (double x : g.top) if (std::isfinite(x)) { lo = std::min(lo, x); hi = std::max(hi, x); ints = ints && x == std::floor(x); ++nf; }
        if (nf && ints && hi - lo < std::max<double>(1 << 20, 4.0 * (double)nf)) {
            size_t range = (size_t)(hi - lo) + 1;
            std::vector<uint32_t> used(range, 0);
            for (double x : g.top) if (std::isfinite(x)) used[(size_t)(x - lo)] = 1;
            std::vector<uint32_t> idx(range); uint32_t m = 0;
            for (size_t r = 0; r < range; ++r) { idx[r] = m; if (used[r]) { vals.push_back(lo + (double)r); ++m; } }
            for (size_t t = 0; t < nt; ++t) { double x = g.top[t]; if (std::isfinite(x)) rank[t] = idx[(size_t)(x - lo)]; }
        } else if (nf) {
            // an order-preserving integer image of each finite double (-0 and +0 are one value), radix sorted
            auto enc = [](double x) { if (x == 0) x = 0.0; uint64_t u; std::memcpy(&u, &x, 8); return (u >> 63) ? ~u : (u | (1ULL << 63)); };
            auto dec = [](uint64_t u) { u = (u >> 63) ? (u & ~(1ULL << 63)) : ~u; double x; std::memcpy(&x, &u, 8); return x; };
            std::vector<uint64_t> key; std::vector<uint32_t> id; key.reserve(nf); id.reserve(nf);
            for (size_t t = 0; t < nt; ++t) { double x = g.top[t]; if (std::isfinite(x)) { key.push_back(enc(x)); id.push_back((uint32_t)t); } }
            lap("keys");
            radix64(key, id, 32);                           // by the high 32 bits, then runs of equal high bits
            for (size_t i = 0; i < key.size();) {
                size_t e = i + 1; while (e < key.size() && (key[e] >> 32) == (key[i] >> 32)) ++e;
                if (e - i > 1) {
                    std::vector<std::pair<uint64_t, uint32_t>> run(e - i);
                    for (size_t q = i; q < e; ++q) run[q - i] = {key[q], id[q]};
                    std::sort(run.begin(), run.end());
                    for (size_t q = i; q < e; ++q) { key[q] = run[q - i].first; id[q] = run[q - i].second; }
                }
                i = e;
            }
            lap("radix");
            uint32_t m = 0;
            for (size_t i = 0; i < key.size(); ++i) {
                if (i > 0 && key[i] != key[i - 1]) ++m;
                if (i == 0 || key[i] != key[i - 1]) vals.push_back(dec(key[i]));
                rank[id[i]] = m;
            }
        }
        lap("values");
        cw.assign(N, (uint64_t)ABSENT);
        { int64_t a[16] = {0}; int64_t pos = 0; for (int i = 0; i < d; ++i) pos += stride[i];
          for (size_t t = 0; t < nt; ++t) {
              cw[pos] = rank[t];
              for (int i = 0; i < d; ++i) { pos += 2 * stride[i]; if (++a[i] < ntx[i]) break; pos -= 2 * stride[i] * ntx[i]; a[i] = 0; }
          } }
        std::vector<uint32_t>().swap(rank);
        for (int ax = 0; ax < d; ++ax) {                        // even coordinates take the min of odd neighbours
            const int64_t st = stride[ax], La = L[ax], block = st * La;
            for (int64_t outer = 0; outer < N; outer += block)
                for (int64_t c = 0; c < La; c += 2) {
                    uint64_t* row = &cw[outer + c * st];
                    const uint64_t* pl = c - 1 >= 0 ? &cw[outer + (c - 1) * st] : nullptr;
                    const uint64_t* ph = c + 1 < La ? &cw[outer + (c + 1) * st] : nullptr;
                    if (pl && ph) for (int64_t in = 0; in < st; ++in) row[in] = std::min(pl[in], ph[in]);
                    else if (pl) for (int64_t in = 0; in < st; ++in) row[in] = pl[in];
                    else for (int64_t in = 0; in < st; ++in) row[in] = ph[in];
                }
        }
        lap("propagate");
    }

    // the gradient: for every present top cube, its block pattern, then the matching written cell by cell
    void gradient() {
        cubdmt::LocalCube lc(d);
        const int n = lc.n;
        offpos.assign(n, 0);
        for (int o = 0; o < n; ++o) for (int a = 0; a < d; ++a) offpos[o] += lc.off[o * d + a] * stride[a];
        // births of the top cubes by compact index (axis 0 fastest, the same order as positions): the
        // neighbour comparisons read this small array instead of the whole grid
        size_t nt = 1; int64_t mt[16];
        for (int i = 0; i < d; ++i) { mt[i] = (int64_t)nt; nt *= (size_t)ntx[i]; }
        std::vector<uint32_t> tb(nt);
        { int64_t a2[16] = {0}; int64_t pos = 0; for (int i = 0; i < d; ++i) pos += stride[i];
          for (size_t t = 0; t < nt; ++t) {
              tb[t] = B(pos);
              for (int i = 0; i < d; ++i) { pos += 2 * stride[i]; if (++a2[i] < ntx[i]) break; pos -= 2 * stride[i] * ntx[i]; a2[i] = 0; }
          } }
        std::vector<int> nb; std::vector<int64_t> nbt; std::vector<cubdmt::M128> nbup;
        for (int e = 0; e < n; ++e) if (e != lc.center) {
            int64_t o = 0; for (int a2 = 0; a2 < d; ++a2) o += lc.off[e * d + a2] * mt[a2];
            nb.push_back(e); nbt.push_back(o); nbup.push_back(lc.up[e]);
        }
        const size_t nn = nb.size();
        // the earlier neighbours as a bit word (bit q = neighbour q); A is the union of their shared faces,
        // assembled from one table per byte of the word
        const int nbytes = (int)((nn + 7) / 8);
        std::vector<cubdmt::M128> TT((size_t)nbytes * 256);
        for (int j = 0; j < nbytes; ++j)
            for (int v = 0; v < 256; ++v) {
                cubdmt::M128 m;
                for (int bit = 0; bit < 8; ++bit) if ((v >> bit) & 1 && 8 * j + bit < (int)nn) m |= nbup[8 * j + bit];
                TT[(size_t)j * 256 + v] = m;
            }
        std::vector<uint8_t> before(nn); for (size_t q = 0; q < nn; ++q) before[q] = nbt[q] < 0;   // earlier in position order
        cubdmt::PatternCache cache;
        // blocks in tiles (top cubes T^d at a time), so that the cells a tile writes stay in cache
        const int64_t T = d == 2 ? 64 : d == 3 ? 16 : 8;
        int64_t ntile[16], tc[16] = {0};
        size_t ntl = 1; for (int i = 0; i < d; ++i) { ntile[i] = (ntx[i] + T - 1) / T; ntl *= (size_t)ntile[i]; }
        for (size_t tl = 0; tl < ntl; ++tl) {
            int64_t lo[16], ext[16]; size_t inner = 1;
            for (int i = 0; i < d; ++i) { lo[i] = tc[i] * T; ext[i] = std::min(T, ntx[i] - lo[i]); inner *= (size_t)ext[i]; }
            int64_t c[16], t = 0, pos = 0;                  // c: grid coordinates of the top cube t
            for (int i = 0; i < d; ++i) { c[i] = lo[i]; t += lo[i] * mt[i]; pos += (2 * lo[i] + 1) * stride[i]; }
            for (size_t in = 0; in < inner; ++in) {
                uint32_t b = tb[t];
                if (b != ABSENT) {
                    bool interior = true;
                    for (int i = 0; i < d; ++i) if (c[i] < 1 || c[i] > ntx[i] - 2) { interior = false; break; }
                    uint64_t E0 = 0, E1 = 0;                     // ABSENT is the largest index: never earlier
                    if (interior) {
                        const uint32_t* tp = tb.data() + t;
                        size_t q = 0;
                        for (; q < nn && q < 64; ++q) { uint32_t b2 = tp[nbt[q]]; E0 |= (uint64_t)((b2 < b) | ((b2 == b) & before[q])) << q; }
                        for (; q < nn; ++q) { uint32_t b2 = tp[nbt[q]]; E1 |= (uint64_t)((b2 < b) | ((b2 == b) & before[q])) << (q - 64); }
                    } else {
                        for (size_t q = 0; q < nn; ++q) {
                            int e = nb[q]; bool ok = true;
                            for (int a2 = 0; a2 < d && ok; ++a2) { int64_t x = c[a2] + lc.off[e * d + a2]; if (x < 0 || x >= ntx[a2]) ok = false; }
                            if (!ok) continue;
                            uint32_t b2 = tb[t + nbt[q]];
                            uint64_t bit = (uint64_t)((b2 < b) | ((b2 == b) & before[q]));
                            if (q < 64) E0 |= bit << q; else E1 |= bit << (q - 64);
                        }
                    }
                    cubdmt::M128 A;
                    for (int j = 0; j < nbytes; ++j) { uint64_t v = j < 8 ? (E0 >> (8 * j)) & 0xFF : (E1 >> (8 * (j - 8))) & 0xFF; if (v) A |= TT[(size_t)j * 256 + v]; }
                    int32_t id = cache.get(A);
                    if (id < 0) { id = (int32_t)pats.size(); pats.push_back(lc.solve(A, pool)); cache.put(A, id); }
                    const cubdmt::Pattern& P = pats[id];
                    if (P.critical > P.minimum) ++imperfect_blocks_;
                    for (uint32_t q = P.start, e = P.start + P.len; q < e; ++q) {
                        const cubdmt::PCell& cl = pool[q];
                        cw[pos + offpos[cl.loc]] |= (uint64_t)cl.code << 32 | (uint64_t)cl.loc << 40 | (uint64_t)cl.ord << 48;
                    }
                }
                for (int i = 0; i < d; ++i) {
                    ++c[i]; t += mt[i]; pos += 2 * stride[i];
                    if (c[i] < lo[i] + ext[i]) break;
                    c[i] = lo[i]; t -= ext[i] * mt[i]; pos -= 2 * ext[i] * stride[i];
                }
            }
            for (int i = 0; i < d; ++i) { if (++tc[i] < ntile[i]) break; tc[i] = 0; }
        }
        if (dbg) fprintf(stderr, "    [morse: %zu block patterns, %lld exact searches, %lld imperfect patterns, %lld imperfect blocks]\n",
                         pats.size(), lc.exact_calls, lc.imperfect, imperfect_blocks_);
    }

    // critical cells by dimension in filtration order (radix sort on (key, order) of the cells found in
    // position order), with the inverse map
    void collect() {
        crit.assign(d + 1, {}); cnt_.assign(d + 1, 0); bypos.assign(d + 1, {}); cval.assign(d + 1, {}); ckey.assign(d + 1, {}); cord.assign(d + 1, {});
        std::vector<std::vector<uint32_t>> P(d + 1);
        { int64_t c[16] = {0}; int k = 0;
          for (int64_t pos = 0; pos < N; ++pos) {
              uint64_t w = cw[pos];
              if ((uint32_t)w != ABSENT) { ++cnt_[k]; if (((w >> 32) & 0xFF) == 0) P[k].push_back((uint32_t)pos); }
              for (int i = 0; i < d; ++i) { k -= (int)(c[i] & 1); if (++c[i] < L[i]) { k += (int)(c[i] & 1); break; } c[i] = 0; }
          } }
        for (int kk = 0; kk <= d; ++kk) {
            size_t m = P[kk].size();
            std::vector<uint64_t> key(m); std::vector<uint32_t> j(m);
            for (size_t i = 0; i < m; ++i) { key[i] = ORD(P[kk][i]); j[i] = (uint32_t)i; }
            radix64(key, j);                                        // by order inside the block, then (stable) by block
            for (size_t i = 0; i < m; ++i) key[i] = key_hi(P[kk][j[i]]);
            radix64(key, j);
            bool keys = kk >= 1 && kk <= d - 1 && d >= 3;              // keys are needed in degrees 1..d-2 (and their rows)
            crit[kk].resize(m); bypos[kk].resize(m); cval[kk].resize(m);
            if (keys) { ckey[kk].resize(m); cord[kk].resize(m); }
            for (size_t i = 0; i < m; ++i) {
                uint32_t p = P[kk][j[i]];
                crit[kk][i] = p; bypos[kk][j[i]] = {p, (int32_t)i}; cval[kk][i] = vals[key[i] >> 32];
                if (keys) { ckey[kk][i] = key[i]; cord[kk][i] = ORD(p); }
            }
            std::vector<uint32_t>().swap(P[kk]);
        }
        idx_in_cw = true;
        for (int kk = 0; kk <= d; ++kk) if (crit[kk].size() >= (1ULL << 31)) idx_in_cw = false;
        if (idx_in_cw)
            for (int kk = 0; kk <= d; ++kk)
                for (size_t i = 0; i < crit[kk].size(); ++i) { uint64_t& w = cw[crit[kk][i]]; w = (w & 0xFFFFFFFFULL) | CRIT | (uint64_t)i << 32; }
        if (dbg) { fprintf(stderr, "    [morse: critical cells"); for (int kk = 0; kk <= d; ++kk) fprintf(stderr, " %zu/%zu", crit[kk].size(), cnt_[kk]); fprintf(stderr, "]\n"); }
    }

    // H0: the two descending V-paths of each critical edge end at critical vertices (walked in position
    // order, memoised per vertex)
    void h0_ends() {
        size_t nv = 1; int64_t m[16];
        for (int i = 0; i < d; ++i) { m[i] = (int64_t)nv; nv *= (size_t)nvx[i]; }
        std::vector<int32_t> memo(nv, -1);
        int64_t c[16];
        for (auto& pi : bypos[0]) { decode(pi.first, c); int64_t v = 0; for (int i = 0; i < d; ++i) v += (c[i] / 2) * m[i]; memo[v] = pi.second; }
        std::vector<int64_t> path;
        auto term = [&](int64_t pos, int64_t v) -> int32_t {
            path.clear();
            while (memo[v] < 0) {
                uint8_t code = MT(pos);                    // a non-critical vertex is matched with an edge
                int ax = (code - 1) >> 1; int64_t s = (code - 1) & 1 ? 1 : -1;
                path.push_back(v);
                pos += 2 * s * stride[ax]; v += s * m[ax];
            }
            int32_t r = memo[v];
            for (int64_t u : path) memo[u] = r;
            return r;
        };
        ends0.resize(crit[1].size());
        for (auto& pi : bypos[1]) {
            int64_t e = pi.first; decode(e, c);
            int ax = 0; while (!(c[ax] & 1)) ++ax;
            int64_t v = 0; for (int i = 0; i < d; ++i) v += (i == ax ? (c[i] - 1) / 2 : c[i] / 2) * m[i];
            ends0[pi.second] = {term(e - stride[ax], v), term(e + stride[ax], v + m[ax])};
        }
    }

    // degree d-1: the two ascending V-paths of each critical (d-1)-cell end at critical top cells or at the
    // virtual node outside the complex (-1)
    void top_ends() {
        size_t nt = 1; int64_t m[16];
        for (int i = 0; i < d; ++i) { m[i] = (int64_t)nt; nt *= (size_t)ntx[i]; }
        const int32_t UNSET = -2, VIRT = -1;
        std::vector<int32_t> memo(nt, UNSET);
        auto tidx = [&](const int64_t* cc) { int64_t t = 0; for (int i = 0; i < d; ++i) t += ((cc[i] - 1) / 2) * m[i]; return t; };
        int64_t c[16];
        for (auto& pi : bypos[d]) { decode(pi.first, c); memo[tidx(c)] = pi.second; }
        std::vector<int64_t> path;
        auto term = [&](int64_t pos, int64_t* cc) -> int32_t {      // cc: coordinates of the top cube at pos
            int64_t t = tidx(cc); path.clear(); int32_t r;
            while (true) {
                if (memo[t] != UNSET) { r = memo[t]; break; }
                uint8_t code = MT(pos);                    // a non-critical top cell is matched with a face
                int ax = (code - 1) >> 1; int64_t s = (code - 1) & 1 ? 1 : -1;
                path.push_back(t);
                int64_t x = cc[ax] + 2 * s;
                if (x < 1 || x > L[ax] - 2 || B(pos + 2 * s * stride[ax]) == ABSENT) { r = VIRT; break; }
                cc[ax] = x; pos += 2 * s * stride[ax]; t += s * m[ax];
            }
            for (int64_t u : path) memo[u] = r;
            return r;
        };
        endsd.resize(crit[d - 1].size());
        for (auto& pi : bypos[d - 1]) {
            int64_t p0 = pi.first; decode(p0, c);
            int ax = 0; while (c[ax] & 1) ++ax;           // the one even axis
            int32_t r[2];
            for (int s = 0; s < 2; ++s) {
                int64_t cc[16]; std::memcpy(cc, c, sizeof(int64_t) * d);
                int64_t x = c[ax] + (s ? 1 : -1); int64_t p = p0 + (s ? 1 : -1) * stride[ax];
                if (x < 0 || x >= L[ax] || B(p) == ABSENT) { r[s] = VIRT; continue; }
                cc[ax] = x; r[s] = term(p, cc);
            }
            endsd[pi.second] = {r[0], r[1]};
        }
    }

    // the Morse coboundary of every critical k-cell, 1 <= k <= min(dmax, d-2): integer sums over ascending
    // V-paths, accumulated in filtration-key order (which strictly increases along every path step).
    // Columns are computed in position order (for locality) and then stored in filtration order.
    void morse_columns() {
        cols.assign(d + 1, MCols());
        int kmax = std::min(dmax, d - 2);
        // heap items: a k-cell u on an ascending path, keyed by its partner rho = V(u) (sigma by itself), with
        // rho carried along. Along a step u -> rho' -> f the key of rho' exceeds that of V(u): in one block
        // the modified Hasse diagram has rho' -> u -> V(u); otherwise rho' lies in a later block.
        struct Item { uint64_t hi; uint8_t ord; uint32_t pos, skip; };
        auto later = [](const Item& a, const Item& b) { return a.hi != b.hi ? a.hi > b.hi : a.ord > b.ord; };
        std::vector<Item> heap; cubdmt::PosMap wt, res;
        long long visited = 0; size_t maxreg = 0;
        for (int k = 1; k <= kmax; ++k) {
            MCols& M = cols[k];
            const size_t nc = crit[k].size();
            // columns in position order first
            std::vector<uint64_t> st(nc + 1, 0); std::vector<int32_t> rw; std::vector<int64_t> cf;
            std::vector<double> lognorm(nc, 0);
            std::vector<std::pair<int32_t, int64_t>> col;
            for (size_t jp = 0; jp < nc; ++jp) {
                const uint32_t sigma = bypos[k][jp].first; const int32_t sidx = bypos[k][jp].second;
                heap.clear(); wt.reset(); res.reset();
                bool ins; wt.at(sigma, ins) = 1; heap.push_back({ckey[k][sidx], cord[k][sidx], sigma, 0xFFFFFFFFu});
                size_t reg = 0;
                while (!heap.empty()) {
                    std::pop_heap(heap.begin(), heap.end(), later); Item it = heap.back(); heap.pop_back();
                    int64_t w = wt.at(it.pos, ins); ++reg;
                    if (w == 0) continue;
                    int64_t u = it.pos, cu[16]; decode(u, cu);
                    const int64_t skip = it.skip == 0xFFFFFFFFu ? -1 : (int64_t)it.skip;   // V(u), none for sigma
                    int odd = 0;
                    for (int i = 0; i < d; ++i) {
                        if (cu[i] & 1) { ++odd; continue; }
                        bool neg = odd & 1;
                        for (int s2 = 0; s2 < 2; ++s2) {
                            int64_t x = cu[i] + (s2 ? 1 : -1);
                            if (x < 0 || x >= L[i]) continue;
                            int64_t rho = u + (s2 ? 1 : -1) * stride[i];
                            if (rho == skip) continue;
                            uint64_t wr = cw[rho];
                            if ((uint32_t)wr == ABSENT) continue;
                            int64_t inc = (s2 == 0) != neg ? 1 : -1;      // [rho : u], as in the classic engine
                            if (wr & CRIT) { add_ck(res.at((uint32_t)(wr >> 32) & 0x7FFFFFFFu, ins), w * inc); continue; }
                            uint8_t mc = (uint8_t)(wr >> 32);
                            if (mc == 0) { add_ck(res.at((uint32_t)crit_index(k + 1, (uint32_t)rho), ins), w * inc); continue; }
                            int ax = (mc - 1) >> 1; int64_t sd = (mc - 1) & 1 ? 1 : -1;
                            // rho's coordinate along ax: odd means its mate is a face, even an upward mate (the path ends)
                            int64_t cax = ax == i ? x : cu[ax];
                            if (!(cax & 1)) continue;
                            int64_t f = rho + sd * stride[ax];             // rho = V(f), f a k-cell
                            int oddb = 0; for (int a = 0; a < ax; ++a) oddb += (int)((a == i ? x : cu[a]) & 1);
                            int64_t incf = ((sd > 0) != (oddb & 1)) ? 1 : -1;   // [rho : f]
                            int64_t& fw = wt.at((uint32_t)f, ins);
                            if (ins) {
                                fw = 0; add_ck(fw, -w * inc * incf);
                                uint64_t khi = (wr & 0xFFFFFFFFULL) << 32 | (uint64_t)(rho - offpos[(wr >> 40) & 0xFF]);
                                heap.push_back({khi, (uint8_t)(wr >> 48), (uint32_t)f, (uint32_t)rho}); std::push_heap(heap.begin(), heap.end(), later);
                            }
                            else add_ck(fw, -w * inc * incf);
                        }
                    }
                }
                visited += (long long)reg; maxreg = std::max(maxreg, reg);
                col.clear();
                double nrm = 0;
                for (size_t q = 0; q < res.size(); ++q) {
                    int64_t v = res.val_of(q);
                    if (v) { col.push_back({(int32_t)res.key_of(q), v}); nrm += (double)v * (double)v; }
                }
                std::sort(col.begin(), col.end());
                for (auto& e : col) { rw.push_back(e.first); cf.push_back(e.second); }
                st[jp + 1] = rw.size();
                lognorm[jp] = nrm > 0 ? 0.5 * std::log2(nrm) : 0;
            }
            // filtration order
            std::vector<uint32_t> jp_of(nc);
            for (size_t jp = 0; jp < nc; ++jp) jp_of[bypos[k][jp].second] = (uint32_t)jp;
            M.start.assign(nc + 1, 0); M.row.reserve(rw.size()); M.coef.reserve(cf.size());
            for (size_t j = 0; j < nc; ++j) {
                size_t jp = jp_of[j];
                for (uint64_t q = st[jp]; q < st[jp + 1]; ++q) { M.row.push_back(rw[q]); M.coef.push_back(cf[q]); }
                M.start[j + 1] = M.row.size();
            }
            // apparent pairs (on the integer matrix, hence in every ring): the oldest coface r of column j
            // has coefficient +-1 and j is r's youngest face
            std::vector<int32_t> yf(crit[k + 1].size(), -1);
            for (size_t j = 0; j < nc; ++j) for (uint64_t q = M.start[j]; q < M.start[j + 1]; ++q) yf[M.row[q]] = (int32_t)j;
            M.apparent.assign(nc, 0); size_t napp = 0;
            for (size_t j = 0; j < nc; ++j) {
                if (M.start[j + 1] == M.start[j]) continue;
                int64_t c0 = M.coef[M.start[j]];
                if ((c0 == 1 || c0 == -1) && yf[M.row[M.start[j]]] == (int32_t)j) { M.apparent[j] = 1; ++napp; }
            }
            // Hadamard: every minor is at most the product of its columns' norms
            size_t mm = std::min(nc, crit[k + 1].size());
            std::sort(lognorm.begin(), lognorm.end(), std::greater<double>());
            for (size_t i = 0; i < mm; ++i) M.hbits += lognorm[i];
            if (dbg) fprintf(stderr, "    [morse degree %d: %zu columns (%zu apparent), %zu entries, %lld path cells visited, largest region %zu]\n", k, nc, napp, M.row.size(), visited, maxreg);
        }
    }
    // Path sums are kept within +-path_limit() (2^62 by default), so that every product with an incidence
    // +-1 and every sum of two of them is exact in int64. Beyond it the constructor throws
    // std::overflow_error, and CubicalBackend falls back to the classic engine (or, with --engine morse,
    // reports the error). The Q reduction itself is guarded separately (Q64 overflow -> QHyb, exact).
    static void add_ck(int64_t& a, int64_t b) {
        const int64_t lim = path_limit();
        if (__builtin_add_overflow(a, b, &a) || a > lim || a < -lim) throw std::overflow_error("cubical morse: a gradient-path sum exceeds the int64 bound");
    }
public:
    // test hook: the bound on path sums (lowering it forces the overflow fallback)
    static int64_t& path_limit() { static int64_t lim = (int64_t)1 << 62; return lim; }
private:

    // zero-length bars: every k-cell is a birth in degree k or a death in degree k-1, so per value v
    //   z_k(v) = n_k(v) - (positive births_k at v) - z_{k-1}(v) - (positive deaths_{k-1} at v)
    void add_zero_length(BRun<Clear>& out, int kmax) const {
        size_t nv = vals.size();
        std::vector<std::vector<int64_t>> n(kmax + 2, std::vector<int64_t>(nv, 0));
        int64_t c[16] = {0}; int k = 0;
        for (int64_t pos = 0; pos < N; ++pos) {
            if (B(pos) != ABSENT && k <= kmax) ++n[k][B(pos)];
            for (int i = 0; i < d; ++i) { k -= (int)(c[i] & 1); if (++c[i] < L[i]) { k += (int)(c[i] & 1); break; } c[i] = 0; }
        }
        auto vidx = [&](double v) { return (size_t)(std::lower_bound(vals.begin(), vals.end(), v) - vals.begin()); };
        std::vector<int64_t> zprev(nv, 0), dprev(nv, 0);
        for (int kk = 0; kk <= kmax; ++kk) {
            std::vector<int64_t> bpos(nv, 0), dcur(nv, 0), z(nv, 0);
            for (auto& b : out.bars[kk]) { ++bpos[vidx(b.birth)]; if (std::isfinite(b.death)) ++dcur[vidx(b.death)]; }
            for (size_t v = 0; v < nv; ++v) {
                z[v] = n[kk][v] - bpos[v] - zprev[v] - dprev[v];
                for (int64_t r = 0; r < z[v]; ++r) out.bars[kk].push_back({vals[v], vals[v]});
            }
            std::sort(out.bars[kk].begin(), out.bars[kk].end());
            zprev.swap(z); dprev.swap(dcur);
        }
    }
};

} // namespace perseus
