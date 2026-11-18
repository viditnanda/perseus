// cubical_2d.hpp -- persistent homology of a 2D top-cell (T-construction) filtration by two union-finds on
// pixels, without building the cubical complex and without sorting all pixels.
//
// Order. Pixels are totally ordered by (value, index), index = x + n0*y (axis 0 fastest): the order of the
// other cubical engines. Non-finite values are absent. X_s is the union of the closed squares of the pixels
// up to s.
//
// H0. Closed squares meet iff their pixels are 8-adjacent, so H0(X) is the 8-connected components of the
// pixels present. A pixel p touches the earlier pixels among its 8 neighbours; those form c0(p) local groups
// (connected among themselves inside the 3x3 window, in which two edge neighbours at a right angle share a
// corner). c0 = 0: p starts a component. c0 = 1: p joins the component of its earlier neighbours and changes
// nothing; it gets a union-find parent pointer at once. c0 >= 2: p may merge components; it is "critical".
// Only critical pixels are sorted and processed, by the elder rule, in the filtration order.
//
// H1. By Alexander duality in S^2, H1(X_s) = H~0(S^2 - X_s), naturally in s, so the H1 barcode of X is the
// H0 barcode of the complement swept downwards. The complement of X_s is the open squares of the pixels not
// yet present, the absent pixels and the outside, glued along edges: 4-connectivity (two diagonal complement
// pixels whose other two neighbours are present are separated by the shared corner, which lies in X_s; if
// either of the others is in the complement, it connects them). Going down, p joins the complement and
// touches c1(p) local groups of its later / absent / outside neighbours (groups of the 8-ring connected
// through ring cells, containing an edge neighbour). c1 = 0: p starts a complement component (forward: p
// fills a hole, an H1 death). c1 = 1: nothing changes. c1 >= 2: components may merge (forward: p may split
// the complement, an H1 birth). The components present before any pixel enters (absent clusters and the
// outside) are born at +infinity; when two of them meet at p, the forward bar is [v(p), inf).
//
// The early/late mask of p is one byte E (bit j: ring neighbour j is earlier); the complement mask is ~E,
// because absent and outside cells are never earlier. Both local counts come from 256-entry tables.
//
// Certificate. The barcodes are those of reductions of graph incidence matrices: the 1-skeleton (degree 0)
// and the dual graph with the outside node (degree 1). Incidence matrices of graphs are totally unimodular,
// so every pivot is +-1, and the barcode is the same over Q and every F_p (docs/CERTIFICATE.md).
#pragma once
#include "io.hpp"
#include "ring.hpp"
#include "types.hpp"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <memory>
#include <cstdint>
#ifdef __linux__
#include <sys/mman.h>
#endif

namespace perseus {

namespace planar {

// ring order around a pixel: 0 E, 1 NE, 2 N, 3 NW, 4 W, 5 SW, 6 S, 7 SE (x right, y up; even = edge neighbours)
static const int RDX[8] = {1, 1, 0, -1, -1, -1, 0, 1};
static const int RDY[8] = {0, 1, 1, 1, 0, -1, -1, -1};

struct Local {
    uint8_t count = 0;           // number of local groups
    uint8_t rep[4] = {0, 0, 0, 0};   // one ring position per group
};

// The ring cell of a group nearest its circular mean: for a pixel with one group of earlier (later) neighbours,
// the middle of the arc, which approximates steepest descent (ascent) on smooth data and keeps the pointer
// chains short (the first cell in ring order would follow level sets).
inline uint8_t middle(int members) {
    double sx = 0, sy = 0;
    for (int j = 0; j < 8; ++j) if ((members >> j) & 1) { double a = j * 0.7853981633974483; sx += std::cos(a); sy += std::sin(a); }
    int best = -1; double bd = -1e9;
    for (int j = 0; j < 8; ++j) if ((members >> j) & 1) { double a = j * 0.7853981633974483, d = std::cos(a) * sx + std::sin(a) * sy; if (d > bd + 1e-12) { bd = d; best = j; } }
    return (uint8_t)best;
}

struct Tables {
    Local h0[256], h1[256];
    uint8_t mid0[256], mid1[256];   // the pointer cell when there is one group (else the first group's rep)
    Tables() {
        for (int m = 0; m < 256; ++m) {
            // H0: earlier cells, adjacent if consecutive on the ring or two edge neighbours at a right angle
            int comp[8]; for (int j = 0; j < 8; ++j) comp[j] = -1;
            int nc = 0;
            for (int j = 0; j < 8; ++j) {
                if (!((m >> j) & 1) || comp[j] >= 0) continue;
                int st[8], sn = 0; st[sn++] = j; comp[j] = nc;
                while (sn) {
                    int u = st[--sn];
                    int nb[4] = {(u + 1) & 7, (u + 7) & 7, (u % 2 == 0) ? (u + 2) & 7 : -1, (u % 2 == 0) ? (u + 6) & 7 : -1};
                    for (int v : nb) if (v >= 0 && ((m >> v) & 1) && comp[v] < 0) { comp[v] = nc; st[sn++] = v; }
                }
                h0[m].rep[nc] = (uint8_t)j; ++nc;
            }
            h0[m].count = (uint8_t)nc;
            mid0[m] = nc == 1 ? middle(m) : h0[m].rep[0];
            // H1: complement cells (mask ~m), adjacent if consecutive on the ring; count groups that contain an
            // edge neighbour, represented by one
            int c = (~m) & 0xFF;
            for (int j = 0; j < 8; ++j) comp[j] = -1;
            int g = 0, n1 = 0, members = 0;
            for (int j = 0; j < 8; ++j) {
                if (!((c >> j) & 1) || comp[j] >= 0) continue;
                int st[8], sn = 0; st[sn++] = j; comp[j] = g;
                int edge = -1;
                while (sn) {
                    int u = st[--sn];
                    if (u % 2 == 0 && (edge < 0 || u < edge)) edge = u;
                    for (int v : {(u + 1) & 7, (u + 7) & 7}) if (((c >> v) & 1) && comp[v] < 0) { comp[v] = g; st[sn++] = v; }
                }
                if (edge >= 0) { h1[m].rep[n1++] = (uint8_t)edge; for (int j2 = 0; j2 < 8; ++j2) if (comp[j2] == g) members |= 1 << j2; }
                ++g;
            }
            h1[m].count = (uint8_t)n1;
            mid1[m] = n1 == 1 ? middle(members) : h1[m].rep[0];
        }
    }
};
inline const Tables& tables() { static const Tables T; return T; }

// The earlier-neighbour byte of interior pixels from the rows below (rm), at (r0) and above (rq): bit j is set
// if ring neighbour j comes before the pixel in the order (value, index). NaN and +inf never compare as earlier;
// ties go to the neighbours of smaller index (W, SW, S, SE). erow4 does four pixels with vector compares.
inline uint8_t erow1(const double* r0, const double* rm, const double* rq, int64_t x) {
    const double vp = r0[x];
    return (uint8_t)((unsigned)(r0[x + 1] < vp) | (unsigned)(rq[x + 1] < vp) << 1 | (unsigned)(rq[x] < vp) << 2
            | (unsigned)(rq[x - 1] < vp) << 3 | (unsigned)(r0[x - 1] <= vp) << 4 | (unsigned)(rm[x - 1] <= vp) << 5
            | (unsigned)(rm[x] <= vp) << 6 | (unsigned)(rm[x + 1] <= vp) << 7);
}
#if defined(__GNUC__)
typedef double d2v __attribute__((vector_size(16)));      // SSE2 / NEON width: no ABI change without AVX
typedef long long l2v __attribute__((vector_size(16)));
inline d2v ld2(const double* p) { d2v r; std::memcpy(&r, p, sizeof r); return r; }
inline l2v erow2(const double* r0, const double* rm, const double* rq, int64_t x) {
    const d2v c = ld2(r0 + x);
    return ((ld2(r0 + x + 1) < c) & 1) | ((ld2(rq + x + 1) < c) & 2) | ((ld2(rq + x) < c) & 4) | ((ld2(rq + x - 1) < c) & 8)
         | ((ld2(r0 + x - 1) <= c) & 16) | ((ld2(rm + x - 1) <= c) & 32) | ((ld2(rm + x) <= c) & 64) | ((ld2(rm + x + 1) <= c) & 128);
}
inline void erow4(const double* r0, const double* rm, const double* rq, int64_t x, uint8_t* out) {
    const l2v a = erow2(r0, rm, rq, x), b = erow2(r0, rm, rq, x + 2);
    out[0] = (uint8_t)a[0]; out[1] = (uint8_t)a[1]; out[2] = (uint8_t)b[0]; out[3] = (uint8_t)b[1];
}
#else
inline void erow4(const double* r0, const double* rm, const double* rq, int64_t x, uint8_t* out) {
    for (int k = 0; k < 4; ++k) out[k] = erow1(r0, rm, rq, x + k);
}
#endif

// large uninitialised buffers for the random-access arrays, on transparent huge pages where the kernel offers
// them on request (Linux): the union-find touches them at random, and 4 KiB pages make every access a TLB miss
struct FreeDeleter { void operator()(void* p) const { std::free(p); } };
template <class T> using Buf = std::unique_ptr<T[], FreeDeleter>;
template <class T> Buf<T> big(size_t n) {
    const size_t HP = size_t(1) << 21;
    size_t bytes = std::max<size_t>(n, 1) * sizeof(T);
    void* p = nullptr;
#if defined(__linux__) && defined(MADV_HUGEPAGE)
    if (bytes >= HP) {
        bytes = (bytes + HP - 1) / HP * HP;
        p = std::aligned_alloc(HP, bytes);
        if (p) madvise(p, bytes, MADV_HUGEPAGE);
    }
#else
    (void)HP;
#endif
    if (!p) p = std::malloc(bytes);
    if (!p) throw std::bad_alloc();
    return Buf<T>(static_cast<T*>(p));
}

// an order-preserving integer image of a finite double (-0 and +0 are one value)
inline uint64_t enc(double x) {
    x = x + 0.0;                                     // -0 + 0 = +0
    uint64_t u; std::memcpy(&u, &x, 8);
    return u ^ ((uint64_t)((int64_t)u >> 63) | (1ULL << 63));
}

// stable LSD radix sort of n (key, payload) pairs by the key bits from `from` up, 11-bit digits, skipping
// constant digits, alternating between (key, pay) and (k2, p2); returns true if the result is in (k2, p2)
inline bool radix(uint64_t* key, uint32_t* pay, uint64_t* k2, uint32_t* p2, size_t n, int from) {
    bool sw = false; if (n < 2) return sw;
    uint64_t allor = 0, alland = ~0ULL; for (size_t i = 0; i < n; ++i) { allor |= key[i]; alland &= key[i]; }
    std::vector<size_t> cnt(2049);
    for (int sh = from; sh < 64; sh += 11) {
        if ((((allor ^ alland) >> sh) & 0x7FF) == 0) continue;
        std::fill(cnt.begin(), cnt.end(), 0);
        for (size_t i = 0; i < n; ++i) cnt[((key[i] >> sh) & 0x7FF) + 1]++;
        for (int i = 0; i < 2048; ++i) cnt[i + 1] += cnt[i];
        for (size_t i = 0; i < n; ++i) { size_t t = cnt[(key[i] >> sh) & 0x7FF]++; k2[t] = key[i]; p2[t] = pay[i]; }
        std::swap(key, k2); std::swap(pay, p2); sw = !sw;
    }
    return sw;
}
inline double dec(uint64_t u) { u = (u >> 63) ? (u & ~(1ULL << 63)) : ~u; double x; std::memcpy(&x, &u, 8); return x; }
// sort n (key, payload) pairs by (key, payload & mask), given in increasing payload order: a stable radix sort
// of the high halves of the keys, then the short runs of equal high halves; returns true if the result is in
// (k2, p2)
inline bool sort_keys(uint64_t* key, uint32_t* pay, uint64_t* k2, uint32_t* p2, size_t n, uint32_t mask) {
    const bool sw = radix(key, pay, k2, p2, n, 32);
    if (sw) { std::swap(key, k2); std::swap(pay, p2); }
    std::vector<std::pair<uint64_t, uint32_t>> run;
    for (size_t i = 0; i < n;) {
        size_t e = i + 1; bool same = true;
        while (e < n && (key[e] >> 32) == (key[i] >> 32)) { same &= key[e] == key[i]; ++e; }
        if (same) { i = e; continue; }                 // equal keys: already in payload order (stable radix)
        if (e - i <= 16) {
            for (size_t q = i + 1; q < e; ++q) {
                uint64_t k = key[q]; uint32_t jj = pay[q]; size_t t = q;
                while (t > i && (key[t - 1] > k || (key[t - 1] == k && (pay[t - 1] & mask) > (jj & mask)))) { key[t] = key[t - 1]; pay[t] = pay[t - 1]; --t; }
                key[t] = k; pay[t] = jj;
            }
        } else if (e - i > 16) {
            run.resize(e - i);
            for (size_t q = i; q < e; ++q) run[q - i] = {key[q], pay[q]};
            std::sort(run.begin(), run.end(), [mask](const std::pair<uint64_t, uint32_t>& a, const std::pair<uint64_t, uint32_t>& b) {
                return a.first < b.first || (a.first == b.first && (a.second & mask) < (b.second & mask)); });
            for (size_t q = i; q < e; ++q) { key[q] = run[q - i].first; pay[q] = run[q - i].second; }
        }
        i = e;
    }
    return sw;
}

} // namespace planar

class CubicalPlanar {
public:
    typedef std::vector<uint32_t> Clear;     // unused: the degrees are independent
    bool keep_zero_length = false;           // must be passed to the constructor (zero_length) to take effect

    static bool supports(const CubicalGrid& g) {
        if (g.ext.size() != 2) return false;
        for (size_t i = 0; i < g.periodic.size(); ++i) if (g.periodic[i]) return false;
        return g.ext[0] > 0 && g.ext[1] > 0 && (uint64_t)(g.ext[0] + 2) * (uint64_t)(g.ext[1] + 2) < (1ULL << 31);   // padded indices in 31 bits
    }

    CubicalPlanar(const CubicalGrid& g, int maxdim, bool zero_length = false) {
        if (!supports(g)) throw std::runtime_error("cubical planar: needs a 2D grid without periodic axes");
        dmax = std::min(maxdim < 0 ? 2 : maxdim, 2);
        keep_zero_length = zero_length; zero_ready = zero_length;
        bool dbg = getenv("PERSEUS_DEBUG") != nullptr;
        auto T0 = std::chrono::steady_clock::now();
        auto lap = [&](const char* w) { if (!dbg) return; auto t = std::chrono::steady_clock::now(); fprintf(stderr, "    [planar %s %.3f s]\n", w, std::chrono::duration<double>(t - T0).count()); T0 = t; };
        const int64_t n0 = g.ext[0], n1 = g.ext[1], W = n0 + 2;   // padded width: a frame of outside pixels
        const int64_t NP = W * (n1 + 2);
        const double* v = g.data();
        const planar::Tables& TB = planar::tables();
        // rows holding -inf (absent, but it would compare as earlier in the branch-free path), found one row ahead
        auto ninf_row = [&](int64_t y) { if (y < 0 || y >= n1) return false; const double* r = v + y * n0; int64_t c = 0;
            for (int64_t x = 0; x < n0; ++x) c += r[x] == -std::numeric_limits<double>::infinity();
            return c > 0; };
        bool ninf_prev = false, ninf_cur = ninf_row(0), ninf_next = false;
        // --- streaming classification ------------------------------------------------------------------------
        // Per present pixel p (padded index pp): the earlier-neighbour byte Eb[pp] and the H0 pointer P[pp]: if
        // the earlier neighbours of p form one local group, P[pp] points at one of them (p joins its component
        // and changes nothing). The pixels that must be processed in the filtration order (H0 births, H0 merge
        // candidates, H1 merge candidates) are collected with their value keys, in index order.
        planar::Buf<uint8_t> Ebuf = planar::big<uint8_t>((size_t)NP); uint8_t* Eb = Ebuf.get();   // set for present pixels
        planar::Buf<int32_t> Pbuf = planar::big<int32_t>((size_t)NP); int32_t* P = Pbuf.get();    // union-find, H1 then H0
        int32_t off[8]; for (int j = 0; j < 8; ++j) off[j] = (int32_t)(planar::RDX[j] + planar::RDY[j] * W);
        int64_t ooff[8]; for (int j = 0; j < 8; ++j) ooff[j] = planar::RDX[j] + planar::RDY[j] * n0;
        bool any_absent = false;
        // inS: births (c0 = 0) and merge candidates (c0 >= 2, equivalently c1 >= 2). The H1 pointer of a pixel
        // points at a later neighbour in one of its complement groups (for a candidate, until it is processed:
        // every group lies in its final component); a hole fill is a root
        struct Cls { int32_t d1, r1; uint32_t cflag; uint8_t inS, fill, special; };
        Cls cls[256]; int32_t d0[256];
        for (int m = 0; m < 256; ++m) {
            Cls& c = cls[m];
            c.inS = (uint8_t)(TB.h0[m].count != 1 || TB.h1[m].count >= 2);
            c.cflag = TB.h1[m].count >= 2 ? 0x80000000u : 0u;
            c.d1 = TB.h1[m].count >= 1 ? off[TB.mid1[m]] : 0; c.r1 = TB.h1[m].count >= 1 ? 0 : -1;
            c.fill = dmax >= 1 && TB.h1[m].count == 0;
            c.special = c.inS | c.fill;
            d0[m] = TB.h0[m].count >= 1 ? off[TB.mid0[m]] : 0;
        }
        // the hole fills (c1 = 0) are the finite H1 roots: P = -4 - (their number), values kept compactly in fv
        planar::Buf<double> fvbuf = planar::big<double>((size_t)(n0 * n1) + 1); double* fv = fvbuf.get();
        int32_t nfill = 0;
        // the processed pixels: value keys and padded indices (bit 31: merge candidate), appended without branches
        planar::Buf<uint64_t> keyb = planar::big<uint64_t>((size_t)(n0 * n1) + 1);
        planar::Buf<uint32_t> sppb = planar::big<uint32_t>((size_t)(n0 * n1) + 1);
        uint64_t* key = keyb.get(); uint32_t* spp = sppb.get();
        size_t ns = 0;
        bool sparse = false;                              // previous row had few special pixels: branch on them
        for (int64_t y = 0; y < n1; ++y) {
            const double* r0 = v + y * n0;
            const int32_t rowp = (int32_t)(1 + (y + 1) * W);
            auto emit = [&](int64_t x, double vp, unsigned E) {
                const int32_t pp = rowp + (int32_t)x; const Cls& c = cls[E];
                const int32_t isf = c.fill;
                const int32_t msk = -isf;                    // without a branch: a fill is a root, else a pointer
                Eb[pp] = (uint8_t)E; fv[nfill] = vp; P[pp] = ((-4 - nfill) & msk) | ((c.r1 | (pp + c.d1)) & ~msk); nfill += isf;
                key[ns] = planar::enc(vp); spp[ns] = (uint32_t)pp | c.cflag; ns += c.inS;
            };
            auto edge = [&](int64_t x) {                   // a pixel on the border of the grid, or any pixel near -inf
                const int64_t o = x + y * n0;
                const double vp = r0[x];
                if (!std::isfinite(vp)) { any_absent = true; Eb[rowp + (int32_t)x] = 0; return; }
                unsigned E = 0;
                for (int j = 0; j < 8; ++j) {
                    int64_t xx = x + planar::RDX[j], yy = y + planar::RDY[j];
                    if (xx < 0 || yy < 0 || xx >= n0 || yy >= n1) continue;
                    double vq = v[o + ooff[j]];
                    if (std::isfinite(vq) && ((vq < vp) || (vq == vp && ooff[j] < 0))) E |= 1u << j;
                }
                emit(x, vp, E);
            };
            ninf_next = ninf_row(y + 1);
            const bool fast = y > 0 && y + 1 < n1 && !ninf_prev && !ninf_cur && !ninf_next && n0 > 2;
            ninf_prev = ninf_cur; ninf_cur = ninf_next;
            if (fast) {
                // NaN and +inf never compare as earlier; ties go to the earlier index (W, SW, S, SE come first)
                const double *rm = r0 - n0, *rq = r0 + n0;
                uint8_t* Er = Eb + rowp;                     // the row's bytes, written in place
                int64_t x = 1;
                for (; x + 4 < n0; x += 4) planar::erow4(r0, rm, rq, x, Er + x);   // 4 pixels at a time
                for (; x + 1 < n0; ++x) Er[x] = planar::erow1(r0, rm, rq, x);
                edge(0);
                const size_t before = ns + (size_t)nfill;
                if (sparse) {
                    for (x = 1; x + 1 < n0; ++x) {
                        const double vp = r0[x];
                        if (!std::isfinite(vp)) { any_absent = true; Er[x] = 0; continue; }
                        const Cls& c = cls[Er[x]];
                        if (__builtin_expect(c.special, 0)) emit(x, vp, Er[x]);
                        else P[rowp + (int32_t)x] = rowp + (int32_t)x + c.d1;
                    }
                } else {
                    for (x = 1; x + 1 < n0; ++x) {
                        const double vp = r0[x];
                        if (!std::isfinite(vp)) { any_absent = true; Er[x] = 0; continue; }
                        emit(x, vp, Er[x]);
                    }
                }
                sparse = (ns + (size_t)nfill - before) * 8 < (size_t)n0;
                edge(n0 - 1);
            } else for (int64_t x = 0; x < n0; ++x) edge(x);
        }
        lap("classification");
        size_t ncand = 0; for (size_t si = 0; si < ns; ++si) ncand += spp[si] >> 31;
        // Compression pays when chains are long and candidates sparse (smooth images); in noise-like images the
        // chains are short and the finds in the processing loops are cheaper than a separate pass.
        bool precompress = ncand * 16 < (size_t)(n0 * n1) && n0 * n1 >= (1 << 21);   // and P outgrows the caches
        if (const char* pc = getenv("PERSEUS_PLANAR_COMPRESS")) precompress = pc[0] == '1';   // for testing both paths
        if (dbg) fprintf(stderr, "    [planar %zu processed, %zu candidates, %d hole fills, of %lld pixels]\n", ns, ncand, (int)nfill, (long long)(n0 * n1));
        auto find = [P](int32_t x) {       // path halving; roots are negative
            int32_t p;
            while ((p = P[x]) >= 0) { int32_t q = P[p]; if (q < 0) return p; P[x] = q; x = q; }
            return x;
        };
        // Full path compression from the neighbours of the candidates, in index order, before the processing in
        // the filtration order: the pointer chains are spatially local, so this is cache friendly, and each pixel
        // on a chain is compressed once; afterwards a find costs one step until components merge.
        std::vector<int32_t> path(64);
        auto compress_from = [&](int32_t q) {
            if (P[q] < 0) return;
            size_t np = 0;
            while (P[P[q]] >= 0) { path[np++] = q; q = P[q]; if (np == path.size()) path.resize(2 * np); }
            const int32_t r = P[q];
            for (size_t i = 0; i < np; ++i) P[path[i]] = r;
        };
        // H1 setup: the outside V (all frame cells, root at padded index 0) and the absent clusters
        if (dmax >= 1) {
            for (int64_t x = 0; x < W; ++x) { P[x] = 0; P[x + (n1 + 1) * W] = 0; }
            for (int64_t y = 0; y < n1 + 2; ++y) { P[y * W] = 0; P[y * W + n0 + 1] = 0; }
            P[0] = -3;
            if (any_absent) {                                // absent clusters (4-connected); the border joins them to the outside
                for (int64_t y = 0; y < n1; ++y) for (int64_t x = 0; x < n0; ++x) {
                    if (std::isfinite(v[x + y * n0])) continue;
                    const int32_t pp = (int32_t)((x + 1) + (y + 1) * W);
                    int32_t a = (x > 0 && !std::isfinite(v[x - 1 + y * n0])) ? find(pp - 1) : -1;
                    int32_t b = (y > 0 && !std::isfinite(v[x + (y - 1) * n0])) ? find(pp - (int32_t)W) : -1;
                    int32_t r;
                    if (a < 0 && b < 0) { P[pp] = -2; r = pp; }
                    else if (a < 0 || b < 0) { r = a < 0 ? b : a; P[pp] = r; }
                    else { if (a != b) { if (a == 0) std::swap(a, b); P[a] = b; } r = b; P[pp] = r; }
                    if ((x == 0 || y == 0 || x + 1 == n0 || y + 1 == n1) && r != 0) P[r] = 0;
                }
            }
            if (precompress) for (size_t si = 0; si < ns; ++si) {
                const uint32_t u = spp[si]; if (!(u >> 31)) continue;
                const int32_t pp = (int32_t)(u & 0x7FFFFFFFu); const planar::Local& b = TB.h1[Eb[pp]];
                for (int k = 0; k < b.count; ++k) compress_from(pp + off[b.rep[k]]);
            }
            lap("H1 setup");
        }
        {
            planar::Buf<uint64_t> k2 = planar::big<uint64_t>(ns); planar::Buf<uint32_t> p2 = planar::big<uint32_t>(ns);
            if (planar::sort_keys(key, spp, k2.get(), p2.get(), ns, 0x7FFFFFFFu)) { keyb.swap(k2); sppb.swap(p2); key = keyb.get(); spp = sppb.get(); }
        }                                              // (value, index): padded indices are in index order
        auto val = [&](size_t si) { return planar::dec(key[si]); };
        lap("sort");
        bars.assign(3, {});
        const double INF = std::numeric_limits<double>::infinity();
        const size_t PD = 12;                            // prefetch distance in the processing loops
        const bool PF = NP > (1 << 19);                  // prefetch only when P (4 bytes a cell) outgrows the L2 cache
        auto fix_runs = [](std::vector<Bar>& B) {        // equal births: order by death
            std::vector<uint64_t> dk, dk2; std::vector<uint32_t> pz, pz2;
            for (size_t i = 0; i < B.size();) {
                size_t e = i + 1; while (e < B.size() && B[e].birth == B[i].birth) ++e;
                const size_t n = e - i;
                if (n > 1 && n <= 32) {                  // insertion sort on the deaths
                    for (size_t q = i + 1; q < e; ++q) { double d = B[q].death; size_t t = q; while (t > i && B[t - 1].death > d) { B[t].death = B[t - 1].death; --t; } B[t].death = d; }
                } else if (n > 32) {                     // radix sort of the encoded deaths (+inf is the largest)
                    dk.resize(n); dk2.resize(n); pz.assign(n, 0); pz2.resize(n);
                    for (size_t q = 0; q < n; ++q) dk[q] = planar::enc(B[i + q].death);
                    if (planar::radix(dk.data(), pz.data(), dk2.data(), pz2.data(), n, 0)) dk.swap(dk2);
                    for (size_t q = 0; q < n; ++q) B[i + q].death = planar::dec(dk[q]);
                }
                i = e;
            }
        };
        // --- H1 first: going down through the complement. Roots: -4 - f a finite component (born at hole fill f,
        // value fv[f]), -2 an absent cluster, -3 the outside (at padded index 0); the last two are born at +infinity.
        // A candidate p with c local groups of each kind gives r1 - 1 H1 births (r1 = distinct complement
        // components it meets) and c - r1 H0 merges, since adding p changes the Euler characteristic by 1 - c.
        // When c = r1, p is marked and the H0 pass treats it as a non-critical pixel.
        const uint32_t SKIP = 0xFFFFFFFFu, PIX = 0x7FFFFFFFu;
        if (dmax >= 1) {
            auto kind = [](int32_t c) { return c <= -4 ? -1 : c; };   // -3 outside, -2 absent, -1 finite: eldest first
            bars[1].reserve(ncand + 1);
            for (size_t si = ns; si-- > 0;) {                // merge candidates in reverse order
                if (PF && si >= 2 * PD) { const uint32_t u = spp[si - 2 * PD]; if (u >> 31) { const int32_t q = (int32_t)(u & PIX); __builtin_prefetch(Eb + q); __builtin_prefetch(P + q - W); __builtin_prefetch(P + q); __builtin_prefetch(P + q + W); } }
                if (PF && si >= PD) { const uint32_t u = spp[si - PD]; if (u >> 31) {   // second stage: the pointer targets
                    const int32_t q = (int32_t)(u & PIX); const planar::Local& bq = TB.h1[Eb[q]];
                    for (int k = 0; k < bq.count; ++k) { const int32_t t = P[q + off[bq.rep[k]]]; if (t >= 0) __builtin_prefetch(P + t); else if (t <= -4) __builtin_prefetch(fv + (-4 - t)); } } }
                const uint32_t u = spp[si];
                if (!(u >> 31)) continue;                    // a birth
                const int32_t pp = (int32_t)(u & PIX);
                const uint8_t E = Eb[pp];
                const planar::Local& b = TB.h1[E];
                const double vb = val(si);
                if (b.count == 2) {                          // the common case, without the general bookkeeping
                    int32_t x = find(pp + off[b.rep[0]]), y = find(pp + off[b.rep[1]]);
                    if (x == y) { P[pp] = x; continue; }
                    const int32_t cx = P[x], cy = P[y];
                    if (cx <= -4 && cy <= -4) {
                        double fx = fv[-4 - cx], fy = fv[-4 - cy];
                        if (fx > fy || (fx == fy && x > y)) { std::swap(x, y); std::swap(fx, fy); }   // y: the eldest (latest)
                        if (fx > vb) bars[1].push_back({vb, fx});
                    } else {
                        if (kind(cx) < kind(cy)) std::swap(x, y);   // y: the eldest (outside, then absent)
                        const int32_t c = P[x];
                        if (c > -4) bars[1].push_back({vb, INF});
                        else { const double dx = fv[-4 - c]; if (dx > vb) bars[1].push_back({vb, dx}); }
                    }
                    P[x] = y; P[pp] = y;
                    if (TB.h0[E].count == 2) spp[si] = SKIP;
                    continue;
                }
                int32_t rr[4]; double rv[4]; int nr = 0, e = 0;
                for (int k = 0; k < b.count; ++k) {
                    int32_t x = find(pp + off[b.rep[k]]); bool dup = false;
                    for (int t = 0; t < nr; ++t) dup |= rr[t] == x;
                    if (dup) continue;
                    rr[nr] = x; const int32_t cx = kind(P[x]); rv[nr] = cx == -1 ? fv[-4 - P[x]] : INF;
                    if (nr > 0) {                            // the eldest going down: the outside, an absent cluster, the latest pixel
                        const int32_t y = rr[e]; const int32_t cy = kind(P[y]);
                        if (cx != cy ? cx < cy : (cx == -1 && (rv[nr] > rv[e] || (rv[nr] == rv[e] && x > y)))) e = nr;
                    }
                    ++nr;
                }
                for (int t = 0; t < nr; ++t) {
                    const int32_t x = rr[t]; if (t == e) continue;
                    if (rv[t] == INF) bars[1].push_back({vb, INF});
                    else if (rv[t] > vb) bars[1].push_back({vb, rv[t]});
                    P[x] = rr[e];
                }
                P[pp] = rr[e];
                if (nr == TB.h0[E].count) spp[si] = SKIP;
            }
            std::reverse(bars[1].begin(), bars[1].end());
            fix_runs(bars[1]);
            lap("H1");
        }
        // --- H0: the elder rule; a root is a birth, stored as -1 - (its rank). P is reset to the descending
        // pointers (one earlier neighbour; candidates and births are set when processed) ------------------------
        {
            int32_t e0[256]; for (int m = 0; m < 256; ++m) e0[m] = TB.h0[m].count ? d0[m] : INT32_MIN;
            for (int64_t y = 0; y < n1; ++y) {               // absent pixels (Eb = 0) are never reached from present ones
                const int32_t rowp = (int32_t)(1 + (y + 1) * W);
                for (int32_t x = 0; x < (int32_t)n0; ++x) {
                    const int32_t pp = rowp + x; const int32_t d = e0[Eb[pp]];
                    P[pp] = d == INT32_MIN ? INT32_MIN : pp + d;   // a birth: set when processed
                }
            }
            if (precompress) for (int64_t y = 0; y < n1; ++y) {
                const int32_t rowp = (int32_t)(1 + (y + 1) * W);
                for (int64_t x = 0; x < n0; ++x) {
                    const int32_t pp = rowp + (int32_t)x; const planar::Local& a = TB.h0[Eb[pp]];   // absent: Eb = 0
                    if (a.count >= 2) for (int k = 0; k < a.count; ++k) compress_from(pp + off[a.rep[k]]);
                }
            }
            lap("H0 pointers");
            planar::Buf<int32_t> diesb = planar::big<int32_t>(ns); int32_t* dies = diesb.get();   // rank of a birth -> rank of its killer
            std::fill(dies, dies + ns, -2);
            bars[0].reserve(ns - ncand + 1);
            for (size_t si = 0; si < ns; ++si) {
                if (PF && si + 2 * PD < ns) { const uint32_t u = spp[si + 2 * PD]; if (u != SKIP) { const int32_t q = (int32_t)(u & PIX); __builtin_prefetch(Eb + q); __builtin_prefetch(P + q - W); __builtin_prefetch(P + q); __builtin_prefetch(P + q + W); } }
                if (PF && si + PD < ns) { const uint32_t u = spp[si + PD]; if (u != SKIP) {   // second stage: the pointer targets
                    const int32_t q = (int32_t)(u & PIX); const planar::Local& aq = TB.h0[Eb[q]];
                    for (int k = 0; k < aq.count; ++k) { const int32_t t = P[q + off[aq.rep[k]]]; if (t >= 0) __builtin_prefetch(P + t); } } }
                const uint32_t u = spp[si];
                if (u == SKIP) continue;
                const int32_t pp = (int32_t)(u & PIX);
                const planar::Local& a = TB.h0[Eb[pp]];
                if (a.count == 0) { P[pp] = -1 - (int32_t)si; dies[si] = -1; continue; }
                int32_t rr[4], rk[4]; int nr = 0, e = 0;
                for (int k = 0; k < a.count; ++k) {
                    int32_t x = find(pp + off[a.rep[k]]); bool dup = false;
                    for (int t = 0; t < nr; ++t) dup |= rr[t] == x;
                    if (!dup) { rr[nr] = x; rk[nr] = -1 - P[x]; if (rk[nr] < rk[e]) e = nr; ++nr; }
                }
                for (int t = 0; t < nr; ++t) if (t != e) { P[rr[t]] = rr[e]; dies[rk[t]] = (int32_t)si; }
                P[pp] = rr[e];
            }
            for (size_t si = 0; si < ns; ++si) {             // bars in birth order
                if (dies[si] == -1) bars[0].push_back({val(si), INF});
                else if (dies[si] >= 0) { double b = val(si), d = val((size_t)dies[si]); if (d > b) bars[0].push_back({b, d}); }
            }
            fix_runs(bars[0]);
            lap("H0");
        }
        if (zero_length) add_zero_length(g);
    }

    int maxdim() const { return dmax; }
    bool unimodular() const { return true; }

    template <class R>
    BRun<Clear> run(const R& ring, int kmin, int kmax, const Clear*) const {
        if (keep_zero_length && !zero_ready) throw std::logic_error("cubical planar: zero-length bars must be requested at construction");
        BRun<Clear> out; out.kmin = kmin;
        kmax = std::min(kmax, dmax);
        int nd = std::max(0, kmax + 1);
        out.bars.assign(nd, {}); out.clear_out.assign(nd, {});
        for (int l = 0; l < R::lanes; ++l) out.moduli.push_back(ring.modulus(l));
        for (int k = std::max(kmin, 0); k <= kmax; ++k) {
            out.bars[k] = bars[k];
            if (keep_zero_length) { out.bars[k].insert(out.bars[k].end(), zero[k].begin(), zero[k].end()); std::sort(out.bars[k].begin(), out.bars[k].end()); }
        }
        return out;
    }

private:
    int dmax = 2; bool zero_ready = false;
    Barcode bars, zero;

    // zero-length bars: every k-cell is a birth in degree k or a death in degree k-1, so per value v
    //   z_k(v) = n_k(v) - (positive births_k at v) - z_{k-1}(v) - (positive deaths_{k-1} at v);
    // cells of the doubled grid take the minimum value of their top cofaces
    void add_zero_length(const CubicalGrid& g) {
        const int64_t n0 = g.ext[0], n1 = g.ext[1];
        const double* v = g.data();
        std::vector<double> vals; for (int64_t i = 0; i < n0 * n1; ++i) if (std::isfinite(v[i])) vals.push_back(v[i] == 0 ? 0.0 : v[i]);
        std::sort(vals.begin(), vals.end()); vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
        size_t nv = vals.size();
        auto idx = [&](double x) { return (size_t)(std::lower_bound(vals.begin(), vals.end(), x == 0 ? 0.0 : x) - vals.begin()); };
        std::vector<std::vector<int64_t>> n(3, std::vector<int64_t>(nv, 0));
        const double INF = std::numeric_limits<double>::infinity();
        auto at = [&](int64_t x, int64_t y) { return (x < 0 || y < 0 || x >= n0 || y >= n1 || !std::isfinite(v[x + y * n0])) ? INF : v[x + y * n0]; };
        for (int64_t Y = 0; Y <= 2 * n1; ++Y) for (int64_t X = 0; X <= 2 * n0; ++X) {
            double m = INF;   // min over the pixels (top cofaces) of the cell with doubled coordinates (X, Y)
            for (int64_t dx : {-1, 1}) for (int64_t dy : {-1, 1}) {
                int64_t px = (X & 1) ? (X - 1) / 2 : (dx < 0 ? X / 2 - 1 : X / 2);
                int64_t py = (Y & 1) ? (Y - 1) / 2 : (dy < 0 ? Y / 2 - 1 : Y / 2);
                m = std::min(m, at(px, py));
            }
            if (!std::isfinite(m)) continue;
            ++n[(X & 1) + (Y & 1)][idx(m)];
        }
        zero.assign(3, {});
        std::vector<int64_t> zprev(nv, 0), dprev(nv, 0);
        for (int k = 0; k <= dmax; ++k) {
            std::vector<int64_t> bpos(nv, 0), dcur(nv, 0);
            for (auto& b : bars[k]) { ++bpos[idx(b.birth)]; if (std::isfinite(b.death)) ++dcur[idx(b.death)]; }
            std::vector<int64_t> z(nv, 0);
            for (size_t i = 0; i < nv; ++i) {
                z[i] = n[k][i] - bpos[i] - zprev[i] - dprev[i];
                for (int64_t r = 0; r < z[i]; ++r) zero[k].push_back({vals[i], vals[i]});
            }
            zprev.swap(z); dprev.swap(dcur);
        }
    }
};

} // namespace perseus
