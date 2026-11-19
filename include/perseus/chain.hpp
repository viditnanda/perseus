// chain.hpp -- filtered based chain complexes over Z: validation, algebraic Morse reduction, and the
// integral outputs (torsion of H_*(K_inf), and the image of H_*(K_s; Z) -> H_*(K_t; Z)).
//
// Algebraic Morse reduction (Forman; Skoldberg; Joellenbeck-Welker), done one pair at a time. A pair
// (sigma < tau) with [tau : sigma] = +-1 and birth(sigma) = birth(tau) is cancelled by the change of basis
//     rho' = rho - [rho : sigma][tau : sigma] tau        for the other cofaces rho of sigma,
//     sigma' = d tau                                     (unit coefficient on sigma).
// Both are filtration-preserving (every term has birth <= that of the cell it replaces: birth(tau) =
// birth(sigma) <= birth(rho), and d tau is born no later than tau), so the complex splits, as a filtered
// complex over Z, into the residual complex plus the acyclic piece tau -> d tau. Consequently every
// homology group H_*(K_s; Z), every map H_*(K_s) -> H_*(K_t), and every barcode over every field are
// those of the residual complex; the pair contributes one zero-length bar [b, b) in degree dim sigma.
// In the residual complex d rho' = d rho - [rho:sigma][tau:sigma] d tau, and tau drops out of every
// boundary (its coefficient there becomes 0 because d d = 0).
#pragma once
#include "complex.hpp"
#include "zlin.hpp"
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <queue>
#include <stdexcept>

namespace perseus {

// faces exist, have dimension d-1, are born no later; d d = 0 exactly (128-bit accumulation, overflow checked)
inline std::string validate_chain_complex(const FilteredComplex& C) {
    std::string e = C.validate(false);
    if (!e.empty()) return e;
    std::unordered_map<int32_t, __int128> acc;
    for (size_t c = 0; c < C.size(); ++c) {
        acc.clear();
        bool over = false;
        for (int64_t a = C.off[c]; a < C.off[c + 1] && !over; ++a) {
            int32_t f = C.face[a];
            for (int64_t b = C.off[f]; b < C.off[f + 1]; ++b)
                if (__builtin_add_overflow(acc[C.face[b]], (__int128)C.coef[a] * C.coef[b], &acc[C.face[b]])) { over = true; break; }
        }
        if (over) {                                // huge coefficients: redo this cell in BigInt
            std::map<int32_t, BigInt> big;
            for (int64_t a = C.off[c]; a < C.off[c + 1]; ++a) {
                int32_t f = C.face[a];
                for (int64_t b = C.off[f]; b < C.off[f + 1]; ++b) {
                    auto it = big.find(C.face[b]); BigInt cur = it == big.end() ? zlin::zbig((int64_t)0) : it->second;
                    big[C.face[b]] = zlin::zadd(cur, zlin::zmul(zlin::zbig(C.coef[a]), zlin::zbig(C.coef[b])));
                }
            }
            for (auto& kv : big) if (!zlin::zzero(kv.second)) return "boundary of boundary is nonzero at cell " + std::to_string(c) + " (on cell " + std::to_string(kv.first) + ")";
            continue;
        }
        for (auto& kv : acc) if (kv.second != 0) return "boundary of boundary is nonzero at cell " + std::to_string(c) + " (on cell " + std::to_string(kv.first) + ")";
    }
    return "";
}

struct ReducedComplex {
    FilteredComplex C;                         // the residual complex
    std::vector<int32_t> orig;                 // residual cell -> original cell
    std::vector<std::pair<int, double>> zero;  // the cancelled pairs, as (degree, birth): zero-length bars
    size_t pairs = 0, skipped_overflow = 0, pairs_phase1 = 0;   // phase 1: collapses and coreductions (no fill-in)
};

// cancel unit pairs with equal births until none is left (greedy; a face with the smallest coboundary is
// preferred, to limit fill-in). Coefficient growth is checked: a pair whose update would overflow int64 is
// left alone.
inline ReducedComplex algebraic_morse(const FilteredComplex& C) {
    size_t n = C.size();
    std::vector<std::vector<std::pair<int32_t, int64_t>>> bd(n);
    std::vector<std::vector<int32_t>> cob(n);
    for (size_t c = 0; c < n; ++c) {
        std::map<int32_t, int64_t> m;
        for (int64_t a = C.off[c]; a < C.off[c + 1]; ++a) m[C.face[a]] += C.coef[a];
        for (auto& kv : m) if (kv.second) { bd[c].push_back(kv); cob[kv.first].push_back((int32_t)c); }
    }
    std::vector<char> alive(n, 1);
    std::vector<int32_t> ncob(n, 0);           // live cofaces whose boundary contains the cell
    for (size_t c = 0; c < n; ++c) for (auto& e : bd[c]) ++ncob[e.first];
    ReducedComplex R;
    // process cells by decreasing dimension: cancelling top pairs first keeps the lower boundaries sparse
    std::vector<int32_t> ord(n); for (size_t i = 0; i < n; ++i) ord[i] = (int32_t)i;
    std::stable_sort(ord.begin(), ord.end(), [&](int32_t a, int32_t b) { return C.dim[a] > C.dim[b]; });
    auto live_cob = [&](int32_t s) {           // compact the coboundary list of s, return its size
        auto& L = cob[s]; size_t w = 0;
        std::sort(L.begin(), L.end()); L.erase(std::unique(L.begin(), L.end()), L.end());
        for (int32_t r : L) {
            if (!alive[r]) continue;
            bool has = false; for (auto& e : bd[r]) if (e.first == s) { has = true; break; }
            if (has) L[w++] = r;
        }
        L.resize(w); return w;
    };
    std::vector<std::pair<int32_t, int64_t>> tmp;
    // cancel (sigma, tau), [tau : sigma] = ts = +-1: rho -= [rho : sigma] ts tau for the other cofaces rho of
    // sigma (applied only if every coefficient fits in int64), then tau leaves the boundaries of its cofaces
    auto cancel = [&](int32_t sigma, int32_t tau, int64_t ts) -> bool {
        std::vector<std::pair<int32_t, std::vector<std::pair<int32_t, int64_t>>>> upd;
        std::vector<std::pair<int32_t, int32_t>> fresh;   // (face, rho): faces new to rho's boundary
        bool ok = true;
        for (int32_t rho : cob[sigma]) {
            if (rho == tau || !alive[rho]) continue;
            int64_t c = 0; for (auto& e : bd[rho]) if (e.first == sigma) { c = e.second; break; }
            if (!c) continue;
            int64_t f; if (__builtin_mul_overflow(c, ts, &f)) { ok = false; break; }   // rho -= f * tau
            tmp.clear(); size_t a = 0, b = 0; const auto& X = bd[rho]; const auto& Y = bd[tau];
            while (ok && (a < X.size() || b < Y.size())) {
                if (b == Y.size() || (a < X.size() && X[a].first < Y[b].first)) { tmp.push_back(X[a++]); continue; }
                int64_t m; if (__builtin_mul_overflow(f, Y[b].second, &m)) { ok = false; break; }
                if (a == X.size() || Y[b].first < X[a].first) { if (m) { tmp.push_back({Y[b].first, -m}); fresh.push_back({Y[b].first, rho}); } ++b; continue; }
                int64_t v; if (__builtin_sub_overflow(X[a].second, m, &v)) { ok = false; break; }
                if (v) tmp.push_back({X[a].first, v});
                ++a; ++b;
            }
            if (!ok) break;
            upd.push_back({rho, tmp});
        }
        if (!ok) { ++R.skipped_overflow; return false; }
        for (auto& u : upd) {
            for (auto& e : bd[u.first]) --ncob[e.first];
            for (auto& e : u.second) ++ncob[e.first];
            bd[u.first].swap(u.second);
        }
        for (auto& fr : fresh) cob[fr.first].push_back(fr.second);
        for (int32_t pi : cob[tau]) {
            if (!alive[pi]) continue;
            auto& B = bd[pi];
            B.erase(std::remove_if(B.begin(), B.end(), [&](const std::pair<int32_t, int64_t>& e) { return e.first == tau; }), B.end());
        }
        for (auto& e : bd[tau]) --ncob[e.first];
        for (auto& e : bd[sigma]) --ncob[e.first];
        alive[sigma] = alive[tau] = 0;
        R.zero.push_back({C.dim[sigma], C.birth[sigma]});
        ++R.pairs;
        return true;
    };
    auto unit_pair = [&](int32_t sigma, int32_t tau, int64_t e) { return (e == 1 || e == -1) && C.birth[sigma] == C.birth[tau]; };
    // phase 1, no fill-in: collapses (sigma has the single live coface tau) and coreductions (tau has the
    // single face sigma), driven by a work list
    std::vector<int32_t> work(ord.rbegin(), ord.rend());
    std::vector<char> queued(n, 1);
    auto push = [&](int32_t c) { if (alive[c] && !queued[c]) { queued[c] = 1; work.push_back(c); } };
    auto phase1 = [&](bool collapses) {
        while (!work.empty()) {
            int32_t c = work.back(); work.pop_back(); queued[c] = 0;
            if (!alive[c]) continue;
            int32_t sigma = -1, tau = -1; int64_t e = 0;
            if (bd[c].size() == 1 && unit_pair(bd[c][0].first, c, bd[c][0].second)) { sigma = bd[c][0].first; tau = c; e = bd[c][0].second; }
            if (sigma < 0 && collapses && ncob[c] == 1) {
                live_cob(c);
                int32_t t = cob[c][0];
                for (auto& x : bd[t]) if (x.first == c) e = x.second;
                if (unit_pair(c, t, e)) { sigma = c; tau = t; }
            }
            if (sigma < 0) continue;
            live_cob(sigma); live_cob(tau);
            std::vector<int32_t> nb;
            for (auto& x : bd[tau]) nb.push_back(x.first);
            for (auto& x : bd[sigma]) nb.push_back(x.first);
            for (int32_t r : cob[tau]) nb.push_back(r);
            for (int32_t r : cob[sigma]) nb.push_back(r);
            if (cancel(sigma, tau, e)) { ++R.pairs_phase1; for (int32_t r : nb) push(r); }
        }
    };
    auto requeue = [&]() { for (int32_t c : ord) push(c); };
    for (int pass = 0; pass < 4; ++pass) {
        phase1(false); requeue(); phase1(true);
        size_t before = R.pairs;
        // phase 2: Gaussian elimination of the remaining unit pairs, least fill-in first (Markowitz cost
        // (|cob sigma| - 1)(|bd tau| - 1), from a lazily updated heap); every elimination may unlock phase 1
        auto best_face = [&](int32_t tau, int64_t& ts, uint64_t& cost) -> int32_t {
            int32_t sigma = -1; cost = UINT64_MAX;
            if (!alive[tau]) return -1;
            uint64_t bt = bd[tau].empty() ? 0 : bd[tau].size() - 1;
            for (auto& e : bd[tau]) {
                if (!alive[e.first] || !unit_pair(e.first, tau, e.second)) continue;
                uint64_t c = (uint64_t)std::max(0, ncob[e.first] - 1) * bt;
                if (c < cost) { cost = c; sigma = e.first; ts = e.second; }
            }
            return sigma;
        };
        typedef std::pair<uint64_t, int32_t> HE;
        std::priority_queue<HE, std::vector<HE>, std::greater<HE>> heap;
        for (int32_t tau : ord) { int64_t ts; uint64_t c; if (best_face(tau, ts, c) >= 0) heap.push({c, tau}); }
        while (!heap.empty()) {
            HE top = heap.top(); heap.pop();
            int64_t ts; uint64_t c; int32_t sigma = best_face(top.second, ts, c);
            if (sigma < 0) continue;
            if (c > top.first && !heap.empty() && c > heap.top().first) { heap.push({c, top.second}); continue; }
            int32_t tau = top.second;
            live_cob(sigma); live_cob(tau);
            std::vector<int32_t> nb;
            for (auto& x : bd[tau]) nb.push_back(x.first);
            for (int32_t r : cob[tau]) nb.push_back(r);
            for (int32_t r : cob[sigma]) nb.push_back(r);
            if (!cancel(sigma, tau, ts)) continue;
            for (int32_t r : nb) push(r);
            phase1(true);
            for (int32_t r : nb) { int64_t t2; uint64_t c2; if (best_face(r, t2, c2) >= 0) heap.push({c2, r}); }
        }
        if (R.pairs == before) break;
    }
    // the residual complex, cells in their original order
    std::vector<int32_t> nid(n, -1);
    std::vector<std::pair<int32_t, int64_t>> b2;
    for (size_t c = 0; c < n; ++c) if (alive[c]) { nid[c] = (int32_t)R.orig.size(); R.orig.push_back((int32_t)c); }
    for (int32_t c : R.orig) {
        b2.clear();
        for (auto& e : bd[c]) { if (!alive[e.first]) throw std::logic_error("algebraic_morse: dead face in a residual boundary"); b2.push_back({nid[e.first], e.second}); }
        R.C.add_cell(C.dim[c], C.birth[c], b2);
    }
    return R;
}

// ---- integral outputs on a (residual) complex -------------------------------------------------------
struct IntegralHomology {                       // per degree: H_k = Z^free + sum Z/torsion
    std::vector<size_t> free;
    std::vector<std::vector<BigInt>> torsion;
};

namespace detail {
// matrix of d_{k+1} restricted to rows (k-cells) and columns ((k+1)-cells), given as cell lists
inline zlin::Mat<int64_t> boundary_block(const FilteredComplex& C, const std::vector<int32_t>& rows, const std::vector<int32_t>& cols) {
    std::unordered_map<int32_t, size_t> rix; for (size_t i = 0; i < rows.size(); ++i) rix[rows[i]] = i;
    zlin::Mat<int64_t> M(rows.size(), std::vector<int64_t>(cols.size(), 0));
    for (size_t j = 0; j < cols.size(); ++j)
        for (int64_t a = C.off[cols[j]]; a < C.off[cols[j] + 1]; ++a) {
            auto it = rix.find(C.face[a]);
            if (it == rix.end()) throw std::logic_error("boundary_block: face outside the row set");
            int64_t v; if (__builtin_add_overflow(M[it->second][j], C.coef[a], &v)) throw std::overflow_error("boundary_block: coefficient overflow");
            M[it->second][j] = v;
        }
    return M;
}
inline void check_dense(size_t r, size_t c, size_t cap) {
    if ((double)r * (double)c > (double)cap) throw std::runtime_error("integral output: the residual complex is too large for dense Smith forms (" + std::to_string(r) + " x " + std::to_string(c) + ")");
}
}

// integral homology of the whole complex (all cells), degrees 0..maxdim
inline IntegralHomology integral_homology(const FilteredComplex& C, int maxdim, size_t cap = 40000000) {
    int D = C.max_dim();
    std::vector<std::vector<int32_t>> cells(std::max(D, maxdim) + 2);
    for (size_t c = 0; c < C.size(); ++c) cells[C.dim[c]].push_back((int32_t)c);
    std::vector<size_t> rank(cells.size() + 1, 0);
    std::vector<std::vector<BigInt>> tors(cells.size() + 1);
    for (int k = 1; k < (int)cells.size(); ++k) {           // d_k: k-cells -> (k-1)-cells
        if (k > maxdim + 1 || cells[k].empty() || cells[k - 1].empty()) continue;
        detail::check_dense(cells[k - 1].size(), cells[k].size(), cap);
        auto d = zlin::invariant_factors(detail::boundary_block(C, cells[k - 1], cells[k]));
        rank[k] = d.size();
        for (auto& x : d) if (!(x == BigInt(1))) tors[k].push_back(x);
    }
    IntegralHomology H;
    for (int k = 0; k <= maxdim; ++k) {
        size_t nk = k < (int)cells.size() ? cells[k].size() : 0;
        size_t rk = rank[k], rk1 = k + 1 < (int)rank.size() ? rank[k + 1] : 0;
        H.free.push_back(nk - rk - rk1);
        H.torsion.push_back(k + 1 < (int)tors.size() ? tors[k + 1] : std::vector<BigInt>());
    }
    return H;
}

// the image of H_k(K_s; Z) -> H_k(K_t; Z), s <= t, for k = 0..maxdim: (Z_k(K_s) + B_k(K_t)) / B_k(K_t)
inline IntegralHomology integral_image(const FilteredComplex& C, double s, double t, int maxdim, size_t cap = 40000000) {
    if (!(s <= t)) throw std::invalid_argument("integral image: need s <= t");
    int D = C.max_dim();
    auto cells_at = [&](int k, double x) { std::vector<int32_t> v; for (size_t c = 0; c < C.size(); ++c) if (C.dim[c] == k && C.birth[c] <= x) v.push_back((int32_t)c); return v; };
    IntegralHomology H;
    for (int k = 0; k <= maxdim; ++k) {
        std::vector<int32_t> kt = k <= D ? cells_at(k, t) : std::vector<int32_t>();
        std::vector<int32_t> ks = k <= D ? cells_at(k, s) : std::vector<int32_t>();
        std::vector<int32_t> km = k >= 1 && k - 1 <= D ? cells_at(k - 1, s) : std::vector<int32_t>();
        std::vector<int32_t> kp = k + 1 <= D ? cells_at(k + 1, t) : std::vector<int32_t>();
        size_t n = kt.size();
        detail::check_dense(km.size() + 1, ks.size() + 1, cap); detail::check_dense(n + 1, kp.size() + ks.size() + 1, cap);
        // Z_k(K_s): kernel of d_k on the k-cells of K_s, written in the coordinates of the k-cells of K_t
        std::vector<std::vector<BigInt>> Z;
        if (k == 0 || km.empty()) { for (size_t j = 0; j < ks.size(); ++j) { std::vector<BigInt> v(ks.size()); v[j] = BigInt(1); Z.push_back(v); } }
        else Z = zlin::kernel_basis(detail::boundary_block(C, km, ks), ks.size());
        std::unordered_map<int32_t, size_t> tix; for (size_t i = 0; i < n; ++i) tix[kt[i]] = i;
        std::vector<std::vector<BigInt>> S, Bg;
        for (auto& z : Z) { std::vector<BigInt> v(n); for (size_t j = 0; j < ks.size(); ++j) v[tix.at(ks[j])] = z[j]; S.push_back(v); }
        // B_k(K_t): boundaries of the (k+1)-cells of K_t
        if (!kp.empty()) {
            auto M = detail::boundary_block(C, kt, kp);
            for (size_t j = 0; j < kp.size(); ++j) { std::vector<BigInt> v(n); bool nz = false; for (size_t i = 0; i < n; ++i) if (M[i][j]) { v[i] = BigInt(M[i][j]); nz = true; } if (nz) Bg.push_back(v); }
        }
        for (auto& b : Bg) S.push_back(b);
        zlin::GroupStructure G = zlin::quotient(S, Bg, n);
        H.free.push_back(G.free_rank); H.torsion.push_back(G.torsion);
    }
    return H;
}

} // namespace perseus

// ---- a certify() backend on the reduced complex -------------------------------------------------------
#include "explicit_engine.hpp"
namespace perseus {
// Runs the explicit engine on the residual complex of an algebraic Morse reduction; with zero-length bars
// requested, the cancelled pairs are added back as [b, b) in their degree.
struct ReducedBackend {
    typedef ExplicitBackend::Clear Clear;
    const ReducedComplex& R;
    const Prepared& P;
    ExplicitBackend inner;
    bool keep_zero_length = false;
    int dmax;
    ReducedBackend(const ReducedComplex& r, const Prepared& p, int maxdim, int direction = ExplicitBackend::COHOMOLOGY)
        : R(r), P(p), inner(p, std::max(0, std::min(maxdim, p.D)), direction), dmax(maxdim) {}
    int maxdim() const { return dmax; }
    template <class Rg>
    BRun<Clear> run(const Rg& ring, int kmin, int kmax, const Clear* in) const {
        kmax = std::min(kmax, dmax);
        BRun<Clear> out;
        if (P.D >= 0 && kmin <= std::min(kmax, P.D)) {
            ExplicitBackend e = inner; e.keep_zero_length = keep_zero_length;
            out = e.run(ring, kmin, std::min(kmax, P.D), in);
        } else { out.kmin = kmin; for (int l = 0; l < Rg::lanes; ++l) out.moduli.push_back(ring.modulus(l)); }
        int nd = std::max(0, kmax + 1);
        if ((int)out.bars.size() < nd) { out.bars.resize(nd); out.clear_out.resize(nd); }
        if (keep_zero_length)
            for (auto& z : R.zero) if (z.first >= kmin && z.first <= kmax) out.bars[z.first].push_back({z.second, z.second});
        for (int k = kmin; k < nd; ++k) std::sort(out.bars[k].begin(), out.bars[k].end());
        return out;
    }
};
}
