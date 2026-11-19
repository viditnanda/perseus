// certify.hpp -- persistence over Q with a certificate of every prime p for which the F_p barcode
// differs ("bad primes"), for any backend. See docs/CERTIFICATE.md for the proofs.
//
// Lemma 2: if p divides no numerator or denominator of a pivot value of the reduction over Q, that
// reduction taken mod p is a valid F_p reduction with the same pivots, hence the same barcode. So the bad
// primes are among the prime divisors of the pivot values; each candidate is settled by a rerun mod p.
// Pivot values do not depend on the reduction algorithm: the pivot of the pair (i, j) is
// +-det D[I, J] / det D[I', J'], with I, J the pairs nested in (i, j) (I', J' without (i, j) itself).
//
// Two ways to obtain the Q reduction:
//
//  * proved: exact rational arithmetic (int64 with overflow checks, BigInt for any degree that overflows).
//    The Q barcode and the pivot values are exact; the certificate is a proof, up to the factoring of
//    the pivot values (any factor that could not be split is reported as unresolved).
//
//  * fast (default): the reduction is carried out modulo K random primes from [2^61, 2^62) at once. Each
//    lane computes the Q reduction mod its prime until the first Q-pivot divisible by that prime; if some
//    but not all lanes vanish there, the run is detected as unlucky and repeated with fresh primes. It can
//    go wrong undetectably only if a nonzero integer fixed by the input (a pivot numerator, or numerator -+
//    denominator for a pivot that would be taken for a unit, or the discrepancy of a wrong reconstruction)
//    is divisible by every lane prime. These integers are bounded via Hadamard, so for every input
//        Pr[wrong barcode or wrong candidate list] <= sum_k 3 P_k (m_k / N)^K + (reconstruction terms),
//    with P_k pairs in degree k, m_k = floor((h_k + 2) / 61) the number of prime factors >= 2^61 an integer
//    below 2^(h_k + 2) can have, and N >= 5.39e16 primes in the window. The bound is reported.
//
// A Backend provides
//   typedef ... Clear;                     clearing data passed from degree k to degree k+1
//   int maxdim() const;
//   template <class R> BRun<Clear> run(const R& ring, int kmin, int kmax, const Clear* in) const;
#pragma once
#include "ring.hpp"
#include "numtheory.hpp"
#include "types.hpp"
#include <chrono>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <thread>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <limits>
#include <algorithm>
#include <random>
#include <cmath>
#include <cstdlib>
#include <cstdio>

namespace perseus {

// run over the field F_p with the cheapest ring
template <class B>
BRun<typename B::Clear> run_field(const B& be, uint64_t p, int kmin, int kmax, const typename B::Clear* in) {
    if (p == 2) return be.run(GF2(), kmin, kmax, in);
    if (p < (1ULL << 31)) return be.run(Mont<uint32_t, 1>({p}), kmin, kmax, in);
    return be.run(Mont<uint64_t, 1>({p}), kmin, kmax, in);
}
// the Q run modulo the given primes in [2^61, 2^62)
template <class B>
BRun<typename B::Clear> run_lanes(const B& be, const std::vector<uint64_t>& mods, int kmin, int kmax, const typename B::Clear* in) {
    switch (mods.size()) {
        case 2: return be.run(Mont<uint64_t, 2>(mods), kmin, kmax, in);
        case 3: return be.run(Mont<uint64_t, 3>(mods), kmin, kmax, in);
        case 4: return be.run(Mont<uint64_t, 4>(mods), kmin, kmax, in);
        case 8: return be.run(Mont<uint64_t, 8>(mods), kmin, kmax, in);
        case 16: return be.run(Mont<uint64_t, 16>(mods), kmin, kmax, in);
        case 32: return be.run(Mont<uint64_t, 32>(mods), kmin, kmax, in);
        default: throw std::runtime_error("run_lanes: unsupported number of lanes");
    }
}

struct PivotValue {
    int k; uint64_t key; double birth, death;
    bool ok = false; BigInt num, den;
};
struct BadPrime {
    uint64_t p; int k0;
    BigInt big; bool is_big = false;      // p >= 2^63 (then a Miller-Rabin probable prime if >= 2^64)
    bool group = false;                   // big is a number whose prime factors are not known: the entry says
                                          // that EVERY prime factor of big behaves this way
    std::string str() const { return is_big ? big.str() : std::to_string(p); }
    Barcode only_p, only_q;               // bars present over F_p but not over Q, and vice versa
};
struct Certificate {
    std::string mode = "fast";            // "fast", "proved", or "unimodular" (proved by total unimodularity)
    // what is guaranteed about the Q barcode and the pivot values
    bool barcode_ok = true;               // false: the run could not be made reliable (see note)
    std::string note;
    double log2_error = -std::numeric_limits<double>::infinity();   // fast mode: log2 of the error bound
    std::string seed;                     // fast mode: the 256-bit seed, as 64 hex digits
    bool seed_user = false;               // the seed came from --seed: the run reproduces, it does not certify
    int lanes = 0, retries = 0;
    std::vector<uint64_t> moduli;
    std::string arithmetic;               // proved mode: "int64", or which degrees needed bigint
    // the pivot values and what follows from them
    std::vector<PivotValue> values;       // non-unit pivots
    int unreconstructed = 0;
    std::map<uint64_t, int> candidates;   // prime < 2^63 -> first degree where it divides a pivot
    std::vector<std::pair<BigInt, int>> big_candidates;   // primes >= 2^63 (probable primes if >= 2^64)
    std::vector<BigInt> unresolved;       // composite factors we could not split (only with classification off)
    std::vector<BigInt> groups;           // big candidates not known to be prime
    bool classified = false;              // every candidate < 2^63 was settled by a rerun
    std::vector<BadPrime> bad;
    // A fast run certifies only when its seed was drawn afresh: the error bound is a probability over the
    // seed, so it means nothing for a seed chosen before (or together with) the input.
    bool certified() const { return barcode_ok && !seed_user; }
    bool candidates_complete() const { return certified() && unreconstructed == 0 && unresolved.empty(); }
    bool complete() const { return candidates_complete() && classified; }
    double t_q = 0, t_recon = 0, t_classify = 0;
};

inline Barcode bar_diff(const Barcode& a, const Barcode& b, int kmin) {   // multiset a \ b, degrees >= kmin
    Barcode out(a.size());
    for (size_t k = kmin; k < a.size(); ++k) {
        std::vector<Bar> y = k < b.size() ? b[k] : std::vector<Bar>();
        std::set_difference(a[k].begin(), a[k].end(), y.begin(), y.end(), std::back_inserter(out[k]));
    }
    return out;
}
inline bool bars_empty(const Barcode& b) { for (auto& v : b) if (!v.empty()) return false; return true; }

struct CertifyOptions {
    int maxdim = 1 << 20;
    int threads = 1;
    bool classify = true;                 // rerun mod each candidate prime
    std::string mode = "proved";          // "proved" (exact; the default), "auto" (the same), or "fast" (random primes)
    std::string seed;                     // fast mode: hex seed to reproduce a run; empty = fresh from std::random_device
    double target_log2 = -64;             // fast mode: aim for an error bound below 2^target
    int max_lanes = 32;
};

// log2 of sum_k c * P_k * (m_k / N)^K, with m_k = floor((h_k + extra) / 61)
inline double log2_error_terms(const std::vector<double>& h, const std::vector<double>& P, int K, double extra = 2, double c = 3) {
    const double N = primes_in_window(62);
    double best = -std::numeric_limits<double>::infinity(), sum = 0, mx = -1e300;
    std::vector<double> t;
    for (size_t k = 0; k < h.size(); ++k) {
        if (P[k] <= 0) continue;
        double m = std::floor((h[k] + extra) / 61.0);
        if (m <= 0) continue;                                     // integers below 2^61 have no such factor
        double l = std::log2(c * P[k]) + K * (std::log2(m) - std::log2(N));
        t.push_back(l); mx = std::max(mx, l);
    }
    if (t.empty()) return best;
    for (double l : t) sum += std::exp2(l - mx);
    return mx + std::log2(sum);
}
inline double log2_add(double a, double b) {
    if (std::isinf(a) && a < 0) return b;
    if (std::isinf(b) && b < 0) return a;
    double m = std::max(a, b); return m + std::log2(std::exp2(a - m) + std::exp2(b - m));
}

namespace detail {
template <class Clear>
void merge_run(BRun<Clear>& into, BRun<Clear>&& r, int k) {
    if ((int)into.bars.size() <= k) { into.bars.resize(k + 1); into.clear_out.resize(k + 1); }
    if ((int)r.bars.size() > k) into.bars[k] = std::move(r.bars[k]);
    if ((int)r.clear_out.size() > k) into.clear_out[k] = std::move(r.clear_out[k]);
    for (auto& pv : r.nonunit) if (pv.k == k) into.nonunit.push_back(std::move(pv));
    into.additions += r.additions; into.columns += r.columns;
    for (size_t j = 0; j < r.hbits.size(); ++j) into.note_bound((int)j, r.hbits[j], r.pairs[j]);
}
}

// proved mode: the Q run in exact arithmetic, degree by degree, falling back to BigInt where int64 overflows
template <class B>
BRun<typename B::Clear> run_exact(const B& be, int kmax, std::string& arithmetic) {
    typedef typename B::Clear Clear;
    BRun<Clear> out; out.kmin = 0;
    auto T0 = std::chrono::steady_clock::now();
    { Q64 R; BRun<Clear> r = be.run(R, 0, kmax, nullptr); if (!R.overflow) { arithmetic = "int64"; return r; } }
    if (getenv("PERSEUS_DEBUG")) fprintf(stderr, "  [debug: int64 overflow after %.3f s]\n", std::chrono::duration<double>(std::chrono::steady_clock::now() - T0).count());
    std::vector<int> big;
    for (int k = 0; k <= kmax; ++k) {
        const Clear* in = k > 0 ? &out.clear_out[k - 1] : nullptr;
        Q64 R; BRun<Clear> r = be.run(R, k, k, in);
        if (R.overflow) { r = be.run(QHyb(), k, k, in); big.push_back(k); }   // exact: int64 pairs and big rationals
        detail::merge_run(out, std::move(r), k);
    }
    arithmetic = "int64, big rationals in degree"; for (size_t i = 0; i < big.size(); ++i) arithmetic += (i ? ", " : " ") + std::to_string(big[i]);
    if (big.empty()) arithmetic = "int64";
    return out;
}

// run modulo N (odd, possibly composite); returns false and sets `split` to a proper factor of N if a
// pivot shares one with N
template <class B>
bool run_mod_big(const B& be, const BigInt& N, int k0, int kmax, const typename B::Clear* in, BRun<typename B::Clear>& r, BigInt& split) {
    size_t b = N.bits();
    auto go = [&](auto ring) { r = be.run(ring, k0, kmax, in); if (ring.split_found) { split = ring.split; return false; } return true; };
    if (b <= 126) return go(MontN<2>(N));
    if (b <= 190) return go(MontN<3>(N));
    if (b <= 254) return go(MontN<4>(N));
    if (b <= 382) return go(MontN<6>(N));
    if (b <= 510) return go(MontN<8>(N));
    if (b <= 1022) return go(MontN<16>(N));
    return go(MontN<32>(N));
}

// backends whose reductions are of totally unimodular matrices declare unimodular() == true
template <class B> auto backend_unimodular(const B& b, int) -> decltype(b.unimodular()) { return b.unimodular(); }
template <class B> bool backend_unimodular(const B&, long) { return false; }

template <class B>
std::pair<BRun<typename B::Clear>, Certificate> certify(const B& be, const CertifyOptions& opt = CertifyOptions()) {
    typedef typename B::Clear Clear;
    int kmax = std::min(opt.maxdim, be.maxdim());
    Certificate cert; cert.mode = opt.mode;
    if (backend_unimodular(be, 0)) {
        // every pivot of a totally unimodular matrix is +-1 (Lemma 2: pivots are ratios of minors), so the
        // barcode is the same over Q and every F_p: no arithmetic to check, no candidates, no bad primes
        Q64 R; BRun<Clear> q = be.run(R, 0, kmax, nullptr);
        cert.mode = "unimodular"; cert.arithmetic = "none needed: totally unimodular boundary matrices";
        cert.classified = true;
        return {q, cert};
    }
    auto now = [] { return std::chrono::steady_clock::now(); };
    auto t0 = now();
    BRun<Clear> q;
    std::map<std::pair<int, uint64_t>, size_t> where;
    // The default is a proof: exact rationals, int64 where they fit and big rationals in the degrees where
    // they do not (run_exact). Random primes only when asked for ("fast"); "auto" is a synonym of "proved".
    if (opt.mode != "fast") {
        q = run_exact(be, kmax, cert.arithmetic); cert.mode = "proved";
        for (auto& pv : q.nonunit) {
            PivotValue v; v.k = pv.k; v.key = pv.key; v.birth = pv.birth; v.death = pv.death;
            v.ok = true; v.num = pv.num; v.den = pv.den;
            if (v.den.neg) { v.num = -v.num; v.den = -v.den; }
            cert.values.push_back(v);
        }
    } else {
        // 1. the Q run modulo K random 62-bit primes, retried with fresh primes when unlucky
        std::array<uint64_t, 4> sd = opt.seed.empty() ? fresh_seed() : parse_seed(opt.seed);
        cert.seed = seed_hex(sd); cert.seed_user = !opt.seed.empty();
        ChaCha20 rng(sd);
        std::set<uint64_t> used;
        auto fresh = [&](int K) { auto m = random_primes(K, 62, rng, used); for (auto x : m) used.insert(x); return m; };
        int K = 2;
        auto attempt_run = [&](int Kr, int kmin, int kmax2, const Clear* in, BRun<Clear>& r) {
            for (int t = 0; t < 7; ++t) {
                r = run_lanes(be, fresh(Kr), kmin, kmax2, in);
                if (!r.unlucky) return true;
                ++cert.retries;
            }
            return false;
        };
        bool ok = attempt_run(K, 0, kmax, nullptr, q);
        if (getenv("PERSEUS_DEBUG")) fprintf(stderr, "  [debug: first run, %d lanes, %.3f s, %zu non-unit]\n", K, std::chrono::duration<double>(now() - t0).count(), q.nonunit.size());
        double err = log2_error_terms(q.hbits, q.pairs, K);
        if (ok && err > opt.target_log2) {               // more lanes for the requested error bound
            int K2 = K;
            for (int c : {3, 4, 8, 16, 32}) { K2 = c; if (log2_error_terms(q.hbits, q.pairs, c) <= opt.target_log2 || c >= opt.max_lanes) break; }
            if (K2 != K) { K = K2; ok = attempt_run(K, 0, kmax, nullptr, q); err = log2_error_terms(q.hbits, q.pairs, K); }
        }
        if (!ok) { cert.barcode_ok = false; cert.note = "every set of trial primes divided a pivot"; }
        cert.lanes = K; cert.moduli = q.moduli; cert.log2_error = err;
        // 2. reconstruct non-unit pivot values: Wang reconstruction from all lanes but two, checked on those two.
        //    Degrees holding non-unit pivots are rerun with more lanes (fresh primes) as needed; the rerun
        //    must reproduce the pairing of q in that degree.
        for (auto& pv : q.nonunit) {
            PivotValue v; v.k = pv.k; v.key = pv.key; v.birth = pv.birth; v.death = pv.death;
            where[{pv.k, pv.key}] = cert.values.size(); cert.values.push_back(v);
        }
        auto reconstruct = [&](const std::vector<Pivot>& piv, const std::vector<uint64_t>& mods, double& errsum, const std::vector<double>& hb) {
            int L = (int)mods.size(); if (L < 3) return;
            for (auto& pv : piv) {
                auto it = where.find({pv.k, pv.key}); if (it == where.end()) continue;
                PivotValue& v = cert.values[it->second]; if (v.ok) continue;
                BigInt x, M, n, d;
                crt(std::vector<uint64_t>(pv.res.begin(), pv.res.begin() + (L - 2)), std::vector<uint64_t>(mods.begin(), mods.begin() + (L - 2)), x, M);
                if (!ratrecon(x, M, n, d)) continue;
                bool agree = true;
                for (int l = L - 2; l < L && agree; ++l) {
                    uint64_t p = mods[l];
                    uint64_t nm = n.mod_small(p), dm = d.mod_small(p);
                    if (n.neg && nm) nm = p - nm;
                    agree = mulmod64(pv.res[l], dm, p) == nm;
                }
                if (!agree) continue;
                v.ok = true; v.num = n; v.den = d;
                double h = pv.k < (int)hb.size() ? hb[pv.k] : 0;
                double m = std::floor((h + 2 + (double)std::max(n.bits(), d.bits())) / 61.0);
                if (m > 0) errsum = log2_add(errsum, 2 * (std::log2(m) - std::log2(primes_in_window(62))));
            }
        };
        double rec_err = -std::numeric_limits<double>::infinity();
        reconstruct(q.nonunit, q.moduli, rec_err, q.hbits);
        int Kr = K;
        while (cert.barcode_ok && Kr < opt.max_lanes) {
            std::set<int> degs; for (auto& v : cert.values) if (!v.ok) degs.insert(v.k);
            if (degs.empty()) break;
            Kr = Kr < 4 ? 4 : 2 * Kr;
            for (int k : degs) {
                BRun<Clear> r;
                if (!attempt_run(Kr, k, k, k > 0 ? &q.clear_out[k - 1] : nullptr, r)) { cert.barcode_ok = false; cert.note = "every set of trial primes divided a pivot"; break; }
                if (r.bars.size() <= (size_t)k || r.bars[k] != q.bars[k]) { cert.barcode_ok = false; cert.note = "independent runs disagree in degree " + std::to_string(k); break; }
                err = log2_add(err, log2_error_terms(r.hbits, r.pairs, Kr));
                reconstruct(r.nonunit, r.moduli, rec_err, r.hbits);
            }
            cert.lanes = std::max(cert.lanes, Kr);
            if (getenv("PERSEUS_DEBUG")) { int nr = 0; for (auto& v : cert.values) nr += !v.ok; fprintf(stderr, "  [debug: after %d lanes, %.3f s, %d unreconstructed]\n", Kr, std::chrono::duration<double>(now() - t0).count(), nr); }
        }
        cert.log2_error = log2_add(err, rec_err);
    }
    for (auto& v : cert.values) if (!v.ok) ++cert.unreconstructed;
    auto t1 = now();
    bool dbg = getenv("PERSEUS_DEBUG") != nullptr;
    if (dbg) { size_t mx = 0; for (auto& v : cert.values) if (v.ok) mx = std::max(mx, std::max(v.num.bits(), v.den.bits()));
               fprintf(stderr, "  [debug: Q run %.3f s, %zu non-unit pivots, largest %zu bits]\n", std::chrono::duration<double>(t1 - t0).count(), cert.values.size(), mx); }
    if (getenv("PERSEUS_DUMP")) { FILE* f = fopen(getenv("PERSEUS_DUMP"), "w"); for (auto& v : cert.values) if (v.ok) fprintf(f, "%d %s %s %.17g %.17g\n", v.k, v.num.str().c_str(), v.den.str().c_str(), v.birth, v.death); fclose(f); }
    // 3. candidates: prime divisors of numerators and denominators. Pivot values share many factors, so
    //    each value is first divided by the primes already found (smallest values first), and composites
    //    that resist rho are split further by gcds with the other values' cofactors.
    std::vector<size_t> ord; for (size_t i = 0; i < cert.values.size(); ++i) if (cert.values[i].ok) ord.push_back(i);
    std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return std::max(cert.values[a].num.bits(), cert.values[a].den.bits()) < std::max(cert.values[b].num.bits(), cert.values[b].den.bits()); });
    std::vector<BigInt> known;                            // primes found so far (any size)
    std::vector<std::pair<BigInt, int>> hard;             // composite cofactors rho could not split, with degree
    // A number may turn up in several degrees; it must always keep the EARLIEST one, since a rerun modulo
    // its primes starts there and would miss changes in any earlier degree.
    auto add_pair = [](std::vector<std::pair<BigInt, int>>& L, const BigInt& b, int k) {
        for (auto& c : L) if (c.first == b) { c.second = std::min(c.second, k); return; }
        L.push_back({b, k});
    };
    auto add_small = [&](uint64_t p, int k) {
        auto it = cert.candidates.find(p); if (it == cert.candidates.end() || it->second > k) cert.candidates[p] = k;
    };
    for (size_t idx : ord) {
        auto& v = cert.values[idx];
        std::map<uint64_t, int> pr; std::vector<BigInt> bp, unres;
        for (BigInt n : {v.num.abs(), v.den}) {
            for (auto& p : known) { if (n.bits() < p.bits()) continue; BigInt qq, rr; divmod(n, p, qq, rr); while (rr.is_zero()) { n = qq; divmod(n, p, qq, rr); } }
            factor_big(n, pr, unres, bp);
        }
        for (auto& kv : pr) known.push_back(BigInt::from_u64(kv.first));
        for (auto& b : bp) known.push_back(b);
        for (auto& u : unres) add_pair(hard, u, v.k);
        // record the primes of this value (including known ones dividing it)
        for (auto& p : known) {
            bool divides = false;
            for (const BigInt* n : {&v.num, &v.den}) { BigInt qq, rr; BigInt a = n->abs(); if (a.bits() < p.bits()) continue; divmod(a, p, qq, rr); if (rr.is_zero()) divides = true; }
            if (!divides) continue;
            if (p.fits_u64() && p.to_u64() < (1ULL << 63)) pr[p.to_u64()] = 1;
            else { bool dup = false; for (auto& b : bp) if (b == p) dup = true; if (!dup) bp.push_back(p); }
        }
        for (auto& kv : pr) {
            if (kv.first >= (1ULL << 63)) { bp.push_back(BigInt::from_u64(kv.first)); continue; }
            auto it = cert.candidates.find(kv.first); if (it == cert.candidates.end() || it->second > v.k) cert.candidates[kv.first] = v.k;
        }
        for (auto& b : bp) add_pair(cert.big_candidates, b, v.k);
    }
    {   // pairwise gcds between hard composites, until nothing splits. Equal composites were merged above
        // (keeping the earliest degree); a composite dividing another splits the larger one.
        bool progress = true;
        while (progress && !hard.empty()) {
            progress = false;
            for (size_t i = 0; i < hard.size() && !progress; ++i)
                for (size_t j = i + 1; j < hard.size() && !progress; ++j) {
                    BigInt g = gcd(hard[i].first, hard[j].first);
                    if (g == BigInt(1) || (g == hard[i].first && g == hard[j].first)) continue;
                    BigInt a = hard[i].first, b = hard[j].first; int ka = hard[i].second, kb = hard[j].second;
                    hard.erase(hard.begin() + j); hard.erase(hard.begin() + i);
                    // g divides both, so it inherits the earlier degree; each cofactor keeps its own
                    for (auto pk : {std::make_pair(g, std::min(ka, kb)), std::make_pair(a / g, ka), std::make_pair(b / g, kb)}) {
                        const BigInt& part = pk.first; int k = pk.second;
                        if (part == BigInt(1)) continue;
                        std::map<uint64_t, int> pr; std::vector<BigInt> bp, unres;
                        factor_big(part, pr, unres, bp);
                        for (auto& kv : pr) { if (kv.first < (1ULL << 63)) add_small(kv.first, k); else bp.push_back(BigInt::from_u64(kv.first)); }
                        for (auto& bb : bp) add_pair(cert.big_candidates, bb, k);
                        for (auto& u : unres) add_pair(hard, u, k);
                    }
                    progress = true;
                }
        }
        for (auto& h : hard) {                         // unsplit composites: classified as groups below
            add_pair(cert.big_candidates, h.first, h.second);
            bool dup = false; for (auto& g : cert.groups) if (g == h.first) dup = true;
            if (!dup) cert.groups.push_back(h.first);
        }
    }
    auto t2 = now();
    if (dbg) fprintf(stderr, "  [debug: factoring %.3f s, %zu + %zu candidates, %zu unresolved]\n", std::chrono::duration<double>(t2 - t1).count(), cert.candidates.size(), cert.big_candidates.size(), cert.unresolved.size());
    // 4. classify candidates by reruns mod p from their first degree
    if (opt.classify && cert.barcode_ok) {
        // big candidates first, sequentially: a run modulo N either classifies every prime factor of N at
        // once, or splits N; split parts that are small primes join the small candidates
        std::vector<std::pair<BigInt, int>> queue = cert.big_candidates, done;
        auto tb0 = now(); int nbig_runs = 0;
        auto is_group = [&](const BigInt& x) { for (auto& g : cert.groups) if (g == x) return true; return false; };
        // the runs are independent: a pool of threads takes jobs from the queue; a split puts its two parts
        // back, so a thread finding the queue empty waits while other runs may still add to it
        std::mutex mu; std::condition_variable cv; int active = 0;
        auto big_worker = [&]() {
            std::unique_lock<std::mutex> lk(mu);
            while (true) {
                cv.wait(lk, [&] { return !queue.empty() || active == 0; });
                if (queue.empty()) return;                      // and nothing running: done
                auto job = queue.back(); queue.pop_back();
                const BigInt N = job.first; const int k0 = job.second;
                if (N.fits_u64() && N.to_u64() < (1ULL << 63)) {    // a small prime or composite: factor it outright
                    if (is_prime64(N.to_u64())) add_small(N.to_u64(), k0);
                    else { std::map<uint64_t, int> f; factor64(N.to_u64(), f); for (auto& kv : f) add_small(kv.first, k0); }
                    continue;
                }
                ++active; ++nbig_runs;
                lk.unlock();
                BRun<Clear> r; BigInt g;
                const bool ok = run_mod_big(be, N, k0, kmax, k0 > 0 ? &q.clear_out[k0 - 1] : nullptr, r, g);
                BadPrime bp;
                if (ok) { bp.only_p = bar_diff(r.bars, q.bars, k0); bp.only_q = bar_diff(q.bars, r.bars, k0); }
                lk.lock();
                --active;
                if (!ok) {
                    BigInt h = N / g;
                    bool grp = is_group(N) || !is_probable_prime(N);
                    for (BigInt part : {g, h}) {
                        add_pair(queue, part, k0);
                        if (grp && !is_probable_prime(part) && !is_group(part)) cert.groups.push_back(part);
                    }
                } else {
                    bp.p = 0; bp.k0 = k0; bp.big = N; bp.is_big = true; bp.group = is_group(N) && !is_probable_prime(N);
                    if (!bars_empty(bp.only_p) || !bars_empty(bp.only_q)) cert.bad.push_back(std::move(bp));
                    done.push_back(job);
                }
                cv.notify_all();
            }
        };
        {
            int nb = std::max(1, std::min<int>(opt.threads, (int)queue.size()));
            std::vector<std::thread> pool; for (int t = 1; t < nb; ++t) pool.emplace_back(big_worker);
            big_worker(); for (auto& th : pool) th.join();
        }
        std::sort(done.begin(), done.end(), [](const std::pair<BigInt, int>& x, const std::pair<BigInt, int>& y) { return x.first < y.first; });
        cert.big_candidates = done;
        auto tb1 = now();
        if (dbg) fprintf(stderr, "  [debug: big candidates %.3f s, %d runs modulo multi-word N, %zu settled]\n", std::chrono::duration<double>(tb1 - tb0).count(), nbig_runs, done.size());
        std::vector<std::pair<uint64_t, int>> cand(cert.candidates.begin(), cert.candidates.end());
        std::vector<double> tcand(cand.size(), 0);
        std::vector<BadPrime> found(cand.size()); std::vector<char> isbad(cand.size(), 0);
        std::atomic<size_t> next(0);
        auto worker = [&]() {
            for (size_t i; (i = next++) < cand.size();) {
                uint64_t p = cand[i].first; int k0 = cand[i].second;
                const Clear* in = k0 > 0 ? &q.clear_out[k0 - 1] : nullptr;
                BadPrime bp; bp.p = p; bp.k0 = k0;
                auto tc0 = std::chrono::steady_clock::now();
                BRun<Clear> r = run_field(be, p, k0, kmax, in);
                tcand[i] = std::chrono::duration<double>(std::chrono::steady_clock::now() - tc0).count();
                bp.only_p = bar_diff(r.bars, q.bars, k0); bp.only_q = bar_diff(q.bars, r.bars, k0);
                if (!bars_empty(bp.only_p) || !bars_empty(bp.only_q)) { found[i] = std::move(bp); isbad[i] = 1; }
            }
        };
        int nt = std::max(1, std::min<int>(opt.threads, (int)cand.size()));
        std::vector<std::thread> pool; for (int t = 1; t < nt; ++t) pool.emplace_back(worker);
        worker(); for (auto& th : pool) th.join();
        for (size_t i = 0; i < cand.size(); ++i) if (isbad[i]) cert.bad.push_back(std::move(found[i]));
        if (dbg) {        // where each small candidate first divides a pivot, against the cost of its rerun
            double tsum = 0; for (double x : tcand) tsum += x;
            fprintf(stderr, "  [debug: small candidates %.3f s wall, %.3f s summed over %zu reruns, %d threads]\n", std::chrono::duration<double>(now() - tb1).count(), tsum, cand.size(), nt);
            if (const char* path = getenv("PERSEUS_CLASSIFY_DUMP")) {
                FILE* f = fopen(path, "w");
                std::map<int, uint64_t> kmaxkey; for (auto& v : cert.values) kmaxkey[v.k] = std::max(kmaxkey[v.k], v.key);
                for (size_t i = 0; i < cand.size(); ++i) {
                    uint64_t p = cand[i].first; int k0 = cand[i].second; uint64_t first = UINT64_MAX, last = 0; int nvals = 0;
                    for (auto& v : cert.values) {
                        if (!v.ok || v.k != k0) continue;
                        bool d = false;
                        for (const BigInt* x : {&v.num, &v.den}) { BigInt a = x->abs(); if (a.mod_small(p) == 0) d = true; }
                        if (d) { first = std::min(first, v.key); last = std::max(last, v.key); ++nvals; }
                    }
                    fprintf(f, "%llu %d %llu %llu %llu %d %.6f %d\n", (unsigned long long)p, k0, (unsigned long long)first, (unsigned long long)last, (unsigned long long)kmaxkey[k0], nvals, tcand[i], (int)isbad[i]);
                }
                fclose(f);
            }
        }
        std::stable_sort(cert.bad.begin(), cert.bad.end(), [](const BadPrime& a, const BadPrime& b) {
            if (a.is_big != b.is_big) return !a.is_big; return a.is_big ? a.big < b.big : a.p < b.p; });
        cert.classified = true;
    }
    auto t3 = now();
    cert.t_q = std::chrono::duration<double>(t1 - t0).count();
    cert.t_recon = std::chrono::duration<double>(t2 - t1).count();
    cert.t_classify = std::chrono::duration<double>(t3 - t2).count();
    return {q, cert};
}

} // namespace perseus
