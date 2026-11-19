// torsion.hpp -- bad primes, their F_p barcodes and the torsion of H_q(K_t; Z) for every t, from Howell
// sweeps of the boundary matrices over Z/m (howell.hpp). Explicit complexes only (a Prepared).
//
// Notation: pairs of degree k are the persistence pairs (row, column) of d_{k+1}, rows the k-cells and columns
// the (k+1)-cells, by local index in filtration order.
//
// * Reference. A prime l in [2^61, 2^63) dividing no pivot value of the Q run gives, by Lemma 2 of
//   docs/CERTIFICATE.md, the same pairs at the level of values as Q; its sweep is the Q side of every
//   comparison below.
// * Classification. Small candidates are packed into squarefree moduli below 2^63, one batch per degree in
//   which they divide a pivot value; one sweep per batch gives the F_p pairs of every prime in it. In the
//   other degrees a prime has the reference pairs (Lemma 2 again). The F_p barcode follows, and the prime is
//   bad iff it differs from the Q barcode.
// * Torsion (Theorem 4's Steps 2-3 with a = -1). The number of cyclic summands of the p-part of
//   tors H_q(K_t; Z) is r_Q(t) - r_{F_p}(t), the ranks of the columns of d_{q+1} born by t. More generally,
//   sum_i min(a_i, k) = k r_Q(t) - log_p |span_t mod p^k| for the p-primary exponents a_i, so sweeps modulo
//   p^k for k = 2, 3, ... give all exponents, stopping at the first k with no exponent >= k.
#pragma once
#include "howell.hpp"
#include "explicit_engine.hpp"
#include "chain.hpp"
#include "certify.hpp"
#include "numtheory.hpp"
#include <thread>
#include <atomic>
#include <mutex>
#include <functional>
#include <set>
#include <map>
#include <cmath>

namespace perseus {
namespace torsion {

using howell::ModU64; using howell::ModBig; using howell::Sweep;
typedef std::vector<std::pair<int32_t, int32_t>> Pairs;

struct PrimeResult {
    Pairs pairs;                                     // modulus exponent 1: the F_p pairs
    std::vector<std::pair<int32_t, int>> dlen;       // (column, change of log_p |span mod p^k|)
    bool split = false;                              // a composite "prime" (group) was split by a pivot ideal
};

template <class M>
static std::vector<typename Sweep<M>::Ent> column(const Prepared& P, int j, int32_t col, const M& R) {
    std::vector<typename Sweep<M>::Ent> x;
    for (int64_t e = P.bd_off[j][col]; e < P.bd_off[j][col + 1]; ++e) {
        auto v = R.from_int(P.bd_coef[j][e]);
        if (!M::is_zero(v)) x.push_back({P.bd_row[j][e], v});
    }
    return x;
}

// sweep d_j over Z/m, m = prod p_i^{k_i}; one result per (p_i, k_i). check_groups: the p_i may be composite
template <class M, class PT>
std::vector<PrimeResult> sweep(const Prepared& P, int j, const M& R, const std::vector<std::pair<PT, int>>& pk, bool check_groups = false) {
    Sweep<M> S(R, P.count(j - 1));
    for (size_t c = 0; c < P.count(j); ++c) S.insert(column(P, j, (int32_t)c, R), (int32_t)c);
    std::vector<PrimeResult> out(pk.size());
    for (auto& ev : S.events)
        for (size_t i = 0; i < pk.size(); ++i) {
            const PT& p = pk[i].first; int k = pk[i].second;
            int vo = M::val(ev.dold, p, k), vn = M::val(ev.dnew, p, k);
            if (check_groups && vn < k) {             // the part of d prime to p^vn must be prime to p
                typename M::T rest = ev.dnew; for (int t = 0; t < vn; ++t) rest = M::div(rest, p);
                if (!M::is_one(M::gcd(rest, p))) out[i].split = true;
            }
            if (k == 1 && vo == 1 && vn == 0) out[i].pairs.push_back({ev.row, ev.col});
            if (vo != vn) out[i].dlen.push_back({ev.col, vo - vn});
        }
    return out;
}

// first-fit decreasing: products below 2^63
inline std::vector<std::vector<std::pair<uint64_t, int>>> pack(std::vector<std::pair<uint64_t, int>> pk) {
    auto lg = [](const std::pair<uint64_t, int>& x) { return x.second * std::log2((double)x.first); };
    std::sort(pk.begin(), pk.end(), [&](auto& a, auto& b) { return lg(a) > lg(b); });
    std::vector<std::vector<std::pair<uint64_t, int>>> bins; std::vector<u128> prod;
    for (auto& x : pk) {
        u128 q = 1; for (int t = 0; t < x.second; ++t) q *= x.first;
        size_t b = 0; for (; b < bins.size(); ++b) if (prod[b] * q < ((u128)1 << 63)) break;
        if (b == bins.size()) { bins.push_back({}); prod.push_back(1); }
        bins[b].push_back(x); prod[b] *= q;
    }
    return bins;
}
inline uint64_t modulus_of(const std::vector<std::pair<uint64_t, int>>& bin) {
    uint64_t m = 1; for (auto& x : bin) for (int t = 0; t < x.second; ++t) m *= x.first; return m;
}

inline void run_jobs(std::vector<std::function<void()>>& jobs, int threads) {
    std::atomic<size_t> next(0);
    auto w = [&]() { for (size_t i; (i = next++) < jobs.size();) jobs[i](); };
    int nt = std::max(1, std::min<int>(threads, (int)jobs.size()));
    std::vector<std::thread> pool; for (int t = 1; t < nt; ++t) pool.emplace_back(w);
    w(); for (auto& t : pool) t.join();
}

// degrees in which a number divides the numerator or the denominator of some pivot value
inline std::set<int> degrees_of(const Certificate& cert, const BigInt& p) {
    std::set<int> D;
    for (auto& v : cert.values) {
        if (!v.ok || D.count(v.k)) continue;
        for (const BigInt* x : {&v.num, &v.den}) {
            BigInt a = x->abs();
            if (p.fits_u64() ? a.mod_small(p.to_u64()) == 0 : (a % p).is_zero()) { D.insert(v.k); break; }
        }
    }
    return D;
}

// a prime in [2^61, 2^63) dividing no pivot value
inline uint64_t reference_prime(const Certificate& cert) {
    for (uint64_t l = (1ULL << 62) + 135;; l += 2) {
        if (!is_prime64(l)) continue;
        bool ok = true;
        for (auto& v : cert.values) { if (!v.ok) continue; if (v.num.abs().mod_small(l) == 0 || v.den.mod_small(l) == 0) { ok = false; break; } }
        if (ok) return l;
    }
}

struct Cache {
    uint64_t ell = 0;
    std::map<int, Pairs> ref;                         // degree -> reference (Q) pairs
    std::map<uint64_t, std::map<int, Pairs>> fp;      // prime -> degree -> F_p pairs, where they differ in general
};

// barcode in degrees k0..kmax (lower degrees empty) from the pairs of degrees k0 - 1 .. kmax
inline Barcode assemble(const Prepared& P, int k0, int kmax, const std::function<const Pairs&(int)>& pairs_of, bool keep0) {
    Barcode B(kmax + 1);
    const double INF = std::numeric_limits<double>::infinity();
    for (int k = k0; k <= kmax && k <= P.D; ++k) {
        std::vector<uint8_t> used(P.count(k), 0);
        if (k < P.D) for (auto& pr : pairs_of(k)) {
            used[pr.first] = 1;
            double b = P.birth[k][pr.first], d = P.birth[k + 1][pr.second];
            if (d > b || keep0) B[k].push_back({b, d});
        }
        if (k >= 1) for (auto& pr : pairs_of(k - 1)) used[pr.second] = 1;
        for (size_t s = 0; s < P.count(k); ++s) if (!used[s]) B[k].push_back({P.birth[k][s], INF});
        std::sort(B[k].begin(), B[k].end());
    }
    return B;
}

inline void ensure_ref(const Prepared& P, Cache& C, const Certificate& cert, const std::set<int>& degs, int threads) {
    if (!C.ell) C.ell = reference_prime(cert);
    std::vector<int> need; for (int k : degs) if (k >= 0 && k < P.D && !C.ref.count(k)) need.push_back(k);
    std::vector<Pairs> res(need.size());
    std::vector<std::function<void()>> jobs;
    for (size_t i = 0; i < need.size(); ++i) jobs.push_back([&, i]() {
        res[i] = sweep(P, need[i] + 1, ModU64(C.ell), std::vector<std::pair<uint64_t, int>>{{C.ell, 1}})[0].pairs; });
    run_jobs(jobs, threads);
    for (size_t i = 0; i < need.size(); ++i) C.ref[need[i]] = std::move(res[i]);
}

// Classify the small candidates (cert.candidates) by sweeps. Returns false if the reference disagrees with the
// Q barcode (then nothing is claimed and the caller reruns). Bad primes get their F_p pairs cached.
inline bool classify(const Prepared& P, const Certificate& cert, const Barcode& qbars, int kmax, int threads, bool keep0,
                     std::vector<BadPrime>& bad, Cache& C, std::string* diag = nullptr) {
    P.build_boundaries();
    kmax = std::min(kmax, P.D);
    std::map<uint64_t, std::set<int>> Dp;
    std::set<int> degs;
    for (auto& kv : cert.candidates) {
        std::set<int> D = degrees_of(cert, BigInt::from_u64(kv.first));
        for (auto it = D.begin(); it != D.end();) { if (*it > kmax || *it >= P.D) it = D.erase(it); else ++it; }
        if (D.empty()) continue;                       // divides no pivot value of a degree we report
        Dp[kv.first] = D;
        int k0 = *D.begin();
        for (int k = std::max(0, k0 - 1); k <= kmax; ++k) degs.insert(k);
    }
    if (Dp.empty()) return true;
    auto t0 = std::chrono::steady_clock::now();
    ensure_ref(P, C, cert, degs, threads);
    {   // the reference must reproduce the Q barcode in every degree it covers
        int lo = *degs.begin(); int k1 = lo == 0 ? 0 : lo + 1;
        Barcode R = assemble(P, k1, kmax, [&](int k) -> const Pairs& { return C.ref.at(k); }, keep0);
        for (int k = k1; k <= kmax; ++k) {
            std::vector<Bar> q = k < (int)qbars.size() ? qbars[k] : std::vector<Bar>();
            if (R[k] != q) { if (diag) *diag = "reference sweep disagrees with the Q barcode in degree " + std::to_string(k); return false; }
        }
    }
    // batches per degree
    struct Job { int k; std::vector<std::pair<uint64_t, int>> bin; std::vector<PrimeResult> res; };
    std::vector<Job> jobs;
    std::map<int, std::vector<std::pair<uint64_t, int>>> per;
    for (auto& kv : Dp) for (int k : kv.second) per[k].push_back({kv.first, 1});
    for (auto& kv : per) for (auto& b : pack(kv.second)) jobs.push_back({kv.first, b, {}});
    std::vector<std::function<void()>> fns;
    for (auto& J : jobs) fns.push_back([&P, &J]() { J.res = sweep(P, J.k + 1, ModU64(modulus_of(J.bin)), J.bin); });
    run_jobs(fns, threads);
    std::map<uint64_t, std::map<int, Pairs>> fp;
    for (auto& J : jobs) for (size_t i = 0; i < J.bin.size(); ++i) fp[J.bin[i].first][J.k] = std::move(J.res[i].pairs);
    for (auto& kv : Dp) {
        uint64_t p = kv.first; int k0 = *kv.second.begin();
        auto& mine = fp[p];
        auto pairs_of = [&](int k) -> const Pairs& { auto it = mine.find(k); return it != mine.end() ? it->second : C.ref.at(k); };
        Barcode B = assemble(P, k0, kmax, pairs_of, keep0);
        BadPrime bp; bp.p = p; bp.k0 = cert.candidates.at(p);
        bp.only_p = bar_diff(B, qbars, k0); bp.only_q = bar_diff(qbars, B, k0);
        if (!bars_empty(bp.only_p) || !bars_empty(bp.only_q)) { bad.push_back(std::move(bp)); C.fp[p] = std::move(mine); }
    }
    if (diag) { char b[160]; snprintf(b, sizeof b, "%zu sweeps for %zu candidates in %.3f s", jobs.size() + C.ref.size(), Dp.size(),
                                      std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()); *diag = b; }
    return true;
}

// ---- torsion of H_q(K_t; Z) for every t --------------------------------------------------------------------
struct Factor { BigInt p; bool group = false; bool split = false; };
struct Step { double value; std::vector<std::pair<int, std::vector<int>>> primary; };   // factor index, exponents descending
struct Trajectory { int degree; std::vector<Step> steps; };
struct TorsionResult {
    std::vector<Factor> factors;
    std::vector<Trajectory> degrees;
    double seconds = 0; int sweeps = 0;
    bool complete = true; std::string note;
};

// invariant factors (ascending, each dividing the next) of a step
inline std::vector<BigInt> invariant_factors(const TorsionResult& T, const Step& s) {
    size_t n = 0; for (auto& pe : s.primary) n = std::max(n, pe.second.size());
    std::vector<BigInt> d(n, BigInt(1));          // d[0] largest
    for (auto& pe : s.primary) for (size_t i = 0; i < pe.second.size(); ++i) for (int e = 0; e < pe.second[i]; ++e) d[i] = d[i] * T.factors[pe.first].p;
    std::reverse(d.begin(), d.end());
    return d;
}

// primes: the numbers to examine (bad primes, or every candidate if classification did not run), with group flags
inline TorsionResult trajectory(const Prepared& P, const Certificate& cert, int kmax, int threads, Cache& C,
                                const std::vector<std::pair<BigInt, bool>>& primes) {
    auto t0 = std::chrono::steady_clock::now();
    P.build_boundaries();
    TorsionResult T;
    kmax = std::min(kmax, P.D - 1);                  // H_D is free
    std::map<int, std::vector<int>> per;             // degree -> factor indices
    for (auto& pg : primes) {
        std::set<int> D = degrees_of(cert, pg.first);
        int fi = -1;
        for (int q : D) if (q <= kmax) { if (fi < 0) { fi = (int)T.factors.size(); T.factors.push_back({pg.first, pg.second, false}); } per[q].push_back(fi); }
    }
    std::set<int> degs; for (auto& kv : per) degs.insert(kv.first);
    if (degs.empty()) { T.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); return T; }
    ensure_ref(P, C, cert, degs, threads);
    std::mutex mu;
    for (auto& kv : per) {
        const int q = kv.first; const size_t nc = P.count(q + 1);
        std::vector<int> rQ(nc, 0);
        for (auto& pr : C.ref.at(q)) rQ[pr.second]++;
        for (size_t c = 1; c < nc; ++c) rQ[c] += rQ[c - 1];
        // n[f][e-1][c]: number of summands of exponent >= e of factor f after column c
        std::map<int, std::vector<std::vector<int>>> n;
        std::map<int, std::vector<int>> Sprev;       // S_{k-1}(c)
        std::vector<int> active;
        // k = 1: from cached F_p pairs or sweeps
        {
            std::vector<std::pair<uint64_t, int>> small; std::vector<int> small_f, big_f;
            for (int f : kv.second) {
                const BigInt& p = T.factors[f].p;
                if (p.fits_u64() && p.to_u64() < (1ULL << 63) && !T.factors[f].group) {
                    auto it = C.fp.find(p.to_u64());
                    if (it != C.fp.end() && it->second.count(q)) {
                        std::vector<int> r(nc, 0); for (auto& pr : it->second.at(q)) r[pr.second]++;
                        for (size_t c = 1; c < nc; ++c) r[c] += r[c - 1];
                        std::vector<int> S(nc); for (size_t c = 0; c < nc; ++c) S[c] = rQ[c] - r[c];
                        Sprev[f] = S;
                    } else { small.push_back({p.to_u64(), 1}); small_f.push_back(f); }
                } else big_f.push_back(f);
            }
            std::vector<std::function<void()>> jobs;
            auto bins = pack(small);
            std::vector<std::vector<PrimeResult>> sres(bins.size()); std::vector<PrimeResult> bres(big_f.size());
            for (size_t b = 0; b < bins.size(); ++b) jobs.push_back([&, b]() { sres[b] = sweep(P, q + 1, ModU64(modulus_of(bins[b])), bins[b]); });
            for (size_t i = 0; i < big_f.size(); ++i) jobs.push_back([&, i]() {
                const Factor& F = T.factors[big_f[i]];
                bres[i] = sweep(P, q + 1, ModBig(F.p), std::vector<std::pair<BigInt, int>>{{F.p, 1}}, F.group)[0]; });
            run_jobs(jobs, threads); T.sweeps += (int)jobs.size();
            auto from_pairs = [&](int f, const PrimeResult& R) {
                std::vector<int> r(nc, 0); for (auto& pr : R.pairs) r[pr.second]++;
                for (size_t c = 1; c < nc; ++c) r[c] += r[c - 1];
                std::vector<int> S(nc); for (size_t c = 0; c < nc; ++c) S[c] = rQ[c] - r[c];
                Sprev[f] = S; if (R.split) T.factors[f].split = true;
            };
            for (size_t b = 0; b < bins.size(); ++b) for (size_t i = 0; i < bins[b].size(); ++i) {
                int f = -1; for (size_t t = 0; t < small.size(); ++t) if (small[t].first == bins[b][i].first) f = small_f[t];
                from_pairs(f, sres[b][i]);
            }
            for (size_t i = 0; i < big_f.size(); ++i) from_pairs(big_f[i], bres[i]);
            for (int f : kv.second) {
                bool any = false; for (int x : Sprev[f]) if (x) { any = true; break; }
                if (any) { n[f].push_back(Sprev[f]); active.push_back(f); }
            }
        }
        // k = 2, 3, ...: lengths modulo p^k
        for (int k = 2; !active.empty() && k <= 64; ++k) {
            std::vector<std::pair<uint64_t, int>> small; std::vector<int> small_f, big_f;
            for (int f : active) {
                const BigInt& p = T.factors[f].p;
                bool fits = !T.factors[f].group && p.fits_u64() && std::log2((double)p.to_u64()) * k < 62.9;
                if (fits) { small.push_back({p.to_u64(), k}); small_f.push_back(f); } else big_f.push_back(f);
            }
            auto bins = pack(small);
            std::vector<std::vector<PrimeResult>> sres(bins.size()); std::vector<PrimeResult> bres(big_f.size());
            std::vector<std::function<void()>> jobs;
            for (size_t b = 0; b < bins.size(); ++b) jobs.push_back([&, b]() { sres[b] = sweep(P, q + 1, ModU64(modulus_of(bins[b])), bins[b]); });
            for (size_t i = 0; i < big_f.size(); ++i) jobs.push_back([&, i, k]() {
                const Factor& F = T.factors[big_f[i]]; BigInt pk(1); for (int t = 0; t < k; ++t) pk = pk * F.p;
                bres[i] = sweep(P, q + 1, ModBig(pk), std::vector<std::pair<BigInt, int>>{{F.p, k}}, F.group)[0]; });
            run_jobs(jobs, threads); T.sweeps += (int)jobs.size();
            std::vector<int> next_active;
            auto use = [&](int f, const PrimeResult& R) {
                std::vector<int> L(nc, 0); for (auto& e : R.dlen) L[e.first] += e.second;
                for (size_t c = 1; c < nc; ++c) L[c] += L[c - 1];
                std::vector<int> S(nc), nk(nc); bool any = false;
                for (size_t c = 0; c < nc; ++c) { S[c] = k * rQ[c] - L[c]; nk[c] = S[c] - Sprev[f][c]; if (nk[c]) any = true; }
                if (R.split) T.factors[f].split = true;
                Sprev[f] = S;
                if (any) { n[f].push_back(nk); next_active.push_back(f); }
            };
            for (size_t b = 0; b < bins.size(); ++b) for (size_t i = 0; i < bins[b].size(); ++i) {
                int f = -1; for (size_t t = 0; t < small.size(); ++t) if (small[t].first == bins[b][i].first) f = small_f[t];
                use(f, sres[b][i]);
            }
            for (size_t i = 0; i < big_f.size(); ++i) use(big_f[i], bres[i]);
            active.swap(next_active);
        }
        if (!active.empty()) { T.complete = false; T.note = "exponent above 64"; }
        // steps at the level of values: the state after the last column of each value
        Trajectory tr; tr.degree = q;
        std::vector<std::pair<int, std::vector<int>>> last;  // the group is 0 before the first column
        for (size_t c = 0; c < nc; ++c) {
            if (c + 1 < nc && P.birth[q + 1][c + 1] == P.birth[q + 1][c]) continue;
            std::vector<std::pair<int, std::vector<int>>> st;
            for (auto& fv : n) {
                std::vector<int> ex;                  // exponents, descending
                const auto& L = fv.second; int E = (int)L.size();
                for (int e = E; e >= 1; --e) {
                    int exact = L[e - 1][c] - (e < E ? L[e][c] : 0);
                    for (int t = 0; t < exact; ++t) ex.push_back(e);
                }
                if (!ex.empty()) st.push_back({fv.first, ex});
            }
            if (st != last) { tr.steps.push_back({P.birth[q + 1][c], st}); last = st; }
        }
        std::lock_guard<std::mutex> lk(mu);
        T.degrees.push_back(std::move(tr));
    }
    for (auto& F : T.factors) if (F.split) { T.complete = false; T.note = "a composite factor was split by the sweep; its exponents are not reliable"; }
    T.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return T;
}

// the numbers whose torsion to follow: the bad primes once classified (Theorem 4: a prime with absolute torsion
// is bad), every candidate otherwise; with group flags
inline std::vector<std::pair<BigInt, bool>> primes_to_follow(const Certificate& cert) {
    std::vector<std::pair<BigInt, bool>> pr;
    auto isgrp = [&](const BigInt& x) { for (auto& g : cert.groups) if (g == x) return true; return false; };
    if (cert.classified) for (auto& b : cert.bad) pr.push_back({b.is_big ? b.big : BigInt::from_u64(b.p), b.group});
    else { for (auto& kv : cert.candidates) pr.push_back({BigInt::from_u64(kv.first), false});
           for (auto& kv : cert.big_candidates) pr.push_back({kv.first, isgrp(kv.first)}); }
    return pr;
}

} // namespace torsion

// the explicit complex behind a backend, for the Howell sweeps
inline const Prepared* prepared_of(const ExplicitBackend& b) { return &b.P; }
inline const Prepared* prepared_of(const ReducedBackend& b) { return &b.P; }
template <class B> inline const Prepared* prepared_of(const B&) { return nullptr; }

} // namespace perseus
