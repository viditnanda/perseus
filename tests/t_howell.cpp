// t_howell: the Howell sweep of howell.hpp against first principles on random integer matrices.
//  * squarefree m: for every prime p | m, the events "p stops dividing d_r at column c" are exactly the
//    persistence pairs (low, column) of a plain column reduction over F_p;
//  * m = prod p^k: after every column, sum_r (k - min(k, v_p(d_r))) = log_p |span mod p^k|, the span
//    enumerated by brute force.
// Both the one-word (ModU64) and the BigInt (ModBig) arithmetic are tested.
#include "perseus/howell.hpp"
#include <cstdio>
#include <cstdlib>
#include <random>
#include <set>
#include <map>
#include <unordered_set>
#include <cmath>

using namespace perseus;
using namespace perseus::howell;

typedef std::vector<std::vector<std::pair<int, int64_t>>> Cols;   // per column: (row, coef), rows ascending

static std::set<std::pair<int, int>> fp_pairs(const Cols& A, int nrows, uint64_t p) {
    std::vector<int> owner(nrows, -1); std::vector<std::map<int, uint64_t>> red(A.size());
    std::set<std::pair<int, int>> out;
    auto md = [&](int64_t c) { int64_t r = c % (int64_t)p; return (uint64_t)(r < 0 ? r + (int64_t)p : r); };
    for (size_t j = 0; j < A.size(); ++j) {
        std::map<int, uint64_t> c; for (auto& e : A[j]) { uint64_t v = md(e.second); if (v) c[e.first] = v; }
        while (!c.empty()) {
            int r = c.rbegin()->first; if (owner[r] < 0) break;
            auto& w = red[owner[r]]; uint64_t f = (uint64_t)((u128)c[r] * ModU64::inv_mod(w.at(r), p) % p);
            for (auto& kv : w) { uint64_t s = (uint64_t)((c.count(kv.first) ? c[kv.first] : 0) + p - (u128)f * kv.second % p) % p;
                                 if (s) c[kv.first] = s; else c.erase(kv.first); }
        }
        if (!c.empty()) { int r = c.rbegin()->first; owner[r] = (int)j; red[j] = c; out.insert({r, (int)j}); }
    }
    return out;
}

// log_p |span of the first c columns mod q = p^k|, by enumeration (q^nrows small)
static int span_log(const Cols& A, size_t c, int nrows, uint64_t p, int k) {
    uint64_t q = 1; for (int i = 0; i < k; ++i) q *= p;
    auto enc = [&](const std::vector<uint64_t>& v) { uint64_t x = 0; for (int i = nrows; i-- > 0;) x = x * q + v[i]; return x; };
    std::unordered_set<uint64_t> S; S.insert(0);
    for (size_t j = 0; j < c; ++j) {
        std::vector<uint64_t> g(nrows, 0);
        for (auto& e : A[j]) { int64_t r = e.second % (int64_t)q; g[e.first] = (uint64_t)(r < 0 ? r + (int64_t)q : r); }
        std::vector<uint64_t> cur(S.begin(), S.end());
        for (uint64_t s : cur) {
            std::vector<uint64_t> v(nrows); uint64_t x = s; for (int i = 0; i < nrows; ++i) { v[i] = x % q; x /= q; }
            for (uint64_t i = 1; i < q; ++i) { for (int t = 0; t < nrows; ++t) v[t] = (v[t] + g[t]) % q; S.insert(enc(v)); }
        }
    }
    int l = 0; size_t n = S.size(); while (n > 1) { n /= p; ++l; } return l;
}

template <class M, class Conv>
static int check(const M& R, Conv conv, const Cols& A, int nrows, const std::vector<std::pair<uint64_t, int>>& pk, bool squarefree) {
    Sweep<M> S(R, nrows);
    int bad = 0;
    std::vector<std::vector<typename M::T>> dcur(pk.size());
    for (size_t j = 0; j < A.size(); ++j) {
        std::vector<typename Sweep<M>::Ent> x;
        for (auto& e : A[j]) { auto v = R.from_int(e.second); if (!M::is_zero(v)) x.push_back({e.first, v}); }
        S.insert(x, (int32_t)j);
        if (!squarefree) {
            for (size_t i = 0; i < pk.size(); ++i) {
                uint64_t p = pk[i].first; int k = pk[i].second; int L = 0;
                for (int r = 0; r < nrows; ++r) L += S.piv[r].empty() ? 0 : k - M::val(S.d[r], conv(p), k);
                int want = span_log(A, j + 1, nrows, p, k);
                if (L != want) { ++bad; if (bad < 5) fprintf(stderr, "length mismatch p=%llu k=%d col %zu: %d vs %d\n", (unsigned long long)p, k, j, L, want); }
            }
        }
    }
    if (squarefree) {
        for (auto& pr : pk) {
            uint64_t p = pr.first; std::set<std::pair<int, int>> got;
            for (auto& ev : S.events) if (M::val(ev.dold, conv(p), 1) == 1 && M::val(ev.dnew, conv(p), 1) == 0) got.insert({ev.row, ev.col});
            auto want = fp_pairs(A, nrows, p);
            if (got != want) { ++bad; if (bad < 5) fprintf(stderr, "pair mismatch p=%llu: %zu vs %zu pairs\n", (unsigned long long)p, got.size(), want.size()); }
        }
    }
    return bad;
}

int main(int argc, char** argv) {
    int iters = argc > 1 ? atoi(argv[1]) : 300;
    std::mt19937_64 rng(12345);
    const uint64_t primes[] = {2, 3, 5, 7, 11, 13, 2305843009213693951ULL};   // the last is 2^61 - 1
    int bad = 0, nsq = 0, npk = 0;
    for (int it = 0; it < iters; ++it) {
        bool sq = it % 2 == 0;
        int nrows = sq ? 2 + rng() % 7 : 1 + rng() % 4, ncols = 1 + rng() % 10;
        Cols A(ncols);
        for (int j = 0; j < ncols; ++j) {
            int scale = (int)(rng() % 3 == 0 ? 1 + rng() % 12 : 1);
            for (int r = 0; r < nrows; ++r) if (rng() % 3 == 0) { int64_t v = (int64_t)(rng() % 9) - 4; if (v) A[j].push_back({r, v * scale}); }
            if (j > 1 && rng() % 4 == 0) {            // a combination of earlier columns: a relation over Q
                std::map<int, int64_t> c; for (int t = 0; t < 2; ++t) { int jj = (int)(rng() % j); int64_t f = (int64_t)(rng() % 5) - 2; for (auto& e : A[jj]) c[e.first] += f * e.second; }
                A[j].clear(); for (auto& kv : c) if (kv.second) A[j].push_back(kv);
            }
        }
        std::vector<std::pair<uint64_t, int>> pk;
        if (sq) {
            u128 m = 1; for (uint64_t p : primes) if (rng() % 2 && (p < 100 || m < 100)) { if (m * p >= ((u128)1 << 62) && p > 100) continue; pk.push_back({p, 1}); m *= p; }
            if (pk.empty()) { pk.push_back({2, 1}); m = 2; }
            if (m < ((u128)1 << 62)) bad += check(ModU64((uint64_t)m), [](uint64_t p) { return p; }, A, nrows, pk, true);
            BigInt mb(1); for (auto& x : pk) mb = mb * BigInt::from_u64(x.first);
            bad += check(ModBig(mb), [](uint64_t p) { return BigInt::from_u64(p); }, A, nrows, pk, true);
            ++nsq;
        } else {
            uint64_t m = 1;
            for (uint64_t p : {2ULL, 3ULL, 5ULL}) {
                if (rng() % 2 == 0) continue;
                int k = 1 + (int)(rng() % (p == 2 ? 3 : 2));
                uint64_t q = 1;
                for (int i = 0; i < k; ++i) q *= p;
                if (std::pow((double)q, nrows) > 2e5) continue;
                pk.push_back({p, k}); m *= q;
            }
            if (pk.empty()) { pk.push_back({2, 2}); m = 4; }
            bad += check(ModU64(m), [](uint64_t p) { return p; }, A, nrows, pk, false);
            bad += check(ModBig(BigInt::from_u64(m)), [](uint64_t p) { return BigInt::from_u64(p); }, A, nrows, pk, false);
            ++npk;
        }
    }
    printf("t_howell: %d squarefree and %d prime-power cases, %d mismatches\n", nsq, npk, bad);
    return bad ? 1 : 0;
}
