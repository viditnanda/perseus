// magnitude.hpp -- integral magnitude homology of finite graphs (Hepworth-Willerton).
//
// For a graph G with path metric d (directed or undirected; d = infinity between unreachable vertices),
// MC_{k,l}(G) is free on the tuples (x_0, ..., x_k) with x_i != x_{i+1}, all d(x_i, x_{i+1}) finite and
// sum d(x_i, x_{i+1}) = l. The differential is d = sum_{i=1}^{k-1} (-1)^i d_i, where d_i deletes x_i when
// d(x_{i-1}, x_{i+1}) = d(x_{i-1}, x_i) + d(x_i, x_{i+1}) and is 0 otherwise. Deleting an interior point keeps
// the endpoints and the length, so MC_{*,l} = sum over (s, t) of complexes MC_{*,l}(s, t); each is reduced
// separately (algebraic Morse reduction over Z, then Smith forms of the residual), in parallel.
#pragma once
#include "chain.hpp"
#include <vector>
#include <map>
#include <unordered_map>
#include <thread>
#include <atomic>
#include <mutex>
#include <chrono>
#include <deque>
#include <functional>

namespace perseus {

struct MagGraph {
    int n = 0; bool directed = false;
    std::vector<std::vector<int>> out;
    std::vector<int> dist;                       // n x n, -1 = unreachable
    MagGraph(int n_, const std::vector<std::pair<int, int>>& edges, bool directed_) : n(n_), directed(directed_), out(n_) {
        for (auto& e : edges) {
            if (e.first < 0 || e.second < 0 || e.first >= n || e.second >= n) throw std::invalid_argument("magnitude: edge endpoint out of range");
            if (e.first == e.second) continue;
            out[e.first].push_back(e.second);
            if (!directed) out[e.second].push_back(e.first);
        }
        for (auto& v : out) { std::sort(v.begin(), v.end()); v.erase(std::unique(v.begin(), v.end()), v.end()); }
        dist.assign((size_t)n * n, -1);
        std::deque<int> q;
        for (int s = 0; s < n; ++s) {
            int* D = &dist[(size_t)s * n]; D[s] = 0; q.assign(1, s);
            while (!q.empty()) { int u = q.front(); q.pop_front(); for (int v : out[u]) if (D[v] < 0) { D[v] = D[u] + 1; q.push_back(v); } }
        }
    }
    int d(int a, int b) const { return dist[(size_t)a * n + b]; }
};

struct MagEntry {                                // MH_{k,l} summed over the pairs computed
    size_t generators = 0, free = 0;
    std::map<std::string, size_t> torsion;       // invariant factor (decimal) -> multiplicity
    size_t residual = 0;                         // cells left after the Morse reduction
};
struct MagResult {
    int lmax = 0;
    std::vector<std::vector<MagEntry>> table;    // table[l][k]
    double seconds = 0;
    size_t pairs = 0, max_complex = 0;
};

// MC_{*,l}(s, t) for l = 0..lmax: tuples from s to t, stored flat by length and size: byl[l][k] holds the
// (k+1)-tuples of length l back to back
typedef std::vector<std::vector<std::vector<int32_t>>> FlatTuples;
inline void magnitude_tuples(const MagGraph& G, int s, int t, int lmax, FlatTuples& byl) {
    byl.assign(lmax + 1, std::vector<std::vector<int32_t>>());
    std::vector<int32_t> cur{(int32_t)s};
    std::function<void(int, int)> go = [&](int u, int used) {
        if (u == t) {
            auto& L = byl[used]; size_t k = cur.size() - 1;
            if (L.size() <= k) L.resize(k + 1);
            L[k].insert(L[k].end(), cur.begin(), cur.end());
        }
        for (int v = 0; v < G.n; ++v) {
            if (v == u) continue;
            int duv = G.d(u, v), dvt = G.d(v, t);
            if (duv <= 0 || dvt < 0 || used + duv + dvt > lmax) continue;
            cur.push_back(v); go(v, used + duv); cur.pop_back();
        }
    };
    if (G.d(s, t) >= 0 && G.d(s, t) <= lmax) go(s, 0);
}

// the chain complex of the tuples of one length l and one pair (s, t), given flat by size; cells grouped by
// k and sorted lexicographically within a size, so a face is found by binary search
inline FilteredComplex magnitude_complex(const MagGraph& G, std::vector<std::vector<int32_t>> byk) {
    std::vector<size_t> base(byk.size() + 1, 0);
    for (size_t k = 0; k < byk.size(); ++k) {                  // sort the (k+1)-tuples of each size
        size_t w = k + 1, m = byk[k].size() / w;
        std::vector<uint32_t> ix(m); for (size_t i = 0; i < m; ++i) ix[i] = (uint32_t)i;
        const int32_t* D = byk[k].data();
        std::sort(ix.begin(), ix.end(), [&](uint32_t a, uint32_t b) { return std::lexicographical_compare(D + a * w, D + a * w + w, D + b * w, D + b * w + w); });
        std::vector<int32_t> S; S.reserve(byk[k].size());
        for (uint32_t i : ix) S.insert(S.end(), D + (size_t)i * w, D + (size_t)i * w + w);
        byk[k].swap(S);
        base[k + 1] = base[k] + m;
    }
    FilteredComplex C;
    std::vector<std::pair<int32_t, int64_t>> bd; std::vector<int32_t> f;
    for (size_t k = 0; k < byk.size(); ++k) {
        size_t w = k + 1, m = byk[k].size() / w;
        for (size_t j = 0; j < m; ++j) {
            const int32_t* x = byk[k].data() + j * w;
            bd.clear();
            for (size_t i = 1; i + 1 <= k; ++i) {
                if (G.d(x[i - 1], x[i]) + G.d(x[i], x[i + 1]) != G.d(x[i - 1], x[i + 1])) continue;
                f.assign(x, x + w); f.erase(f.begin() + i);
                const std::vector<int32_t>& L = byk[k - 1]; size_t lo = 0, hi = L.size() / k;   // faces have size k
                while (lo < hi) { size_t mid = (lo + hi) / 2; if (std::lexicographical_compare(L.data() + mid * k, L.data() + mid * k + k, f.begin(), f.end())) lo = mid + 1; else hi = mid; }
                if (lo == L.size() / k || !std::equal(f.begin(), f.end(), L.data() + lo * k)) throw std::logic_error("magnitude_complex: a face is missing");
                bd.push_back({(int32_t)(base[k - 1] + lo), (i & 1) ? -1 : 1});
            }
            C.add_cell((int)k, 0.0, bd);
        }
    }
    return C;
}

// The lean route for large complexes. The tuples of one (s, t, l) come out of magnitude_tuples in lexicographic
// order within each size (depth-first search with ascending neighbours), so faces are found by binary search
// without sorting. The complex is stored once, in compressed rows (int32 faces, int8 coefficients +-1: deleting
// two different interior points never gives the same tuple) with the compressed coboundary, and the
// fill-in-free part of the algebraic Morse reduction (collapses and coreductions, which only delete cells) runs on
// it with live face and coface counts. Only the survivors are copied into a FilteredComplex for the general
// reduction and the Smith forms. Peak memory is about a third of building the FilteredComplex first.
inline FilteredComplex magnitude_survivors(const MagGraph& G, std::vector<std::vector<int32_t>>& byk, size_t& cells_out) {
    size_t K = byk.size();
    std::vector<size_t> base(K + 1, 0);
    for (size_t k = 0; k < K; ++k) {
        size_t w = k + 1, m = byk[k].size() / w;
        const int32_t* D = byk[k].data();
        bool sorted = true;
        for (size_t i = 1; i < m && sorted; ++i) if (std::lexicographical_compare(D + i * w, D + i * w + w, D + (i - 1) * w, D + (i - 1) * w + w)) sorted = false;
        if (!sorted) {
            std::vector<uint32_t> ix(m); for (size_t i = 0; i < m; ++i) ix[i] = (uint32_t)i;
            std::sort(ix.begin(), ix.end(), [&](uint32_t a, uint32_t b) { return std::lexicographical_compare(D + (size_t)a * w, D + (size_t)a * w + w, D + (size_t)b * w, D + (size_t)b * w + w); });
            std::vector<int32_t> S; S.reserve(byk[k].size());
            for (uint32_t i : ix) S.insert(S.end(), D + (size_t)i * w, D + (size_t)i * w + w);
            byk[k].swap(S);
        }
        base[k + 1] = base[k] + m;
    }
    size_t N = base[K]; cells_out = N;
    if (N > (size_t)INT32_MAX) throw std::runtime_error("magnitude: more than 2^31 tuples in one complex");
    std::vector<uint64_t> off(N + 1, 0); std::vector<int32_t> face; std::vector<int8_t> coef;
    std::vector<int> dimv(N);
    std::vector<int32_t> f;
    for (size_t k = 0; k < K; ++k) {
        size_t w = k + 1, m = byk[k].size() / w;
        for (size_t j = 0; j < m; ++j) {
            const int32_t* x = byk[k].data() + j * w;
            size_t c = base[k] + j; dimv[c] = (int)k;
            for (size_t i = 1; i + 1 <= k; ++i) {
                if (G.d(x[i - 1], x[i]) + G.d(x[i], x[i + 1]) != G.d(x[i - 1], x[i + 1])) continue;
                f.assign(x, x + w); f.erase(f.begin() + i);
                const std::vector<int32_t>& L = byk[k - 1]; size_t lo = 0, hi = L.size() / k;
                while (lo < hi) { size_t mid = (lo + hi) / 2; if (std::lexicographical_compare(L.data() + mid * k, L.data() + mid * k + k, f.begin(), f.end())) lo = mid + 1; else hi = mid; }
                if (lo == L.size() / k || !std::equal(f.begin(), f.end(), L.data() + lo * k)) throw std::logic_error("magnitude_complex: a face is missing");
                face.push_back((int32_t)(base[k - 1] + lo)); coef.push_back((i & 1) ? -1 : 1);
            }
            off[c + 1] = face.size();
        }
    }
    for (auto& v : byk) { std::vector<int32_t>().swap(v); }
    std::vector<uint64_t> coff(N + 1, 0);
    for (int32_t x : face) ++coff[x + 1];
    for (size_t c = 0; c < N; ++c) coff[c + 1] += coff[c];
    std::vector<int32_t> cob(face.size());
    { std::vector<uint64_t> pos(coff.begin(), coff.end() - 1); for (size_t c = 0; c < N; ++c) for (uint64_t a = off[c]; a < off[c + 1]; ++a) cob[pos[face[a]]++] = (int32_t)c; }
    std::vector<char> alive(N, 1), queued(N, 1);
    std::vector<uint32_t> nf(N), nc(N);
    for (size_t c = 0; c < N; ++c) { nf[c] = (uint32_t)(off[c + 1] - off[c]); nc[c] = (uint32_t)(coff[c + 1] - coff[c]); }
    std::vector<int32_t> work(N); for (size_t c = 0; c < N; ++c) work[c] = (int32_t)(N - 1 - c);
    auto push = [&](int32_t c) { if (alive[c] && !queued[c]) { queued[c] = 1; work.push_back(c); } };
    auto kill = [&](int32_t sigma, int32_t tau) {
        alive[sigma] = alive[tau] = 0;
        for (int32_t c : {sigma, tau}) {
            for (uint64_t a = off[c]; a < off[c + 1]; ++a) if (alive[face[a]]) { --nc[face[a]]; push(face[a]); }
            for (uint64_t a = coff[c]; a < coff[c + 1]; ++a) if (alive[cob[a]]) { --nf[cob[a]]; push(cob[a]); }
        }
    };
    while (!work.empty()) {
        int32_t c = work.back(); work.pop_back(); queued[c] = 0;
        if (!alive[c]) continue;
        if (nf[c] == 1) {                              // coreduction: the single live face
            for (uint64_t a = off[c]; a < off[c + 1]; ++a) if (alive[face[a]]) { kill(face[a], c); break; }
        } else if (nc[c] == 1) {                       // collapse: the single live coface
            for (uint64_t a = coff[c]; a < coff[c + 1]; ++a) if (alive[cob[a]]) { kill(c, cob[a]); break; }
        }
    }
    std::vector<int32_t>().swap(cob); std::vector<uint64_t>().swap(coff);
    std::vector<int32_t> nid(N, -1); FilteredComplex C;
    std::vector<std::pair<int32_t, int64_t>> bd;
    for (size_t c = 0; c < N; ++c) {
        if (!alive[c]) continue;
        bd.clear();
        for (uint64_t a = off[c]; a < off[c + 1]; ++a) if (alive[face[a]]) bd.push_back({nid[face[a]], coef[a]});
        nid[c] = (int32_t)C.size(); C.add_cell(dimv[c], 0.0, bd);
    }
    return C;
}

// MH_{k,l}(G) for l <= lmax (all k), over the given (s, t) pairs (all pairs if empty), with `threads` workers
inline MagResult magnitude_homology(const MagGraph& G, int lmax, int threads = 1, std::vector<std::pair<int, int>> pairs = {}) {
    auto T0 = std::chrono::steady_clock::now();
    if (pairs.empty()) for (int s = 0; s < G.n; ++s) for (int t = 0; t < G.n; ++t) if (G.d(s, t) >= 0 && G.d(s, t) <= lmax) pairs.push_back({s, t});
    MagResult R; R.lmax = lmax; R.pairs = pairs.size();
    R.table.assign(lmax + 1, std::vector<MagEntry>(lmax + 1));
    std::mutex mu; std::atomic<size_t> next(0);
    std::string err;
    auto worker = [&]() {
        FlatTuples byl;
        for (size_t i; (i = next++) < pairs.size();) {
            try {
                magnitude_tuples(G, pairs[i].first, pairs[i].second, lmax, byl);
                for (int l = 0; l <= lmax; ++l) {
                    if (byl[l].empty()) continue;
                    std::vector<size_t> gens(lmax + 1, 0);
                    for (size_t k = 0; k < byl[l].size() && k <= (size_t)lmax; ++k) gens[k] = byl[l][k].size() / (k + 1);
                    int top = 0; for (size_t k = 0; k < byl[l].size(); ++k) if (!byl[l][k].empty()) top = (int)k;
                    size_t ncells = 0;
                    FilteredComplex C = magnitude_survivors(G, byl[l], ncells);
                    ReducedComplex Rc = algebraic_morse(C);
                    IntegralHomology H = integral_homology(Rc.C, top);
                    std::vector<size_t> res(lmax + 1, 0); for (int dd : Rc.C.dim) ++res[dd];
                    std::lock_guard<std::mutex> lk(mu);
                    R.max_complex = std::max(R.max_complex, ncells);
                    for (int k = 0; k <= top && k <= lmax; ++k) {
                        MagEntry& E = R.table[l][k];
                        E.generators += gens[k]; E.residual += res[k];
                        E.free += H.free[k];
                        for (auto& x : H.torsion[k]) E.torsion[x.str()]++;
                    }
                }
            } catch (std::exception& e) { std::lock_guard<std::mutex> lk(mu); if (err.empty()) err = e.what(); }
        }
    };
    int nt = std::max(1, threads);
    std::vector<std::thread> pool; for (int i = 1; i < nt; ++i) pool.emplace_back(worker);
    worker(); for (auto& th : pool) th.join();
    if (!err.empty()) throw std::runtime_error(err);
    R.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - T0).count();
    return R;
}

// rank of MH_{k,l} over F_p from the integral answer (universal coefficients)
inline size_t magnitude_rank_mod(const MagResult& R, int k, int l, uint64_t p) {
    auto divp = [&](const std::map<std::string, size_t>& T) {
        size_t c = 0; for (auto& kv : T) { BigInt x; /* parse */ for (char ch : kv.first) x = x * BigInt(10) + BigInt(ch - '0'); if (x.mod_small(p) == 0) c += kv.second; }
        return c;
    };
    const MagEntry& E = R.table[l][k];
    size_t r = E.free + divp(E.torsion);
    if (k > 0) r += divp(R.table[l][k - 1].torsion);
    return r;
}

} // namespace perseus
