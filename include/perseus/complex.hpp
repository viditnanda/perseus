// complex.hpp -- explicit filtered cell complexes with integer boundary coefficients,
// and builders for simplicial and cubical complexes.
#pragma once
#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <limits>
#include <cmath>
#include <stdexcept>

namespace perseus {

struct FilteredComplex {
    std::vector<int> dim;
    std::vector<double> birth;
    std::vector<int64_t> off{0};     // CSR boundary: faces of cell c are face[off[c] .. off[c+1])
    std::vector<int32_t> face;
    std::vector<int64_t> coef;

    size_t size() const { return dim.size(); }
    int max_dim() const { int d = -1; for (int x : dim) d = std::max(d, x); return d; }
    int add_cell(int d, double b, const std::vector<std::pair<int32_t, int64_t>>& bd) {
        int id = (int)dim.size();
        dim.push_back(d); birth.push_back(b);
        for (auto& fc : bd) if (fc.second != 0) { face.push_back(fc.first); coef.push_back(fc.second); }
        off.push_back((int64_t)face.size());
        return id;
    }
    // checks: faces exist, have dimension d-1, are born no later; optionally d^2 = 0 over Z
    std::string validate(bool check_d2 = false) const {
        for (size_t c = 0; c < size(); ++c) {
            for (int64_t e = off[c]; e < off[c + 1]; ++e) {
                int32_t f = face[e];
                if (f < 0 || (size_t)f >= size()) return "cell " + std::to_string(c) + ": face index out of range";
                if (dim[f] != dim[c] - 1) return "cell " + std::to_string(c) + ": face of wrong dimension";
                if (birth[f] > birth[c]) return "cell " + std::to_string(c) + ": born before its face " + std::to_string(f);
            }
        }
        if (check_d2) {
            std::unordered_map<int32_t, int64_t> acc;
            for (size_t c = 0; c < size(); ++c) {
                acc.clear();
                for (int64_t e = off[c]; e < off[c + 1]; ++e) {
                    int32_t f = face[e];
                    for (int64_t g = off[f]; g < off[f + 1]; ++g) acc[face[g]] += coef[e] * coef[g];
                }
                for (auto& kv : acc) if (kv.second != 0) return "boundary of boundary nonzero at cell " + std::to_string(c);
            }
        }
        return "";
    }
};

// ---- simplicial complexes given by (not necessarily maximal) simplices with births ----------
// Faces inherit the minimum birth over all listed simplices containing them (Perseus convention).
struct VecHash {
    size_t operator()(const std::vector<int>& v) const {
        uint64_t h = 1469598103934665603ULL;
        for (int x : v) { h ^= (uint64_t)(uint32_t)x + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2); }
        return (size_t)h;
    }
};

class SimplicialBuilder {
public:
    int maxdim_cap = std::numeric_limits<int>::max();   // drop simplices above this dimension
    // add simplex (vertex labels arbitrary ints) and all its faces
    void add(std::vector<int> v, double b) {
        std::sort(v.begin(), v.end()); v.erase(std::unique(v.begin(), v.end()), v.end());
        int d = (int)v.size() - 1;
        if (d < 0) return;
        if (d > maxdim_cap) {           // add only the faces of dimension maxdim_cap
            std::vector<int> sub; choose_add(v, 0, maxdim_cap + 1, sub, b); return;
        }
        insert_with_faces(v, b);
    }
    FilteredComplex build() const {
        // order: by dimension, then lexicographically (deterministic ids)
        std::vector<std::pair<const std::vector<int>*, double>> all;
        all.reserve(births_.size());
        for (auto& kv : births_) all.push_back({&kv.first, kv.second});
        std::sort(all.begin(), all.end(), [](const std::pair<const std::vector<int>*, double>& x, const std::pair<const std::vector<int>*, double>& y) {
            if (x.first->size() != y.first->size()) return x.first->size() < y.first->size();
            return *x.first < *y.first;
        });
        std::unordered_map<std::vector<int>, int32_t, VecHash> id; id.reserve(all.size() * 2);
        FilteredComplex C;
        std::vector<std::pair<int32_t, int64_t>> bd; std::vector<int> f;
        for (auto& pr : all) {
            const std::vector<int>& s = *pr.first;
            bd.clear();
            if (s.size() > 1) for (size_t i = 0; i < s.size(); ++i) {
                f.assign(s.begin(), s.end()); f.erase(f.begin() + i);
                bd.push_back({id.at(f), (i & 1) ? -1 : 1});
            }
            id[s] = C.add_cell((int)s.size() - 1, pr.second, bd);
        }
        return C;
    }
    size_t size() const { return births_.size(); }
private:
    std::unordered_map<std::vector<int>, double, VecHash> births_;
    void choose_add(const std::vector<int>& v, size_t start, int k, std::vector<int>& cur, double b) {
        if ((int)cur.size() == k) { insert_with_faces(cur, b); return; }
        for (size_t i = start; i < v.size(); ++i) { cur.push_back(v[i]); choose_add(v, i + 1, k, cur, b); cur.pop_back(); }
    }
    void insert_with_faces(const std::vector<int>& v, double b) {
        auto it = births_.find(v);
        if (it != births_.end()) {
            if (it->second <= b) return;           // all faces already at least as old
            it->second = b;
        } else births_.emplace(v, b);
        if (v.size() <= 1) return;
        std::vector<int> f;
        for (size_t i = 0; i < v.size(); ++i) { f.assign(v.begin(), v.end()); f.erase(f.begin() + i); insert_with_faces(f, b); }
    }
};

// ---- cubical complexes on a grid --------------------------------------------------------
// Cells are addressed in doubled ("Khalimsky") coordinates: an elementary cube has coordinate
// c_i odd along the axes where it has extent, even elsewhere. Top cubes have all-odd coordinates.
// ext[i] = number of top cubes along axis i; periodic[i] wraps axis i.
// birth of a cell = min birth of the top cubes containing it; absent top cubes have birth +inf.
inline FilteredComplex build_cubical(const std::vector<int64_t>& ext, const std::vector<bool>& periodic,
                                     const std::vector<double>& top_birth /* lexicographic, axis 0 fastest */) {
    int d = (int)ext.size();
    std::vector<int64_t> L(d), stride(d);
    int64_t total = 1;
    for (int i = 0; i < d; ++i) { L[i] = periodic[i] ? 2 * ext[i] : 2 * ext[i] + 1; stride[i] = total; total *= L[i]; }
    const double INF = std::numeric_limits<double>::infinity();
    std::vector<double> b(total, INF);
    // place top cubes
    {
        std::vector<int64_t> a(d, 0);
        for (size_t t = 0; t < top_birth.size(); ++t) {
            int64_t pos = 0; for (int i = 0; i < d; ++i) pos += (2 * a[i] + 1) * stride[i];
            b[pos] = top_birth[t];
            for (int i = 0; i < d; ++i) { if (++a[i] < ext[i]) break; a[i] = 0; }
        }
    }
    // separable min propagation along each axis (even coordinates take min of odd neighbours)
    for (int ax = 0; ax < d; ++ax) {
        for (int64_t pos = 0; pos < total; ++pos) {
            int64_t c = (pos / stride[ax]) % L[ax];
            if (c & 1) continue;
            double m = INF;
            int64_t lo = c - 1, hi = c + 1;
            if (periodic[ax]) { lo = (lo + L[ax]) % L[ax]; hi = hi % L[ax]; }
            if (lo >= 0) m = std::min(m, b[pos + (lo - c) * stride[ax]]);
            if (hi < L[ax]) m = std::min(m, b[pos + (hi - c) * stride[ax]]);
            b[pos] = m;
        }
    }
    // the propagation above must read already-final odd values: since axis passes only write even
    // coordinates along the current axis and read odd ones along it, one pass per axis suffices.
    // assign ids by dimension
    std::vector<int> cdim(total);
    std::vector<int32_t> id(total, -1);
    std::vector<std::vector<int64_t>> bydim(d + 1);
    for (int64_t pos = 0; pos < total; ++pos) {
        if (b[pos] == INF) continue;
        int k = 0; for (int i = 0; i < d; ++i) k += ((pos / stride[i]) % L[i]) & 1;
        cdim[pos] = k; bydim[k].push_back(pos);
    }
    FilteredComplex C;
    std::vector<std::pair<int32_t, int64_t>> bd;
    for (int k = 0; k <= d; ++k) {
        for (int64_t pos : bydim[k]) {
            bd.clear(); int j = 0;
            for (int i = 0; i < d; ++i) {
                int64_t c = (pos / stride[i]) % L[i];
                if (!(c & 1)) continue;
                int64_t lo = c - 1, hi = c + 1;
                if (periodic[i]) { hi %= L[i]; }
                int64_t plo = pos + (lo - c) * stride[i], phi = pos + (hi - c) * stride[i];
                int64_t s = (j & 1) ? -1 : 1;
                // d Q = sum_j (-1)^j (Q_j^+ - Q_j^-); for a periodic axis of length 1 the two faces coincide
                bd.push_back({id[phi], s}); bd.push_back({id[plo], -s});
                ++j;
            }
            // merge duplicates (periodic axes of extent 1)
            std::sort(bd.begin(), bd.end());
            std::vector<std::pair<int32_t, int64_t>> m;
            for (auto& e : bd) { if (!m.empty() && m.back().first == e.first) m.back().second += e.second; else m.push_back(e); }
            id[pos] = C.add_cell(k, b[pos], m);
        }
    }
    return C;
}

// sparse cubical complex from top-cube anchors (integer coordinates) with births
inline FilteredComplex build_cubical_sparse(int d, const std::vector<std::vector<int64_t>>& anchors, const std::vector<double>& births) {
    std::map<std::vector<int64_t>, double> b;   // doubled coordinates -> birth
    std::vector<int64_t> c(d);
    for (size_t t = 0; t < anchors.size(); ++t) {
        // all 3^d faces of the cube anchored at a: coordinates 2a_i + {0,1,2}
        int64_t n3 = 1; for (int i = 0; i < d; ++i) n3 *= 3;
        for (int64_t s = 0; s < n3; ++s) {
            int64_t r = s; for (int i = 0; i < d; ++i) { c[i] = 2 * anchors[t][i] + r % 3; r /= 3; }
            auto it = b.find(c);
            if (it == b.end()) b.emplace(c, births[t]); else it->second = std::min(it->second, births[t]);
        }
    }
    std::vector<std::vector<const std::pair<const std::vector<int64_t>, double>*>> bydim(d + 1);
    for (auto& kv : b) { int k = 0; for (int i = 0; i < d; ++i) k += (int)(kv.first[i] & 1); bydim[k].push_back(&kv); }
    std::map<std::vector<int64_t>, int32_t> id;
    FilteredComplex C;
    std::vector<std::pair<int32_t, int64_t>> bd;
    for (int k = 0; k <= d; ++k) for (auto* kv : bydim[k]) {
        bd.clear(); int j = 0; std::vector<int64_t> f = kv->first;
        for (int i = 0; i < d; ++i) {
            if (!(kv->first[i] & 1)) continue;
            int64_t s = (j & 1) ? -1 : 1;
            f[i] = kv->first[i] + 1; bd.push_back({id.at(f), s});
            f[i] = kv->first[i] - 1; bd.push_back({id.at(f), -s});
            f[i] = kv->first[i]; ++j;
        }
        id[kv->first] = C.add_cell(k, kv->second, bd);
    }
    return C;
}

} // namespace perseus
