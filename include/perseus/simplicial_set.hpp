// simplicial_set.hpp -- finite filtered simplicial sets and their normalised chain complexes.
//
// A simplicial set is given by its nondegenerate simplices, each with a dimension, a birth and (in
// dimension n >= 1) its n + 1 faces. A face is a simplex s_J y: a nondegenerate y and a degeneracy word,
// kept in the canonical form s_{j1} s_{j2} ... s_{jm} with j1 > j2 > ... > jm (every composite of
// degeneracies has exactly one such form, by s_i s_j = s_{j+1} s_i for i <= j). Faces of degenerate
// simplices follow from
//     d_i s_j = s_{j-1} d_i (i < j),   d_i s_j = id (i = j, j + 1),   d_i s_j = s_j d_{i-1} (i > j + 1).
// The normalised chain complex has the nondegenerate simplices as basis and d x = sum_i (-1)^i [d_i x],
// degenerate faces dropped; repeated faces add, so coefficients 0 and +-2 occur.
#pragma once
#include "complex.hpp"
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <functional>
#include <istream>
#include <sstream>

namespace perseus {

struct SSimp {                       // s_J y, J = (j1 > j2 > ... > jm)
    int32_t nd = -1;
    std::vector<int> J;
    bool operator==(const SSimp& o) const { return nd == o.nd && J == o.J; }
    bool operator<(const SSimp& o) const { return nd != o.nd ? nd < o.nd : J < o.J; }
    bool degenerate() const { return !J.empty(); }
};

// s_k applied to s_J y, in canonical form
inline SSimp apply_degeneracy(int k, SSimp x) {
    // s_k s_{j1} s_{j2} ... : while k <= j_r, s_k s_{j_r} = s_{j_r + 1} s_k
    std::vector<int> out;
    size_t r = 0;
    while (r < x.J.size() && k <= x.J[r]) { out.push_back(x.J[r] + 1); ++r; }
    out.push_back(k);
    for (; r < x.J.size(); ++r) out.push_back(x.J[r]);
    x.J = out; return x;
}
// s_{k1} s_{k2} ... s_{km} y (applied right to left), in canonical form
inline SSimp make_degenerate(int32_t y, const std::vector<int>& ks) {
    SSimp x; x.nd = y;
    for (size_t r = ks.size(); r-- > 0;) x = apply_degeneracy(ks[r], x);
    return x;
}

class SimplicialSet {
public:
    std::vector<int> dim;
    std::vector<double> birth;
    std::vector<std::vector<SSimp>> faces;     // faces[x][i] = d_i x (empty for vertices)
    std::vector<std::string> label;            // optional

    int32_t add(int d, double b, std::vector<SSimp> f, std::string lab = "") {
        if (d < 0) throw std::invalid_argument("simplicial set: negative dimension");
        if ((d == 0 && !f.empty()) || (d > 0 && (int)f.size() != d + 1)) throw std::invalid_argument("simplicial set: a " + std::to_string(d) + "-simplex needs " + std::to_string(d == 0 ? 0 : d + 1) + " faces");
        dim.push_back(d); birth.push_back(b); faces.push_back(std::move(f)); label.push_back(lab);
        return (int32_t)dim.size() - 1;
    }
    size_t size() const { return dim.size(); }
    int sdim(const SSimp& x) const { return dim.at(x.nd) + (int)x.J.size(); }

    // d_i of s_J y
    SSimp face(const SSimp& x, int i) const {
        if (x.J.empty()) {
            if (dim[x.nd] == 0) throw std::logic_error("simplicial set: face of a vertex");
            return faces[x.nd][i];
        }
        int j = x.J[0];
        SSimp rest = x; rest.J.erase(rest.J.begin());
        if (i < j) return apply_degeneracy(j - 1, face(rest, i));
        if (i == j || i == j + 1) return rest;
        return apply_degeneracy(j, face(rest, i - 1));
    }

    // checks: face ids exist, dimensions match, faces are born no later, degeneracy indices are in range,
    // and the simplicial identities d_i d_j = d_{j-1} d_i (i < j) hold on every nondegenerate simplex
    std::string validate() const {
        for (size_t x = 0; x < size(); ++x) {
            for (size_t i = 0; i < faces[x].size(); ++i) {
                const SSimp& f = faces[x][i];
                if (f.nd < 0 || (size_t)f.nd >= size()) return "simplex " + std::to_string(x) + ": face " + std::to_string(i) + " refers to a missing simplex";
                for (size_t r = 0; r < f.J.size(); ++r) {
                    if (r + 1 < f.J.size() && f.J[r] <= f.J[r + 1]) return "simplex " + std::to_string(x) + ": degeneracy word not canonical";
                    if (f.J[r] < 0 || f.J[r] > dim[f.nd] + (int)(f.J.size() - 1 - r)) return "simplex " + std::to_string(x) + ": degeneracy index out of range in face " + std::to_string(i);
                }
                if (sdim(f) != dim[x] - 1) return "simplex " + std::to_string(x) + ": face " + std::to_string(i) + " has the wrong dimension";
                if (birth[f.nd] > birth[x]) return "simplex " + std::to_string(x) + ": born before its face " + std::to_string(i);
            }
        }
        for (size_t x = 0; x < size(); ++x) {
            int n = dim[x]; if (n < 2) continue;
            SSimp X; X.nd = (int32_t)x;
            for (int j = 1; j <= n; ++j) for (int i = 0; i < j; ++i) {
                SSimp a = face(face(X, j), i), b = face(face(X, i), j - 1);
                if (!(a == b)) return "simplex " + std::to_string(x) + ": d_" + std::to_string(i) + " d_" + std::to_string(j) + " != d_" + std::to_string(j - 1) + " d_" + std::to_string(i);
            }
        }
        return "";
    }

    // the normalised chain complex (cells = nondegenerate simplices, same ids)
    FilteredComplex chains() const {
        FilteredComplex C;
        std::map<int32_t, int64_t> acc;
        std::vector<std::pair<int32_t, int64_t>> bd;
        for (size_t x = 0; x < size(); ++x) {
            acc.clear();
            for (size_t i = 0; i < faces[x].size(); ++i) if (!faces[x][i].degenerate()) acc[faces[x][i].nd] += (i & 1) ? -1 : 1;
            bd.clear(); for (auto& kv : acc) if (kv.second) bd.push_back(kv);
            C.add_cell(dim[x], birth[x], bd);
        }
        return C;
    }
};

// ---- constructors --------------------------------------------------------------------------------------

// The nerve of a finite category, up to dimension maxdim. Morphisms m = 0..M-1 with source, target, an
// optional weight and an identity flag (exactly one identity per object); compose(g, f) = g o f for
// composable f: a -> b, g: b -> c. A k-simplex is a chain (f_1, ..., f_k) of composable morphisms
// (f_i: a_{i-1} -> a_i); it is nondegenerate when no f_i is an identity. Faces: d_0 drops f_1, d_k drops f_k,
// d_i composes f_{i+1} o f_i; a chain with identities at positions q_1 < ... < q_m is
// s_{q_m - 1} ... s_{q_1 - 1} of the chain without them. Birth of a chain: the maximum weight of the
// composites f_j o ... o f_i over all runs i <= j (so faces are born no later, for any weights); vertices
// (objects) are born at the least weight of their identity, or at object_birth if given.
struct Category {
    int objects = 0;
    std::vector<int> src, tgt;
    std::vector<double> weight;
    std::vector<int> identity;                    // identity morphism of each object
    std::function<int(int, int)> compose;         // (g, f) -> g o f
    std::vector<double> object_birth;             // optional
};

inline SimplicialSet nerve(const Category& Cat, int maxdim) {
    int M = (int)Cat.src.size();
    std::vector<char> is_id(M, 0); for (int o = 0; o < Cat.objects; ++o) is_id.at(Cat.identity.at(o)) = 1;
    auto w = [&](int m) { return Cat.weight.empty() ? 0.0 : Cat.weight[m]; };
    SimplicialSet S;
    std::map<std::vector<int>, int32_t> id;       // nondegenerate chains (vertices: {-1 - object})
    for (int o = 0; o < Cat.objects; ++o) {
        double b = Cat.object_birth.empty() ? w(Cat.identity[o]) : Cat.object_birth[o];
        id[{-1 - o}] = S.add(0, b, {}, "o" + std::to_string(o));
    }
    // a chain (possibly containing identities) as a simplex of the nerve
    std::function<SSimp(const std::vector<int>&, int, int)> as_simplex = [&](const std::vector<int>& ch, int a0, int dummy) -> SSimp {
        (void)dummy;
        if (ch.empty()) { SSimp v; v.nd = id.at({-1 - a0}); return v; }
        std::vector<int> red, qs;
        for (size_t i = 0; i < ch.size(); ++i) { if (is_id[ch[i]]) qs.push_back((int)i + 1); else red.push_back(ch[i]); }
        SSimp base;
        if (red.empty()) base.nd = id.at({-1 - a0});
        else base.nd = id.at(red);
        std::vector<int> ks; for (size_t r = qs.size(); r-- > 0;) ks.push_back(qs[r] - 1);   // s_{q_m-1} ... s_{q_1-1}
        return make_degenerate(base.nd, ks);
    };
    std::vector<std::vector<int>> level;          // nondegenerate chains of the current length
    for (int m = 0; m < M; ++m) if (!is_id[m]) level.push_back({m});
    for (int k = 1; k <= maxdim && !level.empty(); ++k) {
        for (auto& ch : level) {
            // birth: max weight of the composites of all runs
            double b = -std::numeric_limits<double>::infinity();
            for (size_t i = 0; i < ch.size(); ++i) { int c = ch[i]; b = std::max(b, w(c)); for (size_t j = i + 1; j < ch.size(); ++j) { c = Cat.compose(ch[j], c); b = std::max(b, w(c)); } }
            std::vector<SSimp> f(k + 1);
            int a0 = Cat.src[ch[0]];
            for (int i = 0; i <= k; ++i) {
                std::vector<int> fc; int start = a0;
                if (i == 0) { fc.assign(ch.begin() + 1, ch.end()); start = Cat.tgt[ch[0]]; }
                else if (i == k) fc.assign(ch.begin(), ch.end() - 1);
                else { fc.assign(ch.begin(), ch.begin() + (i - 1)); fc.push_back(Cat.compose(ch[i], ch[i - 1])); fc.insert(fc.end(), ch.begin() + i + 1, ch.end()); }
                f[i] = as_simplex(fc, start, 0);
            }
            id[ch] = S.add(k, b, f);
        }
        if (k == maxdim) break;
        std::vector<std::vector<int>> next;
        for (auto& ch : level) for (int m = 0; m < M; ++m) if (!is_id[m] && Cat.src[m] == Cat.tgt[ch.back()]) { auto c2 = ch; c2.push_back(m); next.push_back(c2); }
        level.swap(next);
    }
    return S;
}

// the normalised bar construction BG of a finite group (the nerve of G as a one-object category), given
// by its multiplication table mult[g][h] = g h and identity e, up to dimension maxdim
inline SimplicialSet bar_construction(const std::vector<std::vector<int>>& mult, int e, int maxdim) {
    Category C; int n = (int)mult.size();
    C.objects = 1; C.src.assign(n, 0); C.tgt.assign(n, 0); C.identity = {e};
    C.compose = [&mult](int g, int f) { return mult[f][g]; };      // f then g: the product f g
    return nerve(C, maxdim);
}

// X / A for a subcomplex A (nondegenerate ids, closed under faces): A becomes one vertex, born at the least
// birth in A. Returns the quotient; new ids: the basepoint first, then the simplices outside A in order.
inline SimplicialSet quotient(const SimplicialSet& X, const std::vector<int32_t>& A_ids) {
    std::vector<char> inA(X.size(), 0); for (int32_t a : A_ids) inA.at(a) = 1;
    if (A_ids.empty()) throw std::invalid_argument("quotient: A is empty");
    for (int32_t a : A_ids) for (auto& f : X.faces[a]) if (!inA[f.nd]) throw std::invalid_argument("quotient: A is not closed under faces");
    SimplicialSet Q;
    double b0 = std::numeric_limits<double>::infinity(); for (int32_t a : A_ids) b0 = std::min(b0, X.birth[a]);
    int32_t star = Q.add(0, b0, {}, "*");
    std::vector<int32_t> nid(X.size(), -1);
    for (size_t x = 0; x < X.size(); ++x) {
        if (inA[x]) continue;
        std::vector<SSimp> f;
        for (auto& g : X.faces[x]) {
            if (inA[g.nd]) { int m = X.sdim(g); std::vector<int> ks; for (int r = m - 1; r >= 0; --r) ks.push_back(r); f.push_back(make_degenerate(star, ks)); }
            else { SSimp h = g; h.nd = nid.at(g.nd); if (h.nd < 0) throw std::invalid_argument("quotient: faces must come before their simplices"); f.push_back(h); }
        }
        nid[x] = Q.add(X.dim[x], X.birth[x], f, X.label[x]);
    }
    return Q;
}

// the cone CX = X * v (cone point last): simplices x, x * v and v; d_i(x * v) = (d_i x) * v for i <= dim x,
// d_{n+1}(x * v) = x; for a vertex x, d_0(x * v) = v. Returns the cone and the ids of X inside it (0..|X|-1).
inline SimplicialSet cone(const SimplicialSet& X) {
    SimplicialSet C = X;
    double b0 = std::numeric_limits<double>::infinity(); for (double b : X.birth) b0 = std::min(b0, b);
    if (X.size() == 0) b0 = 0;
    int32_t v = C.add(0, b0, {}, "v");
    std::vector<int32_t> cid(X.size(), -1);
    for (size_t x = 0; x < X.size(); ++x) {
        int n = X.dim[x];
        std::vector<SSimp> f(n + 2);
        if (n == 0) { f[0].nd = v; f[1].nd = (int32_t)x; }
        else {
            for (int i = 0; i <= n; ++i) {                 // (s_J y) * v = s_J (y * v)
                SSimp g = X.faces[x][i];
                if (cid.at(g.nd) < 0) throw std::invalid_argument("cone: faces must come before their simplices");
                g.nd = cid[g.nd]; f[i] = g;
            }
            f[n + 1].nd = (int32_t)x;
        }
        cid[x] = C.add(n + 1, X.birth[x], f);
    }
    return C;
}
// the suspension CX / X
inline SimplicialSet suspension(const SimplicialSet& X) {
    SimplicialSet C = cone(X);
    std::vector<int32_t> A(X.size()); for (size_t i = 0; i < X.size(); ++i) A[i] = (int32_t)i;
    return quotient(C, A);
}

// the minimal simplicial n-sphere: one vertex and one n-simplex with every face s_{n-2} ... s_0 *
inline SimplicialSet minimal_sphere(int n) {
    SimplicialSet S; int32_t v = S.add(0, 0, {}, "*");
    if (n == 0) { S.add(0, 0, {}, "*'"); return S; }
    std::vector<int> ks; for (int r = n - 2; r >= 0; --r) ks.push_back(r);
    std::vector<SSimp> f(n + 1, make_degenerate(v, ks));
    S.add(n, 0, f, "sigma");
    return S;
}

// text format (perseus sset): one nondegenerate simplex per line, "dim birth face_0 ... face_dim" (vertices:
// "0 birth"); a face is the 0-based line number y of an earlier simplex, or "y@k1.k2...km" for
// s_{k1} s_{k2} ... s_{km} y. Blank lines and lines starting with # do not count.
inline SimplicialSet read_simplicial_set(std::istream& is) {
    SimplicialSet S; std::string line; int ln = 0;
    while (std::getline(is, line)) {
        ++ln;
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line); std::vector<std::string> t; std::string w; while (ss >> w) t.push_back(w);
        if (t.empty()) continue;
        if (t.size() < 2) throw std::runtime_error("sset: line " + std::to_string(ln) + ": expected dim birth faces...");
        int d = std::stoi(t[0]); double b = std::stod(t[1]);
        std::vector<SSimp> f;
        for (size_t i = 2; i < t.size(); ++i) {
            std::string tok = t[i]; std::vector<int> ks; size_t at = tok.find('@');
            if (at != std::string::npos) {
                std::string ds = tok.substr(at + 1); tok = tok.substr(0, at);
                for (char& ch : ds) if (ch == '.' || ch == ',') ch = ' ';
                std::istringstream ks_s(ds); int k; while (ks_s >> k) ks.push_back(k);
            }
            long long y = std::stoll(tok);
            if (y < 0 || y >= (long long)S.size()) throw std::runtime_error("sset: line " + std::to_string(ln) + ": a face must be an earlier simplex");
            f.push_back(make_degenerate((int32_t)y, ks));
        }
        S.add(d, b, f);
    }
    return S;
}

} // namespace perseus
