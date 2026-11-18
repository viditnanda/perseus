// io.hpp -- input formats and output writers.
#pragma once
#include "complex.hpp"
#include "types.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <cstring>
#include <tuple>

namespace perseus {

// A Vietoris-Rips problem: symmetric dissimilarities on n points (lower triangle, row-major:
// d(i,j) for i > j at index i(i-1)/2 + j), maximum homological degree, optional threshold.
struct RipsInput {
    int n = 0;
    std::vector<double> lower;
    int maxdim = 1;
    double threshold = std::numeric_limits<double>::infinity();
    std::vector<double> vertex_birth;   // flag filtrations with vertex births (explicit path); empty = all 0
    // sparse sources (no n x n storage): Euclidean points (row-major, pdim coordinates), or explicit edges
    std::vector<double> pts; int pdim = 0;
    std::vector<std::tuple<int, int, double>> edges; bool from_edges = false;
    double pdist(int i, int j) const {
        double s = 0; const double* a = &pts[(size_t)i * pdim]; const double* b = &pts[(size_t)j * pdim];
        for (int k = 0; k < pdim; ++k) { double t = a[k] - b[k]; s += t * t; }
        return std::sqrt(s);
    }
    double dist(int i, int j) const {
        if (i == j) return 0.0;
        if (lower.empty() && pdim > 0) return pdist(i, j);
        if (i < j) std::swap(i, j); return lower[(size_t)i * (i - 1) / 2 + j];
    }
    bool dense() const { return !lower.empty() || n <= 1; }
    // materialise the lower triangle (for the implicit engine and edge collapse)
    void densify() {
        if (dense()) return;
        lower.assign((size_t)n * (n - 1) / 2, std::numeric_limits<double>::infinity());
        if (pdim > 0) { for (int i = 1; i < n; ++i) for (int j = 0; j < i; ++j) lower[(size_t)i * (i - 1) / 2 + j] = pdist(i, j); }
        else for (auto& e : edges) {
            int a = std::get<0>(e), b = std::get<1>(e); if (a < b) std::swap(a, b);
            double& sl = lower[(size_t)a * (a - 1) / 2 + b]; sl = std::min(sl, std::get<2>(e));
        }
    }
};


// whitespace/comma tokenizer over the whole stream (fast numeric parsing; '#' starts a comment)
struct Tok {
    std::string buf; const char* p; const char* e;
    explicit Tok(std::istream& s) {
        buf.assign(std::istreambuf_iterator<char>(s), std::istreambuf_iterator<char>());
        p = buf.data(); e = p + buf.size();
    }
    void skip() {
        while (p < e) {
            if (*p == '#') { while (p < e && *p != '\n') ++p; }
            else if (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r' || *p == ',') ++p;
            else break;
        }
    }
    bool next(std::string& t) {
        skip(); if (p >= e) return false;
        const char* q = p; while (q < e && !(*q == ' ' || *q == '\t' || *q == '\n' || *q == '\r' || *q == ',')) ++q;
        t.assign(p, q); p = q; return true;
    }
    bool num(double& x) {
        skip(); if (p >= e) return false;
        char* end; x = std::strtod(p, &end);
        if (end == p) throw std::runtime_error("expected a number near '" + std::string(p, std::min<size_t>(20, e - p)) + "'");
        p = end; return true;
    }
    bool inum(long long& x) {
        double d; if (!num(d)) return false;
        x = (long long)d;
        if ((double)x != d) throw std::runtime_error("expected an integer, got " + std::to_string(d));
        return true;
    }
};
static inline std::vector<std::string> split_ws(const std::string& line) {
    std::vector<std::string> out; std::istringstream ss(line); std::string t; while (ss >> t) out.push_back(t); return out;
}


// cubical: dense grid of top-cell values. First line: the sizes n_1 ... n_d. Then n_1*...*n_d values in
// C order (last axis fastest, as numpy writes them); nan or inf marks an absent top cell.
// (Axis order only permutes coordinates, so the array is passed to the builder with its axes reversed.)
// a dense grid of top-cell values, axis 0 fastest (the convention of build_cubical); NaN/inf = absent
struct CubicalGrid {
    std::vector<int64_t> ext; std::vector<bool> periodic; std::vector<double> top;
    const double* view = nullptr;    // optional: values held by the caller instead of top (read by the planar engine)
    const double* data() const { return view ? view : top.data(); }
};
inline CubicalGrid grid_from_array(std::vector<int64_t> shape_c, std::vector<double> vals_c, std::vector<bool> periodic_c) {
    if (periodic_c.empty()) periodic_c.assign(shape_c.size(), false);
    if (periodic_c.size() != shape_c.size()) throw std::runtime_error("cubical: --periodic needs one flag per axis");
    for (auto& v : vals_c) if (std::isnan(v)) v = std::numeric_limits<double>::infinity();
    std::reverse(shape_c.begin(), shape_c.end()); std::reverse(periodic_c.begin(), periodic_c.end());
    return {shape_c, periodic_c, std::move(vals_c)};
}
struct Input {
    bool is_rips = false;
    FilteredComplex complex;
    RipsInput rips;
    bool is_grid = false;                // cubical grid for the implicit cubical engine
    CubicalGrid grid;
};
inline FilteredComplex cubical_from_grid(const CubicalGrid& g) { return build_cubical(g.ext, g.periodic, g.top); }
inline FilteredComplex cubical_from_array(std::vector<int64_t> shape_c, std::vector<double> vals_c, std::vector<bool> periodic_c) {
    return cubical_from_grid(grid_from_array(shape_c, vals_c, periodic_c));
}
inline CubicalGrid read_cubical_grid(std::istream& is, const std::vector<bool>& periodic) {
    std::string line, all; std::vector<int64_t> shape;
    while (std::getline(is, line)) {
        if (line.empty() || line[0] == '#') continue;
        for (char& ch : line) if (ch == ',') ch = ' ';
        for (auto& t : split_ws(line)) { long long v = std::stoll(t); if (v <= 0) throw std::runtime_error("cubical: sizes must be positive"); shape.push_back(v); }
        break;
    }
    if (shape.empty()) throw std::runtime_error("cubical: missing sizes line");
    Tok T(is); int64_t tot = 1; for (auto x : shape) tot *= x;
    std::vector<double> vals(tot);
    for (int64_t i = 0; i < tot; ++i) if (!T.num(vals[i])) throw std::runtime_error("cubical: expected " + std::to_string(tot) + " values");
    return grid_from_array(shape, vals, periodic);
}
inline FilteredComplex read_cubical_text(std::istream& is, const std::vector<bool>& periodic) { return cubical_from_grid(read_cubical_grid(is, periodic)); }
// npy: a numpy array file (any real or integer dtype, C or Fortran order) of top-cell values
struct NpyArray { std::vector<int64_t> shape; std::vector<double> values; bool fortran = false; };
inline NpyArray read_npy_array(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + path);
    char magic[6]; f.read(magic, 6);
    if (std::string(magic, 6) != "\x93NUMPY") throw std::runtime_error("npy: not a numpy file");
    unsigned char ver[2]; f.read((char*)ver, 2);
    uint32_t hlen = 0;
    if (ver[0] == 1) { unsigned char b[2]; f.read((char*)b, 2); hlen = b[0] | (b[1] << 8); }
    else { unsigned char b[4]; f.read((char*)b, 4); hlen = b[0] | (b[1] << 8) | (b[2] << 16) | ((uint32_t)b[3] << 24); }
    std::string h(hlen, ' '); f.read(&h[0], hlen);
    auto field = [&](const std::string& k) { size_t p = h.find("'" + k + "'"); if (p == std::string::npos) throw std::runtime_error("npy: header lacks " + k); return h.substr(h.find(':', p) + 1); };
    std::string descr = field("descr"); descr = descr.substr(descr.find('\'') + 1); descr = descr.substr(0, descr.find('\''));
    bool fortran = field("fortran_order").find("True") < field("fortran_order").find(',');
    std::string sh = field("shape"); sh = sh.substr(sh.find('(') + 1); sh = sh.substr(0, sh.find(')'));
    std::vector<int64_t> shape; for (char& ch : sh) if (ch == ',') ch = ' '; for (auto& t : split_ws(sh)) shape.push_back(std::stoll(t));
    if (shape.empty()) throw std::runtime_error("npy: zero-dimensional array");
    int64_t tot = 1; for (auto x : shape) tot *= x;
    char endian = descr[0], kind = descr[1]; int size = std::stoi(descr.substr(2));
    if (endian == '>') throw std::runtime_error("npy: big-endian arrays are not supported");
    std::vector<char> raw((size_t)tot * size); f.read(raw.data(), raw.size());
    std::vector<double> v(tot);
    for (int64_t i = 0; i < tot; ++i) {
        const char* q = raw.data() + (size_t)i * size;
        if (kind == 'f' && size == 8) { double x; memcpy(&x, q, 8); v[i] = x; }
        else if (kind == 'f' && size == 4) { float x; memcpy(&x, q, 4); v[i] = x; }
        else if (kind == 'i' && size == 8) { int64_t x; memcpy(&x, q, 8); v[i] = (double)x; }
        else if (kind == 'i' && size == 4) { int32_t x; memcpy(&x, q, 4); v[i] = x; }
        else if (kind == 'i' && size == 2) { int16_t x; memcpy(&x, q, 2); v[i] = x; }
        else if ((kind == 'i' || kind == 'u' || kind == 'b') && size == 1) v[i] = kind == 'i' ? (double)(int8_t)*q : (double)(uint8_t)*q;
        else if (kind == 'u' && size == 2) { uint16_t x; memcpy(&x, q, 2); v[i] = x; }
        else if (kind == 'u' && size == 4) { uint32_t x; memcpy(&x, q, 4); v[i] = x; }
        else if (kind == 'u' && size == 8) { uint64_t x; memcpy(&x, q, 8); v[i] = (double)x; }
        else throw std::runtime_error("npy: unsupported dtype " + descr);
    }
    NpyArray a; a.shape = shape; a.values = std::move(v); a.fortran = fortran;
    return a;
}
inline CubicalGrid grid_from_npy(NpyArray a, std::vector<bool> periodic) {
    if (!a.fortran) return grid_from_array(a.shape, std::move(a.values), periodic);
    if (periodic.empty()) periodic.assign(a.shape.size(), false);          // axis 0 fastest already
    if (periodic.size() != a.shape.size()) throw std::runtime_error("npy: --periodic needs one flag per axis");
    for (auto& x : a.values) if (std::isnan(x)) x = std::numeric_limits<double>::infinity();
    return {a.shape, periodic, std::move(a.values)};
}
inline FilteredComplex cubical_from_npy(NpyArray a, std::vector<bool> periodic) { return cubical_from_grid(grid_from_npy(std::move(a), periodic)); }
// cubes: sparse cubical complex, one top cube per line: its d integer anchor coordinates, then its birth.
// A repeated anchor keeps its earliest birth.
inline FilteredComplex read_cubes(std::istream& is) {
    std::map<std::vector<int64_t>, double> top; std::string line; int d = -1;
    while (std::getline(is, line)) {
        if (line.empty() || line[0] == '#') continue;
        for (char& ch : line) if (ch == ',') ch = ' ';
        auto t = split_ws(line); if (t.empty()) continue;
        if (d < 0) d = (int)t.size() - 1;
        if ((int)t.size() != d + 1 || d < 1) throw std::runtime_error("cubes: every line needs the same number of coordinates, then a birth");
        std::vector<int64_t> a(d); for (int i = 0; i < d; ++i) a[i] = std::stoll(t[i]);
        double b = std::stod(t[d]); if (!std::isfinite(b)) continue;
        auto it = top.find(a); if (it == top.end()) top.emplace(a, b); else it->second = std::min(it->second, b);
    }
    std::vector<std::vector<int64_t>> an; std::vector<double> bs;
    for (auto& kv : top) { an.push_back(kv.first); bs.push_back(kv.second); }
    if (an.empty()) throw std::runtime_error("cubes: no cubes");
    return build_cubical_sparse(d, an, bs);
}
// flag: a filtered graph whose clique (flag) complex is the filtration. A line with two tokens "v birth"
// gives a vertex its birth (default 0); a line with three tokens "u v birth" is an edge. Labels are
// arbitrary tokens. An edge may not be born before its endpoints; a repeated edge keeps its earliest birth.
inline RipsInput read_flag(std::istream& is) {
    std::map<std::string, int> lab; std::vector<double> vb;
    std::vector<std::tuple<int, int, double>> ed; std::string line; int ln = 0;
    auto id = [&](const std::string& s) { auto it = lab.find(s); if (it == lab.end()) { it = lab.emplace(s, (int)lab.size()).first; vb.push_back(0.0); } return it->second; };
    while (std::getline(is, line)) {
        ++ln;
        if (line.empty() || line[0] == '#') continue;
        for (char& ch : line) if (ch == ',') ch = ' ';
        auto t = split_ws(line); if (t.empty()) continue;
        if (t.size() == 2) { int v = id(t[0]); vb[v] = std::stod(t[1]); }
        else if (t.size() == 3) { int a = id(t[0]), b = id(t[1]); if (a == b) throw std::runtime_error("flag: loop at line " + std::to_string(ln)); ed.emplace_back(a, b, std::stod(t[2])); }
        else throw std::runtime_error("flag: line " + std::to_string(ln) + ": expected 'v birth' or 'u v birth'");
    }
    RipsInput R; R.n = (int)lab.size(); R.from_edges = true;
    for (auto& e : ed) {
        int a = std::get<0>(e), b = std::get<1>(e); double w = std::get<2>(e);
        if (w < std::max(vb[a], vb[b])) throw std::runtime_error("flag: an edge is born before one of its endpoints");
    }
    R.edges = std::move(ed);
    bool any = false; for (double x : vb) if (x != 0) any = true;
    if (any) R.vertex_birth = vb;
    return R;
}

// ---- point clouds, distances, simplices, cells -------------------------------------------
// points: one point per line (whitespace or comma separated); Euclidean distance.
inline RipsInput read_points(std::istream& is) {
    RipsInput R; Tok T(is); std::string line;
    // dimension from the first data line
    const char* q = T.p;
    while (q < T.e && (*q == '#' || *q == '\n' || *q == '\r')) { if (*q == '#') while (q < T.e && *q != '\n') ++q; else ++q; }
    const char* eol = q; while (eol < T.e && *eol != '\n') ++eol;
    { std::string first(q, eol); for (char& c : first) if (c == ',' || c == ';') c = ' '; R.pdim = (int)split_ws(first).size(); }
    if (R.pdim == 0) return R;
    for (char* c = (char*)T.p; c < T.e; ++c) if (*c == ';') *c = ' ';
    double x; while (T.num(x)) R.pts.push_back(x);
    if (R.pts.size() % R.pdim) throw std::runtime_error("points: inconsistent dimension");
    R.n = (int)(R.pts.size() / R.pdim);
    return R;
}
// lower-distance (Ripser's format): entries of the strictly lower triangle, row by row
inline RipsInput read_lower(std::istream& is) {
    std::vector<double> v; std::string t; Tok T(is);
    while (T.next(t)) { for (char& c : t) if (c == ',') c = ' '; auto u = split_ws(t); for (auto& x : u) v.push_back(std::stod(x)); }
    RipsInput R; R.n = (int)std::llround((1 + std::sqrt(1 + 8.0 * v.size())) / 2);
    if ((size_t)R.n * (R.n - 1) / 2 != v.size()) throw std::runtime_error("lower-distance: not a triangular number of entries");
    R.lower = v; return R;
}
// full distance matrix (whitespace/comma separated rows)
inline RipsInput read_matrix(std::istream& is) {
    std::vector<std::vector<double>> M; std::string line;
    while (std::getline(is, line)) {
        if (line.empty() || line[0] == '#') continue;
        for (char& c : line) if (c == ',' || c == ';') c = ' ';
        auto t = split_ws(line); if (t.empty()) continue;
        std::vector<double> r; for (auto& x : t) r.push_back(std::stod(x)); M.push_back(r);
    }
    RipsInput R; R.n = (int)M.size(); R.lower.resize((size_t)R.n * (R.n - 1) / 2);
    for (int i = 1; i < R.n; ++i) for (int j = 0; j < i; ++j) R.lower[(size_t)i * (i - 1) / 2 + j] = M[i][j];
    return R;
}
// binary lower triangle of float32 (or float64 if `f64`)
inline RipsInput read_binary(const std::string& path, bool f64) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f) throw std::runtime_error("cannot open " + path);
    size_t bytes = (size_t)f.tellg(); f.seekg(0);
    size_t m = bytes / (f64 ? 8 : 4);
    RipsInput R; R.n = (int)std::llround((1 + std::sqrt(1 + 8.0 * m)) / 2); R.lower.resize(m);
    if (f64) f.read((char*)R.lower.data(), m * 8);
    else { std::vector<float> v(m); f.read((char*)v.data(), m * 4); for (size_t i = 0; i < m; ++i) R.lower[i] = v[i]; }
    return R;
}
// simplices: one simplex per line, vertex labels then the birth as the last token; faces are added
// with the minimum birth of their cofaces
inline FilteredComplex read_simplices(std::istream& is) {
    SimplicialBuilder B; std::string line; std::map<std::string, int> lab;
    while (std::getline(is, line)) {
        if (line.empty() || line[0] == '#') continue;
        for (char& c : line) if (c == ',') c = ' ';
        auto t = split_ws(line); if (t.size() < 2) continue;
        std::vector<int> s;
        for (size_t i = 0; i + 1 < t.size(); ++i) { auto it = lab.find(t[i]); if (it == lab.end()) it = lab.emplace(t[i], (int)lab.size()).first; s.push_back(it->second); }
        B.add(s, std::stod(t.back()));
    }
    return B.build();
}
// cells: general complex, one cell per line: dim birth [face coef]... ; faces are 0-based line
// numbers of earlier cells (blank and comment lines do not count). Coefficients are integers.
inline FilteredComplex read_cells(std::istream& is) {
    FilteredComplex C; std::string line; int ln = 0;
    while (std::getline(is, line)) {
        ++ln;
        if (line.empty() || line[0] == '#') continue;
        auto t = split_ws(line); if (t.empty()) continue;
        if (t.size() < 2 || t.size() % 2) throw std::runtime_error("cells: line " + std::to_string(ln) + ": expected dim birth [face coef]...");
        std::vector<std::pair<int32_t, int64_t>> bd;
        for (size_t i = 2; i < t.size(); i += 2) {
            long long f = std::stoll(t[i]);
            if (f < 0 || f >= (long long)C.size()) throw std::runtime_error("cells: line " + std::to_string(ln) + ": face must be an earlier cell");
            bd.push_back({(int32_t)f, std::stoll(t[i + 1])});
        }
        C.add_cell(std::stoi(t[0]), std::stod(t[1]), bd);
    }
    return C;
}

// ---- outputs ---------------------------------------------------------------------------
inline std::string fmt_num(double x) {
    if (std::isinf(x)) return "inf";
    char b[64];                                        // shortest decimal that reads back as exactly x
    for (int prec = 6; prec <= 17; ++prec) { snprintf(b, sizeof b, "%.*g", prec, x); if (std::strtod(b, nullptr) == x) break; }
    return b;
}
// files: prefix_k.txt ("birth death", death -1 if infinite), prefix_betti.txt
inline void write_files(const std::string& prefix, const Barcode& bc) {
    for (size_t k = 0; k < bc.size(); ++k) {
        std::ofstream f(prefix + "_" + std::to_string(k) + ".txt");
        for (auto& b : bc[k]) f << fmt_num(b.birth) << " " << (std::isinf(b.death) ? std::string("-1") : fmt_num(b.death)) << "\n";
    }
    // Betti numbers at every time they may change: sweep over sorted births and deaths per degree
    std::set<double> times;
    for (auto& v : bc) for (auto& b : v) { times.insert(b.birth); if (!std::isinf(b.death)) times.insert(b.death); }
    std::vector<std::vector<double>> B(bc.size()), Dd(bc.size());
    for (size_t k = 0; k < bc.size(); ++k) {
        for (auto& b : bc[k]) { B[k].push_back(b.birth); Dd[k].push_back(b.death); }
        std::sort(B[k].begin(), B[k].end()); std::sort(Dd[k].begin(), Dd[k].end());
    }
    std::vector<size_t> ib(bc.size(), 0), id(bc.size(), 0);
    std::ofstream f(prefix + "_betti.txt");
    for (double t : times) {
        f << "\n " << fmt_num(t) << " ";
        for (size_t k = 0; k < bc.size(); ++k) {
            while (ib[k] < B[k].size() && B[k][ib[k]] <= t) ++ib[k];     // born by t
            while (id[k] < Dd[k].size() && Dd[k][id[k]] <= t) ++id[k];   // dead by t
            f << (long)(ib[k] - id[k]) << " ";
        }
    }
}
inline void print_ripser(FILE* out, const Barcode& bc) {
    for (size_t k = 0; k < bc.size(); ++k) {
        fprintf(out, "persistence intervals in dim %zu:\n", k);
        for (auto& b : bc[k]) {
            if (std::isinf(b.death)) fprintf(out, " [%s, )\n", fmt_num(b.birth).c_str());
            else fprintf(out, " [%s,%s)\n", fmt_num(b.birth).c_str(), fmt_num(b.death).c_str());
        }
    }
}

} // namespace perseus
