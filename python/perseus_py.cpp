// Python bindings for Perseus 6 (pybind11).
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "perseus/io.hpp"
#include "perseus/explicit_engine.hpp"
#include "perseus/rips_engine.hpp"
#include "perseus/rips_sparse.hpp"
#include "perseus/cubical_engine.hpp"
#include "perseus/certify.hpp"
#include "perseus/collapse.hpp"
#include <thread>
#include <cmath>

namespace py = pybind11;
using namespace perseus;

static py::list bars_to_py(const Barcode& bc) {
    py::list out;
    for (auto& v : bc) {
        py::array_t<double> a({(py::ssize_t)v.size(), (py::ssize_t)2});
        auto m = a.mutable_unchecked<2>();
        for (size_t i = 0; i < v.size(); ++i) { m(i, 0) = v[i].birth; m(i, 1) = v[i].death; }
        out.append(a);
    }
    return out;
}
static py::object bigint_to_py(const BigInt& x) { return py::int_(py::str(x.str())); }

struct SolveOpts {
    std::string coeff = "Q", mode = "auto";
    bool classify = true;
    int threads = 0;
    std::vector<uint64_t> fp;
    uint64_t seed = 0;
    double target = 64;
};
static SolveOpts make_opts(py::object coeff, bool classify, int threads, std::vector<uint64_t> fp, std::string mode, py::object seed, double target) {
    SolveOpts o; o.coeff = py::str(coeff).cast<std::string>(); o.classify = classify; o.threads = threads; o.fp = std::move(fp);
    if (mode != "auto" && mode != "fast" && mode != "proved") throw std::invalid_argument("mode must be 'auto', 'fast' or 'proved'");
    o.mode = mode; o.seed = seed.is_none() ? 0 : seed.cast<uint64_t>(); o.target = target;
    return o;
}

template <class B>
static py::dict solve(const B& be, const SolveOpts& o) {
    py::dict res;
    int threads = o.threads > 0 ? o.threads : (int)std::max(1u, std::thread::hardware_concurrency());
    bool qmode = (o.coeff == "Q" || o.coeff == "q" || o.coeff == "0");
    if (qmode) {
        CertifyOptions co; co.threads = threads; co.classify = o.classify; co.mode = o.mode; co.seed = o.seed; co.target_log2 = -o.target;
        std::pair<BRun<typename B::Clear>, Certificate> r;
        { py::gil_scoped_release nogil; r = certify(be, co); }
        const Certificate& c = r.second;
        res["coefficients"] = "Q";
        res["barcode"] = bars_to_py(r.first.bars);
        py::dict cert;
        cert["mode"] = c.mode;
        cert["barcode_certified"] = c.barcode_ok;
        if (!c.note.empty()) cert["note"] = c.note;
        if (c.mode == "fast") {
            cert["log2_error_bound"] = std::isinf(c.log2_error) ? py::object(py::none()) : py::object(py::float_(c.log2_error));
            cert["seed"] = py::int_(c.seed); cert["lanes"] = c.lanes; cert["retries"] = c.retries;
            py::list mods; for (auto m : c.moduli) mods.append(py::int_(m)); cert["moduli"] = mods;
        }
        if (c.mode == "proved") cert["arithmetic"] = c.arithmetic;
        cert["bad_primes_status"] = !c.candidates_complete() ? "incomplete" : (c.classified ? "complete" : "candidates only");
        cert["complete"] = c.complete();
        py::list vals;
        py::object Fraction = py::module_::import("fractions").attr("Fraction");
        for (auto& v : c.values) {
            py::dict d; d["degree"] = v.k; d["birth"] = v.birth; d["death"] = v.death;
            d["value"] = v.ok ? Fraction(bigint_to_py(v.num), bigint_to_py(v.den)) : py::none();
            vals.append(d);
        }
        cert["nonunit_pivots"] = vals;
        py::list cand; for (auto& kv : c.candidates) cand.append(py::int_(kv.first));
        for (auto& kv : c.big_candidates) cand.append(bigint_to_py(kv.first));
        cert["candidates"] = cand;
        py::dict bad;
        for (auto& b : c.bad) {
            py::dict d; d["first_degree"] = b.k0; d["only_over_Fp"] = bars_to_py(b.only_p); d["only_over_Q"] = bars_to_py(b.only_q);
            d["every_prime_factor_of"] = b.group;
            d["probable_prime"] = b.is_big && !b.group && b.big.bits() > 64;
            bad[b.is_big ? bigint_to_py(b.big) : py::object(py::int_(b.p))] = d;
        }
        cert["bad_primes"] = bad;
        py::list unres; for (auto& u : c.unresolved) unres.append(bigint_to_py(u));
        cert["unresolved_factors"] = unres;
        cert["unreconstructed_pivots"] = c.unreconstructed;
        py::dict secs; secs["Q"] = c.t_q; secs["factoring"] = c.t_recon; secs["classification"] = c.t_classify;
        cert["seconds"] = secs;
        res["certificate"] = cert;
    } else {
        uint64_t p = std::stoull(o.coeff);
        if (!is_prime64(p)) throw std::invalid_argument("coeff must be 'Q' or a prime");
        Barcode b; { py::gil_scoped_release nogil; b = run_field(be, p, 0, be.maxdim(), nullptr).bars; }
        res["coefficients"] = py::int_(p);
        res["barcode"] = bars_to_py(b);
    }
    if (!o.fp.empty()) {
        py::dict fb;
        for (uint64_t p : o.fp) {
            if (!is_prime64(p)) throw std::invalid_argument("fp entries must be primes");
            Barcode b; { py::gil_scoped_release nogil; b = run_field(be, p, 0, be.maxdim(), nullptr).bars; }
            fb[py::int_(p)] = bars_to_py(b);
        }
        res["fields"] = fb;
    }
    return res;
}

// the same route choice as the command line: the sparse flag complex for vertex births, for
// engine="sparse", or automatically for n >= 2048 with average degree <= 64; otherwise the implicit
// Rips engine, edge-collapsed for maxdim >= 2
static py::dict solve_rips(RipsInput& R, int maxdim, const std::string& engine, const std::string& collapse, const SolveOpts& o) {
    if (engine != "auto" && engine != "implicit" && engine != "sparse") throw std::invalid_argument("engine must be 'auto', 'implicit' or 'sparse'");
    if (!std::isfinite(R.threshold)) R.threshold = enclosing_radius_any(R);
    std::vector<WEdge> E;
    bool sparse = engine == "sparse" || !R.vertex_birth.empty();
    if (!sparse && engine == "auto" && R.n >= 2048) { E = radius_edges(R, R.threshold); sparse = (double)E.size() <= 32.0 * R.n; }
    if (sparse) {
        if (E.empty()) E = radius_edges(R, R.threshold);
        FilteredComplex C = flag_complex(R.n, E, R.vertex_birth, maxdim + 1);
        Prepared P(C); C = FilteredComplex();
        ExplicitBackend be(P, maxdim, auto_direction(P, maxdim));
        return solve(be, o);
    }
    R.densify();
    if (collapse == "on" || (collapse == "auto" && maxdim >= 2 && R.n >= 32)) {
        CollapseStats cs; RipsInput Cc = edge_collapse(R, 32, cs, 0.02, collapse == "auto" ? 0.10 : 0.0);
        if (cs.used) R.lower.swap(Cc.lower);
    }
    RipsBackend be(R, maxdim);
    return solve(be, o);
}

static RipsInput rips_from_numpy(py::object points, py::object distances) {
    RipsInput R;
    if (!points.is_none()) {
        auto P = py::array_t<double, py::array::c_style | py::array::forcecast>(points);
        if (P.ndim() != 2) throw std::invalid_argument("points must be a 2d array (n, d)");
        R.n = (int)P.shape(0); R.pdim = (int)P.shape(1);
        R.pts.assign(P.data(), P.data() + P.size());
        if (R.pdim == 0) throw std::invalid_argument("points must have at least one coordinate");
    } else if (!distances.is_none()) {
        auto D = py::array_t<double, py::array::c_style | py::array::forcecast>(distances);
        if (D.ndim() == 2) {
            auto a = D.unchecked<2>(); int n = (int)D.shape(0);
            if (D.shape(1) != n) throw std::invalid_argument("distance matrix must be square");
            R.n = n; R.lower.resize((size_t)n * (n - 1) / 2);
            for (int i = 1; i < n; ++i) for (int j = 0; j < i; ++j) R.lower[(size_t)i * (i - 1) / 2 + j] = a(i, j);
        } else if (D.ndim() == 1) {      // condensed lower triangle, row by row (Ripser order)
            auto a = D.unchecked<1>(); size_t m = (size_t)D.shape(0);
            R.n = (int)std::llround((1 + std::sqrt(1 + 8.0 * m)) / 2);
            if ((size_t)R.n * (R.n - 1) / 2 != m) throw std::invalid_argument("condensed distances: length is not triangular");
            R.lower.resize(m); for (size_t i = 0; i < m; ++i) R.lower[i] = a(i);
        } else throw std::invalid_argument("distances must be 1d (lower triangle) or 2d");
    } else throw std::invalid_argument("give points or distances");
    return R;
}

#define COMMON_ARGS py::arg("coeff") = "Q", py::arg("classify") = true, py::arg("threads") = 0, \
    py::arg("fp") = std::vector<uint64_t>(), py::arg("mode") = "auto", py::arg("seed") = py::none(), py::arg("target") = 64.0

PYBIND11_MODULE(_perseus, m) {
    m.doc() = "Perseus 6: persistent homology over Q with certified bad primes, or over any F_p.";

    m.def("rips", [](py::object points, py::object distances, int maxdim, py::object threshold, std::string engine, std::string collapse,
                     bool explicit_engine, py::object coeff, bool classify, int threads, std::vector<uint64_t> fp, std::string mode, py::object seed, double target) {
        SolveOpts o = make_opts(coeff, classify, threads, fp, mode, seed, target);
        RipsInput R = rips_from_numpy(points, distances);
        R.threshold = threshold.is_none() ? std::numeric_limits<double>::infinity() : threshold.cast<double>();
        return solve_rips(R, maxdim, explicit_engine ? "sparse" : engine, collapse, o);
    }, py::arg("points") = py::none(), py::arg("distances") = py::none(), py::arg("maxdim") = 1, py::arg("threshold") = py::none(),
       py::arg("engine") = "auto", py::arg("collapse") = "auto", py::arg("explicit_engine") = false, COMMON_ARGS,
       "Vietoris-Rips persistence of a point cloud (n, d) or a distance matrix (square or condensed).");

    m.def("flag", [](std::vector<std::tuple<int, int, double>> edges, std::vector<double> vertex_births, int maxdim, std::string engine,
                     std::string collapse, py::object coeff, bool classify, int threads, std::vector<uint64_t> fp, std::string mode, py::object seed, double target) {
        SolveOpts o = make_opts(coeff, classify, threads, fp, mode, seed, target);
        int n = (int)vertex_births.size();
        for (auto& e : edges) n = std::max(n, std::max(std::get<0>(e), std::get<1>(e)) + 1);
        std::vector<double> vb(n, 0.0); for (size_t i = 0; i < vertex_births.size(); ++i) vb[i] = vertex_births[i];
        RipsInput R; R.n = n; R.from_edges = true; R.threshold = std::numeric_limits<double>::infinity();
        for (auto& e : edges) {
            int a = std::get<0>(e), b = std::get<1>(e); double w = std::get<2>(e);
            if (a == b || a < 0 || b < 0) throw std::invalid_argument("flag: bad edge");
            if (w < std::max(vb[a], vb[b])) throw std::invalid_argument("flag: an edge is born before one of its endpoints");
        }
        R.edges = std::move(edges);
        bool any = false; for (double x : vb) if (x != 0) any = true;
        if (any) R.vertex_birth = vb;
        return solve_rips(R, maxdim, engine, collapse, o);
    }, py::arg("edges"), py::arg("vertex_births") = std::vector<double>(), py::arg("maxdim") = 1, py::arg("engine") = "auto",
       py::arg("collapse") = "auto", COMMON_ARGS,
       "Clique (flag) filtration of a weighted graph: edges as (u, v, birth); optional vertex births.");

    m.def("cells", [](std::vector<std::tuple<int, double, std::vector<std::pair<int32_t, int64_t>>>> cells, int maxdim, bool zero_length,
                      std::string direction, py::object coeff, bool classify, int threads, std::vector<uint64_t> fp, std::string mode, py::object seed, double target) {
        SolveOpts o = make_opts(coeff, classify, threads, fp, mode, seed, target);
        FilteredComplex C;
        for (auto& t : cells) C.add_cell(std::get<0>(t), std::get<1>(t), std::get<2>(t));
        std::string err = C.validate(false); if (!err.empty()) throw std::invalid_argument(err);
        Prepared P(C); int md_ = maxdim < 0 ? P.D : maxdim;
        int dir = direction == "cohomology" ? 0 : direction == "homology" ? 1 : direction == "forward" ? 2 : auto_direction(P, md_);
        ExplicitBackend be(P, md_, dir); be.keep_zero_length = zero_length;
        return solve(be, o);
    }, py::arg("cells"), py::arg("maxdim") = -1, py::arg("zero_length") = false, py::arg("direction") = "auto", COMMON_ARGS,
       "General filtered cell complex: list of (dim, birth, [(face index, integer coefficient), ...]).");

    m.def("simplices", [](std::vector<std::pair<std::vector<int>, double>> simp, int maxdim, py::object coeff, bool classify, int threads,
                          std::vector<uint64_t> fp, std::string mode, py::object seed, double target) {
        SolveOpts o = make_opts(coeff, classify, threads, fp, mode, seed, target);
        SimplicialBuilder B; for (auto& s : simp) B.add(s.first, s.second);
        FilteredComplex C = B.build();
        Prepared P(C); int md_ = maxdim < 0 ? P.D : maxdim; ExplicitBackend be(P, md_, auto_direction(P, md_));
        return solve(be, o);
    }, py::arg("simplices"), py::arg("maxdim") = -1, COMMON_ARGS,
       "Simplicial complex from (vertices, birth) pairs; faces get the minimum birth of their cofaces.");

    m.def("cubical", [](py::array_t<double, py::array::f_style | py::array::forcecast> top, std::vector<bool> periodic, int maxdim, std::string engine,
                        py::object coeff, bool classify, int threads, std::vector<uint64_t> fp, std::string mode, py::object seed, double target) {
        SolveOpts o = make_opts(coeff, classify, threads, fp, mode, seed, target);
        CubicalGrid g;
        for (py::ssize_t i = 0; i < top.ndim(); ++i) g.ext.push_back(top.shape(i));
        if (periodic.empty()) periodic.assign(g.ext.size(), false);
        if (periodic.size() != g.ext.size()) throw std::invalid_argument("periodic must have one flag per axis");
        g.periodic = periodic;
        g.top.assign(top.data(), top.data() + top.size());       // Fortran order: axis 0 fastest
        for (auto& x : g.top) if (std::isnan(x)) x = std::numeric_limits<double>::infinity();   // NaN = absent cube
        if (engine == "explicit") {
            FilteredComplex C = cubical_from_grid(g);
            Prepared P(C); int md_ = maxdim < 0 ? P.D : maxdim; ExplicitBackend be(P, md_, auto_direction(P, md_));
            return solve(be, o);
        }
        if (engine != "auto" && engine != "implicit") throw std::invalid_argument("engine must be 'auto', 'implicit' or 'explicit'");
        CubicalBackend be(g, maxdim);
        return solve(be, o);
    }, py::arg("top"), py::arg("periodic") = std::vector<bool>(), py::arg("maxdim") = -1, py::arg("engine") = "auto", COMMON_ARGS,
       "Cubical complex from an n-d array of top-cell values (NaN or inf = absent); optional periodic axes.");
}
