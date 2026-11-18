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
#include "perseus/chain.hpp"
#include "perseus/simplicial_set.hpp"
#include "perseus/magnitude.hpp"
#include "perseus/laurent.hpp"
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
    std::string seed;
    double target = 64;
};
static SolveOpts make_opts(py::object coeff, bool classify, int threads, std::vector<uint64_t> fp, std::string mode, py::object seed, double target) {
    SolveOpts o; o.coeff = py::str(coeff).cast<std::string>(); o.classify = classify; o.threads = threads; o.fp = std::move(fp);
    if (mode != "auto" && mode != "fast" && mode != "proved") throw std::invalid_argument("mode must be 'auto', 'fast' or 'proved'");
    o.mode = mode; if (seed.is_none()) o.seed.clear();
    else if (py::isinstance<py::int_>(seed)) o.seed = py::str(py::module_::import("builtins").attr("format")(seed, "x")).cast<std::string>();
    else o.seed = seed.cast<std::string>();
    if (!o.seed.empty()) parse_seed(o.seed); o.target = target;
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
        cert["barcode_certified"] = c.certified();
        if (!c.note.empty()) cert["note"] = c.note;
        if (c.mode == "fast") {
            cert["log2_error_bound"] = std::isinf(c.log2_error) ? py::object(py::none()) : py::object(py::float_(c.log2_error));
            cert["seed"] = c.seed; cert["seed_supplied"] = c.seed_user; cert["lanes"] = c.lanes; cert["retries"] = c.retries;
            py::list mods; for (auto m : c.moduli) mods.append(py::int_(m)); cert["moduli"] = mods;
        }
        if (c.mode == "proved" || c.mode == "unimodular") cert["arithmetic"] = c.arithmetic;
        cert["bad_primes_status"] = c.seed_user && c.barcode_ok ? "conditional on the supplied seed" : !c.candidates_complete() ? "incomplete" : (c.classified ? "complete" : "candidates only");
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


static py::list integral_to_py(const IntegralHomology& H) {
    py::list out;
    for (size_t k = 0; k < H.free.size(); ++k) {
        py::dict d; d["degree"] = k; d["free"] = H.free[k];
        py::list t; for (auto& x : H.torsion[k]) t.append(bigint_to_py(x));
        d["torsion"] = t; out.append(d);
    }
    return out;
}

// a filtered based chain complex over Z (cells, simplicial sets, nerves): exact validation (including
// d o d = 0 over Z), the algebraic Morse reduction of unit pairs with equal births, the certified
// reduction, and on request the integral outputs
static py::dict solve_chain(FilteredComplex C, int maxdim, bool zero_length, const std::string& direction, bool reduce,
                            bool integral, py::object image, const SolveOpts& o) {
    std::string err = validate_chain_complex(C); if (!err.empty()) throw std::invalid_argument(err);
    int md_ = maxdim < 0 ? std::max(0, C.max_dim()) : maxdim;
    ReducedComplex R;
    if (reduce) { py::gil_scoped_release nogil; R = algebraic_morse(C); }
    else { R.C = C; for (size_t i = 0; i < C.size(); ++i) R.orig.push_back((int32_t)i); }
    Prepared P(R.C);
    int dir = direction == "cohomology" ? 0 : direction == "homology" ? 1 : direction == "forward" ? 2 : auto_direction(P, std::min(md_, std::max(P.D, 0)));
    ReducedBackend be(R, P, md_, dir); be.keep_zero_length = zero_length;
    py::dict res = solve(be, o);
    py::dict red; red["cells"] = C.size(); red["after"] = R.C.size(); red["pairs"] = R.pairs; red["skipped_overflow"] = R.skipped_overflow;
    res["reduction"] = red;
    if (integral || !image.is_none()) {
        py::dict I;
        if (integral) {
            IntegralHomology H; { py::gil_scoped_release nogil; H = integral_homology(R.C, md_); }
            I["homology"] = integral_to_py(H);
            bool qmode = (o.coeff == "Q" || o.coeff == "q" || o.coeff == "0");
            if (qmode && res.contains("certificate")) {
                py::dict fb;
                for (auto item : res["certificate"]["bad_primes"].cast<py::dict>()) {
                    if (!py::isinstance<py::int_>(item.first)) continue;
                    py::int_ pi = item.first.cast<py::int_>();
                    if (pi.attr("bit_length")().cast<int>() > 63) continue;
                    uint64_t p = pi.cast<uint64_t>();
                    if (!is_prime64(p)) continue;
                    Barcode b; { py::gil_scoped_release nogil; b = run_field(be, p, 0, be.maxdim(), nullptr).bars; }
                    fb[py::int_(p)] = bars_to_py(b);
                }
                I["bad_prime_barcodes"] = fb;
            }
        }
        if (!image.is_none()) {
            auto st = image.cast<std::pair<double, double>>();
            IntegralHomology H; { py::gil_scoped_release nogil; H = integral_image(R.C, st.first, st.second, md_); }
            py::dict d; d["s"] = st.first; d["t"] = st.second; d["groups"] = integral_to_py(H);
            I["image"] = d;
        }
        res["integral"] = I;
    }
    return res;
}

// faces of a simplicial set from Python: an int y, or (y, [k1, k2, ...]) meaning s_{k1} s_{k2} ... y
static SSimp face_from_py(py::handle f) {
    if (py::isinstance<py::int_>(f)) return make_degenerate(f.cast<int32_t>(), {});
    auto t = f.cast<py::tuple>();
    if (t.size() != 2) throw std::invalid_argument("a face is y or (y, [k1, k2, ...])");
    return make_degenerate(t[0].cast<int32_t>(), t[1].cast<std::vector<int>>());
}
static SimplicialSet sset_from_py(py::list simplices) {
    SimplicialSet S;
    for (auto item : simplices) {
        auto t = item.cast<py::tuple>();
        if (t.size() < 2) throw std::invalid_argument("a simplex is (dim, birth, faces)");
        int d = t[0].cast<int>(); double b = t[1].cast<double>();
        std::vector<SSimp> f;
        if (t.size() > 2) for (auto x : t[2].cast<py::list>()) {
            SSimp y = face_from_py(x);
            if (y.nd < 0 || (size_t)y.nd >= S.size()) throw std::invalid_argument("a face must be an earlier simplex");
            f.push_back(y);
        }
        S.add(d, b, f);
    }
    std::string e = S.validate(); if (!e.empty()) throw std::invalid_argument("invalid simplicial set: " + e);
    return S;
}
static py::dict sset_to_py(const SimplicialSet& S) {
    py::list out;
    for (size_t x = 0; x < S.size(); ++x) {
        py::list f;
        for (auto& y : S.faces[x]) { if (y.J.empty()) f.append(y.nd); else f.append(py::make_tuple(y.nd, y.J)); }
        out.append(py::make_tuple(S.dim[x], S.birth[x], f));
    }
    py::dict d; d["simplices"] = out; return d;
}

// ---- Laurent coefficients ----
static laurent::LPoly lpoly_from_py(py::handle h) {
    if (py::isinstance<py::str>(h)) return laurent::parse_lpoly(h.cast<std::string>());
    if (py::isinstance<py::int_>(h)) return laurent::parse_lpoly(py::str(h).cast<std::string>());
    if (py::isinstance<py::dict>(h)) {
        laurent::LPoly r;
        for (auto kv : h.cast<py::dict>()) {
            int e = kv.first.cast<int>(); laurent::LPoly m = lpoly_from_py(kv.second);
            m.v += e; r = laurent::ladd(r, m);
        }
        return r;
    }
    throw std::invalid_argument("a Laurent coefficient is an int, a string such as '1 - t^-1', or a dict {exponent: coefficient}");
}
static py::list lbars_to_py(const laurent::LBarcode& B) {
    py::list out;
    for (auto& v : B) {
        py::array_t<double> a({(py::ssize_t)v.size(), (py::ssize_t)2});
        auto m = a.mutable_unchecked<2>();
        for (size_t i = 0; i < v.size(); ++i) { m(i, 0) = v[i].first; m(i, 1) = v[i].second; }
        out.append(a);
    }
    return out;
}
static py::list lhom_to_py(const laurent::LHomology& H) {
    py::list out;
    for (size_t k = 0; k < H.free.size(); ++k) {
        py::dict d; d["degree"] = k; d["free"] = H.free[k];
        py::list t; for (auto& f : H.torsion[k]) t.append(laurent::to_string(f));
        d["torsion"] = t; out.append(d);
    }
    return out;
}

// the same route choice as the command line: the sparse flag complex for vertex births, for
// engine="sparse", or automatically when prefer_sparse says the flag complex is small; otherwise the implicit
// Rips engine, edge-collapsed for maxdim >= 2
static py::dict solve_rips(RipsInput& R, int maxdim, const std::string& engine, const std::string& collapse, const SolveOpts& o) {
    if (engine != "auto" && engine != "implicit" && engine != "sparse") throw std::invalid_argument("engine must be 'auto', 'implicit' or 'sparse'");
    if (!std::isfinite(R.threshold)) R.threshold = enclosing_radius_any(R);
    std::vector<WEdge> E;
    bool sparse = engine == "sparse" || !R.vertex_birth.empty();
    if (!sparse && engine == "auto") sparse = prefer_sparse(R, E, maxdim + 1);
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
                      std::string direction, bool reduce, bool integral, py::object image, py::object coeff, bool classify, int threads, std::vector<uint64_t> fp, std::string mode, py::object seed, double target) {
        SolveOpts o = make_opts(coeff, classify, threads, fp, mode, seed, target);
        FilteredComplex C;
        for (auto& t : cells) {
            for (auto& fc : std::get<2>(t)) if (fc.first < 0 || (size_t)fc.first >= C.size()) throw std::invalid_argument("cells: a face must be an earlier cell");
            C.add_cell(std::get<0>(t), std::get<1>(t), std::get<2>(t));
        }
        return solve_chain(std::move(C), maxdim, zero_length, direction, reduce, integral, image, o);
    }, py::arg("cells"), py::arg("maxdim") = -1, py::arg("zero_length") = false, py::arg("direction") = "auto",
       py::arg("reduce") = true, py::arg("integral") = false, py::arg("image") = py::none(), COMMON_ARGS,
       "Filtered based chain complex over Z: list of (dim, birth, [(face index, integer coefficient), ...]).\n"
       "d o d = 0 is checked exactly; unit pairs with equal births are cancelled first (reduce=False to skip).\n"
       "integral=True adds the integral homology of the whole complex and the F_p barcode of each bad prime;\n"
       "image=(s, t) adds the image of H_*(K_s; Z) -> H_*(K_t; Z).");

    m.def("simplicial_set", [](py::list simplices, int maxdim, bool zero_length, bool reduce, bool integral, py::object image,
                               py::object coeff, bool classify, int threads, std::vector<uint64_t> fp, std::string mode, py::object seed, double target) {
        SolveOpts o = make_opts(coeff, classify, threads, fp, mode, seed, target);
        SimplicialSet S = sset_from_py(simplices);
        return solve_chain(S.chains(), maxdim, zero_length, "auto", reduce, integral, image, o);
    }, py::arg("simplices"), py::arg("maxdim") = -1, py::arg("zero_length") = false, py::arg("reduce") = true,
       py::arg("integral") = false, py::arg("image") = py::none(), COMMON_ARGS,
       "Finite filtered simplicial set: list of nondegenerate simplices (dim, birth, [d_0 x, ..., d_n x]); a face\n"
       "is an earlier simplex y or (y, [k1, ..., km]) for s_{k1} ... s_{km} y. Normalised chains.");

    m.def("nerve", [](std::vector<int> src, std::vector<int> tgt, std::vector<int> identities, py::object composition,
                      std::vector<double> weights, int maxdim, bool reduce, bool integral, py::object image, bool return_sset,
                      py::object coeff, bool classify, int threads, std::vector<uint64_t> fp, std::string mode, py::object seed, double target) {
        SolveOpts o = make_opts(coeff, classify, threads, fp, mode, seed, target);
        if (src.size() != tgt.size()) throw std::invalid_argument("src and tgt must have the same length");
        int M = (int)src.size();
        Category Cat; Cat.objects = (int)identities.size(); Cat.src = src; Cat.tgt = tgt; Cat.identity = identities; Cat.weight = weights;
        if (!weights.empty() && (int)weights.size() != M) throw std::invalid_argument("one weight per morphism");
        std::vector<std::vector<int>> table(M, std::vector<int>(M, -1));
        if (py::isinstance<py::dict>(composition)) {
            for (auto kv : composition.cast<py::dict>()) { auto gf = kv.first.cast<std::pair<int, int>>(); table.at(gf.first).at(gf.second) = kv.second.cast<int>(); }
        } else {
            auto t = composition.cast<std::vector<std::vector<int>>>();
            if ((int)t.size() != M) throw std::invalid_argument("composition must be M x M");
            for (int g = 0; g < M; ++g) { if ((int)t[g].size() != M) throw std::invalid_argument("composition must be M x M"); table[g] = t[g]; }
        }
        for (int o2 = 0; o2 < Cat.objects; ++o2) {
            int e = identities[o2]; if (e < 0 || e >= M || src[e] != o2 || tgt[e] != o2) throw std::invalid_argument("bad identity");
            for (int f = 0; f < M; ++f) { if (tgt[f] == o2) table[e][f] = f; if (src[f] == o2) table[f][e] = f; }
        }
        for (int g = 0; g < M; ++g) for (int f = 0; f < M; ++f) if (tgt[f] == src[g]) {
            int h = table[g][f];
            if (h < 0 || h >= M || src[h] != src[f] || tgt[h] != tgt[g]) throw std::invalid_argument("composition (" + std::to_string(g) + ", " + std::to_string(f) + ") missing or with the wrong ends");
        }
        for (int h = 0; h < M; ++h) for (int g = 0; g < M; ++g) for (int f = 0; f < M; ++f)
            if (tgt[f] == src[g] && tgt[g] == src[h] && table[h][table[g][f]] != table[table[h][g]][f]) throw std::invalid_argument("composition is not associative");
        Cat.compose = [table](int g, int f) { return table[g][f]; };
        SimplicialSet S = nerve(Cat, maxdim + 1);
        py::dict res = solve_chain(S.chains(), maxdim, false, "auto", reduce, integral, image, o);
        if (return_sset) res["simplicial_set"] = sset_to_py(S)["simplices"];
        return res;
    }, py::arg("src"), py::arg("tgt"), py::arg("identities"), py::arg("composition"), py::arg("weights") = std::vector<double>(),
       py::arg("maxdim") = 2, py::arg("reduce") = true, py::arg("integral") = false, py::arg("image") = py::none(),
       py::arg("return_sset") = false, COMMON_ARGS,
       "Nerve of a finite category, truncated above maxdim + 1. Morphism m: src[m] -> tgt[m]; identities[o] is the\n"
       "identity of object o; composition[g][f] (or a dict {(g, f): g o f}) is g o f when tgt[f] == src[g].\n"
       "A chain (f_1, ..., f_n) is born at the largest weight of a composite f_j o ... o f_i of consecutive morphisms.");

    m.def("bar_construction", [](std::vector<std::vector<int>> mult, int identity, int maxdim, bool reduce, bool integral,
                                 py::object coeff, bool classify, int threads, std::vector<uint64_t> fp, std::string mode, py::object seed, double target) {
        SolveOpts o = make_opts(coeff, classify, threads, fp, mode, seed, target);
        int n = (int)mult.size();
        for (auto& r : mult) { if ((int)r.size() != n) throw std::invalid_argument("mult must be square"); for (int x : r) if (x < 0 || x >= n) throw std::invalid_argument("mult entry out of range"); }
        for (int a = 0; a < n; ++a) for (int b = 0; b < n; ++b) for (int c = 0; c < n; ++c) if (mult[mult[a][b]][c] != mult[a][mult[b][c]]) throw std::invalid_argument("mult is not associative");
        for (int a = 0; a < n; ++a) if (mult[identity][a] != a || mult[a][identity] != a) throw std::invalid_argument("not an identity");
        SimplicialSet S = bar_construction(mult, identity, maxdim + 1);
        return solve_chain(S.chains(), maxdim, false, "auto", reduce, integral, py::none(), o);
    }, py::arg("mult"), py::arg("identity") = 0, py::arg("maxdim") = 3, py::arg("reduce") = true, py::arg("integral") = true, COMMON_ARGS,
       "B M for a finite monoid M (mult[a][b] = a * b), truncated above maxdim + 1; for a group this is BG, so with\n"
       "integral=True the result carries H_k(G; Z) for k <= maxdim.");

    m.def("laurent", [](py::list cells, int maxdim, bool classify, bool homology) {
        laurent::LaurentComplex X;
        for (auto item : cells) {
            auto t = item.cast<py::tuple>();
            if (t.size() < 2) throw std::invalid_argument("a cell is (dim, birth, [(face, coefficient), ...])");
            std::vector<std::pair<int32_t, laurent::LPoly>> f;
            if (t.size() > 2) for (auto x : t[2].cast<py::list>()) {
                auto fc = x.cast<py::tuple>();
                if (fc.size() != 2) throw std::invalid_argument("a face is (index, coefficient)");
                int32_t y = fc[0].cast<int32_t>();
                if (y < 0 || (size_t)y >= X.size()) throw std::invalid_argument("a face must be an earlier cell");
                f.push_back({y, lpoly_from_py(fc[1])});
            }
            X.add_cell(t[0].cast<int>(), t[1].cast<double>(), f);
        }
        std::string e = X.validate(); if (!e.empty()) throw std::invalid_argument(e);
        int md = maxdim < 0 ? std::max(0, X.max_dim()) : maxdim;
        laurent::LaurentResult L; { py::gil_scoped_release nogil; L = laurent::laurent_persistence(X, md, classify, homology); }
        py::dict res;
        res["coefficients"] = "Z[t, t^-1]"; res["barcode_over"] = "Q(t)";
        res["barcode"] = lbars_to_py(L.bars);
        if (classify) {
            py::list bi;
            for (auto& b : L.bad_irreducibles) {
                py::dict d; d["polynomial"] = laurent::to_string(b.g);
                py::list c; for (auto& x : b.g) c.append(bigint_to_py(x)); d["coefficients"] = c;
                d["irreducible"] = b.irreducible; d["every_irreducible_factor_of"] = !b.irreducible;
                d["barcode"] = lbars_to_py(b.bars); bi.append(d);
            }
            res["bad_irreducibles"] = bi;
            py::dict bp; for (auto& kv : L.bad_primes) bp[py::int_(kv.first)] = lbars_to_py(kv.second);
            res["bad_primes"] = bp;
            py::dict cand; cand["polynomials"] = L.candidate_polys; cand["primes"] = L.candidate_primes; res["candidates"] = cand;
            res["unresolved"] = L.unresolved; res["complete"] = L.complete;
        }
        if (homology) {
            res["homology_over"] = "Q[t, t^-1]";
            res["homology"] = lhom_to_py(L.homology);
            if (classify) {
                py::dict hb; for (auto& kv : L.homology_bad_primes) hb[py::int_(kv.first)] = lhom_to_py(kv.second);
                res["homology_bad_primes"] = hb; res["homology_candidate_primes"] = L.homology_candidate_primes;
            }
        }
        return res;
    }, py::arg("cells"), py::arg("maxdim") = -1, py::arg("classify") = true, py::arg("homology") = true,
       "Filtered based chain complex over Z[t, t^-1]: list of (dim, birth, [(face index, coefficient), ...]), a\n"
       "coefficient being an int, a string such as '1 - t^-1' or a dict {exponent: coefficient}. Returns the barcode\n"
       "over Q(t), the bad irreducibles f (the barcode over Q[t]/(f) differs) and bad primes p (over F_p(t)), the\n"
       "homology of the whole complex over Q[t, t^-1] (free rank and invariant factors) and the primes p where the\n"
       "homology over F_p[t, t^-1] is not the reduction mod p. The Z[t, t^-1]-module itself is not computed.");

    m.def("magnitude", [](std::vector<std::pair<int, int>> edges, int n, bool directed, int lmax, int threads, std::vector<std::pair<int, int>> pairs) {
        for (auto& e : edges) n = std::max(n, std::max(e.first, e.second) + 1);
        MagGraph G(n, edges, directed);
        MagResult R; { py::gil_scoped_release nogil; R = magnitude_homology(G, lmax, threads > 0 ? threads : (int)std::max(1u, std::thread::hardware_concurrency()), pairs); }
        py::dict res, tab;
        for (int l = 0; l <= lmax; ++l) for (int k = 0; k <= l; ++k) {
            const MagEntry& e = R.table[l][k];
            if (!e.generators) continue;
            py::dict d; d["free"] = e.free;
            py::dict t; for (auto& kv : e.torsion) t[py::int_(py::str(kv.first))] = kv.second;
            d["torsion"] = t; d["generators"] = e.generators; d["after_reduction"] = e.residual;
            tab[py::make_tuple(k, l)] = d;
        }
        res["table"] = tab; res["seconds"] = R.seconds; res["pairs"] = R.pairs; res["vertices"] = n; res["directed"] = directed;
        return res;
    }, py::arg("edges"), py::arg("n") = 0, py::arg("directed") = false, py::arg("lmax") = 4, py::arg("threads") = 0,
       py::arg("pairs") = std::vector<std::pair<int, int>>(),
       "Magnitude homology MH_{k,l} over Z of a graph (Hepworth-Willerton), split by (s, t) pairs and lengths:\n"
       "table[(k, l)] = {'free': rank, 'torsion': {order: multiplicity}, 'generators': ..., 'after_reduction': ...}.");

    m.def("simplices", [](std::vector<std::pair<std::vector<int>, double>> simp, int maxdim, py::object coeff, bool classify, int threads,
                          std::vector<uint64_t> fp, std::string mode, py::object seed, double target) {
        SolveOpts o = make_opts(coeff, classify, threads, fp, mode, seed, target);
        SimplicialBuilder B; for (auto& s : simp) B.add(s.first, s.second);
        FilteredComplex C = B.build();
        Prepared P(C); int md_ = maxdim < 0 ? P.D : maxdim; ExplicitBackend be(P, md_, auto_direction(P, md_));
        return solve(be, o);
    }, py::arg("simplices"), py::arg("maxdim") = -1, COMMON_ARGS,
       "Simplicial complex from (vertices, birth) pairs; faces get the minimum birth of their cofaces.");

    m.def("cubical", [](py::array_t<double, py::array::forcecast> top_in, std::vector<bool> periodic, int maxdim, std::string engine,
                        py::object coeff, bool classify, int threads, std::vector<uint64_t> fp, std::string mode, py::object seed, double target) {
        SolveOpts o = make_opts(coeff, classify, threads, fp, mode, seed, target);
        CubicalGrid g;
        const py::ssize_t nd = top_in.ndim();
        if (periodic.empty()) periodic.assign(nd, false);
        if ((py::ssize_t)periodic.size() != nd) throw std::invalid_argument("periodic must have one flag per axis");
        // The engines read axis 0 fastest (Fortran order). A C-ordered array is read as is, with its axes
        // reversed: that is the transposed grid, an isomorphic cubical complex, so the barcode is the same
        // (only ties among equal values are broken in another order, which no positive-length bar sees and
        // zero-length bars are counted per value). This avoids a transposing copy.
        bool fortran = top_in.flags() & py::array::f_style;
        py::array_t<double> top = fortran ? py::array_t<double>(top_in)
                                          : py::array_t<double>(py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(top_in));
        for (py::ssize_t i = 0; i < nd; ++i) {
            py::ssize_t a = fortran ? i : nd - 1 - i;
            g.ext.push_back(top.shape(a)); g.periodic.push_back(periodic[a]);
        }
        bool planar = engine == "planar" || ((engine == "auto" || engine == "implicit") && CubicalPlanar::supports(g));
        if (planar && CubicalPlanar::supports(g)) g.view = top.data();   // no copy: the planar engine reads the array in place
        else g.top.assign(top.data(), top.data() + top.size());
        if (engine == "explicit")                             // NaN = absent cube (the implicit engines test isfinite)
            for (auto& x : g.top) if (std::isnan(x)) x = std::numeric_limits<double>::infinity();
        if (engine == "explicit") {
            FilteredComplex C = cubical_from_grid(g);
            Prepared P(C); int md_ = maxdim < 0 ? P.D : maxdim; ExplicitBackend be(P, md_, auto_direction(P, md_));
            return solve(be, o);
        }
        if (engine != "auto" && engine != "implicit" && engine != "planar" && engine != "morse" && engine != "classic")
            throw std::invalid_argument("engine must be 'auto', 'planar', 'morse', 'classic' or 'explicit'");
        CubicalBackend be(g, maxdim, engine);
        return solve(be, o);
    }, py::arg("top"), py::arg("periodic") = std::vector<bool>(), py::arg("maxdim") = -1, py::arg("engine") = "auto", COMMON_ARGS,
       "Cubical complex from an n-d array of top-cell values (NaN or inf = absent); optional periodic axes.");
}
