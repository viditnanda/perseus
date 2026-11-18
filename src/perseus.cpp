// Perseus 6 -- persistent homology over Q with a certificate for every prime, or over any F_p.
#include "perseus/io.hpp"
#include "perseus/explicit_engine.hpp"
#include "perseus/rips_explicit.hpp"
#include "perseus/rips_engine.hpp"
#include "perseus/certify.hpp"
#include "perseus/collapse.hpp"
#include "perseus/rips_sparse.hpp"
#include "perseus/cubical_engine.hpp"
#include <cstdio>
#include <cstring>
#include <string>
#include <thread>
#include <chrono>

using namespace perseus;

static const char* USAGE =
"Perseus 6: persistent homology over Q with every bad prime certified, or over any F_p.\n"
"\n"
"usage: perseus <format> <input> [<output prefix>] [options]\n"
"\n"
"formats:  points | distance | lower-distance | bin32 | bin64    Vietoris-Rips of a point cloud / distances\n"
"          flag                                                  clique filtration of a weighted graph\n"
"          cubical | npy | cubes                                 cubical complexes (dense text, numpy, sparse)\n"
"          simplices | cells                                     simplicial / general cell complexes\n"
"          (old Perseus 4/5 files: convert them with tools/convert_legacy.py)\n"
"\n"
"options:\n"
"  --coeff Q|p      coefficients: Q (default; certified bad primes) or a prime p\n"
"  --mode m         Q certificate: auto (default: exact int64 arithmetic, a proof; if it overflows, as fast)\n"
"                   | fast (random 62-bit primes; error probability bounded and reported)\n"
"                   | proved (exact rational arithmetic, bigint where int64 overflows)\n"
"  --proved         same as --mode proved\n"
"  --target t       fast mode: use enough primes for an error bound below 2^-t (default 64)\n"
"  --seed hex       fast mode: reproduce a run from its reported seed (up to 64 hex digits).\n"
"                   A run with a supplied seed is not certified: the error bound is a probability\n"
"                   over a seed drawn after the input is fixed (default: fresh 256-bit seed)\n"
"  --fp p,q,...     also report the barcodes over these prime fields\n"
"  --maxdim d       highest homological degree (Rips/flag default 1, otherwise the top dimension)\n"
"  --threshold r    Rips/flag: ignore edges longer than r (default: enclosing radius)\n"
"  --periodic a,b   cubical/npy: periodic axes (0-based)\n"
"  --no-classify    Q mode: list candidate primes but do not rerun to classify them\n"
"  --zero-length    keep zero-length intervals (explicit complexes)\n"
"  --threads n      threads for classifying candidate primes (default: all cores)\n"
"  --collapse c     Rips/flag: auto | on | full | off   iterated edge collapse (Boissonnat-Pritam,\n"
"                   Glisse-Pritam); auto = on for maxdim >= 2 when it removes enough edges\n"
"  --engine e       Rips/flag: auto | implicit | sparse (explicit flag complex from the radius graph;\n"
"                   auto picks sparse when the flag complex is small next to n^2: n >= 2048,\n"
"                   average degree <= 64, at most n^2/50 simplices)\n"
"                   cubical/npy: auto | planar (2D, no periodic axes: union-finds and duality)\n"
"                   | morse (discrete gradient; dimensions 2-4, no periodic axes) | classic\n"
"                   | explicit (build the complex); auto = planar in 2D, else morse where it applies\n"
"  --direction d    explicit complexes: auto | cohomology | homology | forward\n"
"  --print fmt      ripser | perseus | json | none   (stdout; default ripser unless a prefix is given)\n"
"  --stats          timing and engine statistics on stderr\n"
"\n"
"outputs with a prefix: prefix_k.txt (birth death, -1 = infinite), prefix_betti.txt, and in Q mode\n"
"prefix_cert.json; with --fp, prefix_F<p>_k.txt.\n";

struct Opts {
    std::string format, input, prefix, coeff, engine = "auto", print, direction = "auto", collapse = "auto";
    std::vector<uint64_t> fp;
    std::vector<int> periodic;
    int maxdim = -1, threads = 0;
    double threshold = -1, target = 64;
    std::string mode = "auto"; std::string seed;
    bool classify = true, zero = false, stats = false;
};

static std::string json_bars(const Barcode& bc) {
    std::string s = "[";
    for (size_t k = 0; k < bc.size(); ++k) {
        s += (k ? "," : ""); s += "[";
        for (size_t i = 0; i < bc[k].size(); ++i) {
            s += (i ? "," : ""); s += "[" + fmt_num(bc[k][i].birth) + "," + (std::isinf(bc[k][i].death) ? std::string("null") : fmt_num(bc[k][i].death)) + "]";
        }
        s += "]";
    }
    return s + "]";
}
static std::string guarantee_text(const Certificate& c) {
    char b[200];
    if (!c.barcode_ok) return "NOT certified (" + c.note + ")";
    if (c.mode == "proved") return "proved (exact arithmetic: " + c.arithmetic + ")";
    if (c.mode == "unimodular") return "proved (total unimodularity: the barcode is the same over Q and every F_p)";
    if (std::isinf(c.log2_error)) snprintf(b, sizeof b, "probabilistic, error probability 0 (no integer can have a 62-bit prime factor)");
    else snprintf(b, sizeof b, "probabilistic, error probability < 2^%.1f", std::ceil(c.log2_error * 10) / 10);
    if (c.seed_user) return "NOT certified (reproduced from a supplied seed; a fresh seed would give " + std::string(b).substr(15) + ")";
    return std::string(b) + (c.note.empty() ? "" : " (" + c.note + ")");
}
static std::string json_cert(const Barcode& q, const Certificate& c, int coeff_q) {
    std::string s = "{\n  \"coefficients\": \"" + std::string(coeff_q ? "Q" : "F_p") + "\",\n";
    s += "  \"barcode\": " + json_bars(q) + ",\n";
    if (coeff_q) {
        s += "  \"mode\": \"" + c.mode + "\",\n";
        s += "  \"guarantee\": {\"barcode_certified\": " + std::string(c.certified() ? "true" : "false");
        if (c.mode == "fast") {
            char b[64]; snprintf(b, sizeof b, "%.2f", c.log2_error);
            s += ", \"log2_error_bound\": " + (std::isinf(c.log2_error) ? std::string("null") : std::string(b)) + ", \"seed\": \"" + c.seed + "\", \"seed_supplied\": " + (c.seed_user ? "true" : "false") +
                 ", \"lanes\": " + std::to_string(c.lanes) + ", \"retries\": " + std::to_string(c.retries) + ", \"moduli\": [";
            for (size_t i = 0; i < c.moduli.size(); ++i) s += (i ? "," : "") + std::to_string(c.moduli[i]);
            s += "]";
        }
        if (c.mode == "proved" || c.mode == "unimodular") s += ", \"arithmetic\": \"" + c.arithmetic + "\"";
        if (!c.note.empty()) s += ", \"note\": \"" + c.note + "\"";
        s += "},\n";
        s += "  \"bad_primes_status\": \"" + std::string(c.seed_user && c.barcode_ok ? "conditional on the supplied seed" : !c.candidates_complete() ? "incomplete" : (c.classified ? "complete" : "candidates only")) + "\",\n";
    }
    s += "  \"complete\": " + std::string(c.complete() || !coeff_q ? "true" : "false") + ",\n";
    s += "  \"nonunit_pivots\": [";
    for (size_t i = 0; i < c.values.size(); ++i) {
        auto& v = c.values[i];
        s += (i ? "," : ""); s += "\n    {\"degree\": " + std::to_string(v.k) + ", \"birth\": " + fmt_num(v.birth) + ", \"death\": " + fmt_num(v.death) +
             ", \"value\": " + (v.ok ? "\"" + v.num.str() + (v.den == BigInt(1) ? "" : "/" + v.den.str()) + "\"" : std::string("null")) + "}";
    }
    s += "],\n  \"candidates\": [";
    { bool f = true; for (auto& kv : c.candidates) { s += (f ? "" : ","); s += std::to_string(kv.first); f = false; }
      for (auto& kv : c.big_candidates) { s += (f ? "" : ","); s += kv.first.str(); f = false; } }
    s += "],\n  \"bad_primes\": [";
    for (size_t i = 0; i < c.bad.size(); ++i) {
        auto& b = c.bad[i];
        s += (i ? "," : ""); s += "\n    {\"p\": " + b.str() + (b.group ? ", \"every_prime_factor_of\": true" : (b.is_big && b.big.bits() > 64 ? ", \"probable_prime\": true" : "")) + ", \"first_degree\": " + std::to_string(b.k0) +
             ", \"only_over_Fp\": " + json_bars(b.only_p) + ", \"only_over_Q\": " + json_bars(b.only_q) + "}";
    }
    s += "],\n  \"unresolved_factors\": [";
    for (size_t i = 0; i < c.unresolved.size(); ++i) { s += (i ? "," : ""); s += "\"" + c.unresolved[i].str() + "\""; }
    s += "],\n  \"unreconstructed_pivots\": " + std::to_string(c.unreconstructed) + ",\n";
    char tb[200]; snprintf(tb, sizeof tb, "  \"seconds\": {\"Q\": %.4f, \"reconstruction\": %.4f, \"classification\": %.4f}\n}\n", c.t_q, c.t_recon, c.t_classify);
    return s + tb;
}

static int choose_direction(const Prepared& P, const std::string& d0, int maxdim) {
    std::string d = d0;
    if (d == "auto" && getenv("PERSEUS_DIRECTION")) d = getenv("PERSEUS_DIRECTION");   // testing hook
    if (d == "homology") return ExplicitBackend::HOMOLOGY;
    if (d == "cohomology") return ExplicitBackend::COHOMOLOGY;
    if (d == "forward") return ExplicitBackend::FORWARD;
    return auto_direction(P, maxdim);
}

template <class B>
static int solve(const B& be, const Opts& o) {
    auto T0 = std::chrono::steady_clock::now();
    std::string coeff = o.coeff.empty() ? "Q" : o.coeff;
    Barcode main_bars; Certificate cert; bool qmode = (coeff == "Q" || coeff == "q" || coeff == "0");
    if (qmode) {
        CertifyOptions co; co.threads = o.threads; co.classify = o.classify;
        co.mode = o.mode; co.seed = o.seed; co.target_log2 = -o.target;
        auto res = certify(be, co);
        main_bars = res.first.bars; cert = res.second;
    } else {
        uint64_t p = std::stoull(coeff);
        if (!is_prime64(p)) { fprintf(stderr, "--coeff: %llu is not a prime\n", (unsigned long long)p); return 2; }
        main_bars = run_field(be, p, 0, be.maxdim(), nullptr).bars;
    }
    auto T1 = std::chrono::steady_clock::now();
    std::string pr = o.print.empty() ? (o.prefix.empty() ? "ripser" : "none") : o.print;
    if (pr == "ripser") print_ripser(stdout, main_bars);
    else if (pr == "perseus") { for (size_t k = 0; k < main_bars.size(); ++k) for (auto& b : main_bars[k]) printf("%zu %s %s\n", k, fmt_num(b.birth).c_str(), std::isinf(b.death) ? "-1" : fmt_num(b.death).c_str()); }
    else if (pr == "json") fputs(json_cert(main_bars, cert, qmode).c_str(), stdout);
    if (!o.prefix.empty()) {
        write_files(o.prefix, main_bars);
        if (qmode) { FILE* f = fopen((o.prefix + "_cert.json").c_str(), "w"); fputs(json_cert(main_bars, cert, 1).c_str(), f); fclose(f); }
    }
    if (qmode) {
        fprintf(stderr, "Q barcode: %s\n", guarantee_text(cert).c_str());
        fprintf(stderr, "bad primes: ");
        if (!cert.barcode_ok) fprintf(stderr, "not determined");
        else if (cert.seed_user) fprintf(stderr, "(conditional on the supplied seed) ");
        if (!cert.barcode_ok) {}
        else if (!cert.classified) { fprintf(stderr, "not classified; candidates:"); for (auto& kv : cert.candidates) fprintf(stderr, " %llu", (unsigned long long)kv.first); for (auto& kv : cert.big_candidates) fprintf(stderr, " %s", kv.first.str().c_str()); if (cert.candidates.empty() && cert.big_candidates.empty()) fprintf(stderr, " none"); }
        else { if (cert.bad.empty()) fprintf(stderr, "none"); for (auto& b : cert.bad) fprintf(stderr, b.group ? "[every prime factor of %s] " : "%s ", b.str().c_str()); }
        fprintf(stderr, " (%zu non-unit pivots", cert.values.size());
        if (cert.unreconstructed) fprintf(stderr, "; %d values not reconstructed, list INCOMPLETE", cert.unreconstructed);
        if (!cert.unresolved.empty()) fprintf(stderr, "; %zu unresolved factors, list INCOMPLETE", cert.unresolved.size());
        fprintf(stderr, ")\n");
    }
    for (uint64_t p : o.fp) {
        if (!is_prime64(p)) { fprintf(stderr, "--fp: %llu is not a prime\n", (unsigned long long)p); continue; }
        Barcode b = run_field(be, p, 0, be.maxdim(), nullptr).bars;
        if (!o.prefix.empty()) write_files(o.prefix + "_F" + std::to_string(p), b);
        else { printf("# coefficients F_%llu\n", (unsigned long long)p); print_ripser(stdout, b); }
    }
    if (o.stats) {
        double tot = std::chrono::duration<double>(std::chrono::steady_clock::now() - T0).count();
        if (qmode) fprintf(stderr, "  [solve %.3f s: Q run %.3f s, reconstruction %.3f s (%d lanes), classification %.3f s (%zu candidates)]\n",
                           tot, cert.t_q, cert.t_recon, cert.lanes, cert.t_classify, cert.candidates.size());
        else fprintf(stderr, "  [solve %.3f s: field run %.3f s]\n", tot, std::chrono::duration<double>(T1 - T0).count());
    }
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 3) { fputs(USAGE, stderr); return 1; }
    Opts o; std::vector<std::string> pos;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto val = [&]() -> std::string { if (i + 1 >= argc) { fprintf(stderr, "%s needs a value\n", a.c_str()); exit(1); } return argv[++i]; };
        if (a == "--coeff") o.coeff = val();
        else if (a == "--fp") { std::string v = val(); for (char& c : v) if (c == ',') c = ' '; for (auto& t : split_ws(v)) o.fp.push_back(std::stoull(t)); }
        else if (a == "--maxdim") o.maxdim = std::stoi(val());
        else if (a == "--periodic") { std::string v = val(); for (char& ch : v) if (ch == ',') ch = ' '; for (auto& t : split_ws(v)) o.periodic.push_back(std::stoi(t)); }
        else if (a == "--threshold") o.threshold = std::stod(val());
        else if (a == "--no-classify") o.classify = false;
        else if (a == "--mode") { o.mode = val(); if (o.mode != "auto" && o.mode != "fast" && o.mode != "proved") { fprintf(stderr, "--mode: auto, fast or proved\n"); return 1; } }
        else if (a == "--proved") o.mode = "proved";
        else if (a == "--target") o.target = std::stod(val());
        else if (a == "--seed") { o.seed = val(); parse_seed(o.seed); }
        else if (a == "--zero-length") o.zero = true;
        else if (a == "--threads") o.threads = std::stoi(val());
        else if (a == "--engine") o.engine = val();
        else if (a == "--direction") o.direction = val();
        else if (a == "--collapse") o.collapse = val();
        else if (a == "--print") o.print = val();
        else if (a == "--stats") o.stats = true;
        else if (a == "-h" || a == "--help") { fputs(USAGE, stdout); return 0; }
        else pos.push_back(a);
    }
    if (pos.size() < 2) { fputs(USAGE, stderr); return 1; }
    o.format = pos[0]; o.input = pos[1]; if (pos.size() > 2) o.prefix = pos[2];
    for (char& c : o.format) c = (char)tolower(c);
    if (o.threads <= 0) o.threads = std::max(1u, std::thread::hardware_concurrency());
    auto TR0 = std::chrono::steady_clock::now();
    auto lap = [&](const char* what) {
        if (!o.stats) return;
        auto t = std::chrono::steady_clock::now();
        fprintf(stderr, "  [%s %.3f s]\n", what, std::chrono::duration<double>(t - TR0).count()); TR0 = t;
    };
    try {
        Input in; std::string f = o.format;
        bool binary_file = (f == "bin32" || f == "bin64" || f == "npy");
        std::ifstream is; if (!binary_file) { is.open(o.input); if (!is) { fprintf(stderr, "cannot open %s\n", o.input.c_str()); return 2; } }
        auto periodic_flags = [&](int d) {
            std::vector<bool> p(d, false);
            for (int a : o.periodic) { if (a < 0 || a >= d) throw std::runtime_error("--periodic: axis out of range"); p[a] = true; }
            return o.periodic.empty() ? std::vector<bool>() : p;
        };
        if (f == "points") { in.is_rips = true; in.rips = read_points(is); }
        else if (f == "lower-distance" || f == "lower") { in.is_rips = true; in.rips = read_lower(is); }
        else if (f == "distance" || f == "matrix") { in.is_rips = true; in.rips = read_matrix(is); }
        else if (f == "bin32" || f == "bin64") { in.is_rips = true; in.rips = read_binary(o.input, f == "bin64"); }
        else if (f == "flag") { in.is_rips = true; in.rips = read_flag(is); in.rips.threshold = std::numeric_limits<double>::infinity(); }
        else if (f == "cubical") {
            // sizes line decides the dimension; periodic flags are validated against it
            std::streampos pos = is.tellg();
            std::string line; int d = 0;
            while (std::getline(is, line)) { if (line.empty() || line[0] == '#') continue; for (char& ch : line) if (ch == ',') ch = ' '; d = (int)split_ws(line).size(); break; }
            is.clear(); is.seekg(pos);
            CubicalGrid g = read_cubical_grid(is, periodic_flags(d));
            if (o.engine == "explicit") in.complex = cubical_from_grid(g); else { in.is_grid = true; in.grid = std::move(g); }
        }
        else if (f == "npy") {
            NpyArray a = read_npy_array(o.input); int d = (int)a.shape.size();
            CubicalGrid g = grid_from_npy(std::move(a), periodic_flags(d));
            if (o.engine == "explicit") in.complex = cubical_from_grid(g); else { in.is_grid = true; in.grid = std::move(g); }
        }
        else if (f == "cubes") in.complex = read_cubes(is);
        else if (f == "simplices") in.complex = read_simplices(is);
        else if (f == "cells") in.complex = read_cells(is);
        else {
            fprintf(stderr, "unknown format '%s'%s\n\n%s", f.c_str(),
                    (f == "simtop" || f == "nmfsimtop" || f == "cubtop" || f == "scubtop" || f == "cellcomp" || f == "distmat" || f == "corrmat" || f == "rips" || f == "brips")
                        ? " (an old Perseus format: convert the file with tools/convert_legacy.py)" : "", USAGE);
            return 1;
        }
        lap("read and build input");

        if (in.is_rips) {
            RipsInput& R = in.rips;
            if (o.threshold >= 0) R.threshold = o.threshold;
            int md = o.maxdim >= 0 ? o.maxdim : 1;
            if (!std::isfinite(R.threshold)) R.threshold = enclosing_radius_any(R);
            // sparse route: build the radius graph and the flag complex explicitly. Taken for vertex births,
            // for --engine explicit|sparse, and automatically when the flag complex is small next to the
            // implicit engine's n x n tables (see prefer_sparse).
            bool sparse_src = !R.dense();
            std::vector<WEdge> E;
            bool want_sparse = o.engine == "explicit" || o.engine == "sparse" || !R.vertex_birth.empty();
            if (!want_sparse && o.engine == "auto") want_sparse = prefer_sparse(R, E, md + 1);
            if (want_sparse) {
                if (E.empty()) E = radius_edges(R, R.threshold);
                lap("radius graph");
                FilteredComplex C = flag_complex(R.n, E, R.vertex_birth, md + 1);
                std::vector<WEdge>().swap(E); RipsInput().pts.swap(R.pts);
                if (o.stats) fprintf(stderr, "  [sparse flag complex: %zu simplices up to dimension %d]\n", C.size(), md + 1);
                lap("flag complex");
                Prepared P(C); C = FilteredComplex();
                ExplicitBackend be(P, md, choose_direction(P, o.direction, md)); be.keep_zero_length = o.zero;
                lap("prepare");
                return solve(be, o);
            }
            std::vector<WEdge>().swap(E);
            if (sparse_src) R.densify();
            // edge collapse pays for itself in degrees >= 2; for H0/H1 the active H1 engine is cheaper
            if (o.collapse == "on" || o.collapse == "full" || (o.collapse == "auto" && md >= 2 && R.n >= 32)) {
                CollapseStats cs;
                bool aut = (o.collapse == "auto");
                RipsInput C = edge_collapse(R, o.collapse == "full" ? 1000 : 32, cs, o.collapse == "full" ? 0.0 : 0.02, aut ? 0.10 : 0.0);
                if (o.stats) {
                    if (cs.used) fprintf(stderr, "  [edge collapse: %zu -> %zu edges (%zu shifted), %d passes, %.3f s]\n",
                                         cs.edges_before, cs.edges_after, cs.shifted, cs.passes, cs.seconds);
                    else fprintf(stderr, "  [edge collapse not worth it here, abandoned after %.3f s]\n", cs.seconds);
                }
                if (cs.used) R.lower.swap(C.lower);
            }
            RipsBackend be(R, md); be.keep_zero_length = o.zero;
            std::vector<double>().swap(R.lower);          // the backend keeps its own edge data
            lap("rips setup");
            return solve(be, o);
        }
        if (in.is_grid) {
            CubicalBackend be(in.grid, o.maxdim, o.engine, o.zero); be.keep_zero_length = o.zero;
            CubicalGrid().top.swap(in.grid.top);
            lap("cubical setup");
            if (o.stats) {
                fprintf(stderr, "  [cubical engine: %s", be.name());
                if (const CubicalMorse* m = be.morse()) {
                    fprintf(stderr, "; critical cells");
                    for (int k = 0; k <= (int)in.grid.ext.size(); ++k) fprintf(stderr, " %zu/%zu", m->critical(k), m->count(k));
                    fprintf(stderr, "; %zu block patterns, %lld blocks above the minimum", m->patterns(), m->imperfect_blocks());
                }
                fprintf(stderr, "]\n");
            }
            return solve(be, o);
        }
        std::string err = in.complex.validate(false);
        if (!err.empty()) { fprintf(stderr, "invalid complex: %s\n", err.c_str()); return 3; }
        lap("validate");
        Prepared P(in.complex);
        lap("prepare");
        ExplicitBackend be(P, o.maxdim >= 0 ? o.maxdim : P.D, choose_direction(P, o.direction, o.maxdim >= 0 ? o.maxdim : P.D)); be.keep_zero_length = o.zero;
        return solve(be, o);
    } catch (std::exception& e) {
        fprintf(stderr, "error: %s\n", e.what()); return 2;
    }
}
