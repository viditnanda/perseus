// Perseus 6 -- persistent homology over Q with a certificate for every prime, or over any F_p.
#include "perseus/io.hpp"
#include "perseus/explicit_engine.hpp"
#include "perseus/rips_explicit.hpp"
#include "perseus/rips_engine.hpp"
#include "perseus/certify.hpp"
#include "perseus/collapse.hpp"
#include "perseus/rips_components.hpp"
#include "perseus/cubical_engine.hpp"
#include "perseus/chain.hpp"
#include "perseus/simplicial_set.hpp"
#include "perseus/torsion.hpp"
#include "perseus/magnitude.hpp"
#include "perseus/pc_group.hpp"
#include "perseus/sylow.hpp"
#include "perseus/pc_maps.hpp"
#include "perseus/laurent.hpp"
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
"          sset                                                  finite simplicial set (normalised chains)\n"
"          magnitude                                             magnitude homology of a graph (edge list)\n"
"          laurent                                               chain complex over Z[t, t^-1] (Novikov / twisted)\n"
"          pcgroup                                               finite solvable group: H_*(G; Z) and the mod-p\n"
"                                                                cohomology ring (tools/export_pcgroup.g writes these)\n"
"          pcseries                                              a chain of homomorphisms of finite solvable groups: the\n"
"                                                                maps on H_*(-; F_p) and their persistence (Ellis-King;\n"
"                                                                tools/export_pcseries.g writes normal and subgroup series)\n"
"          sylow                                                 any finite group, through its Sylow subgroups: H_*(G; Z)\n"
"                                                                and mod-p cohomology rings (tools/export_sylow.g)\n"
"          (old Perseus 4/5 files: convert them with tools/convert_legacy.py)\n"
"\n"
"options:\n"
"  --coeff Q|p      coefficients: Q (default; certified bad primes) or a prime p\n"
"  --mode m         Q certificate: proved (default: exact rationals, int64 where they fit, big ones\n"
"                   where they do not; a proof) | auto (the same as proved)\n"
"                   | fast (random 62-bit primes: NOT a proof; the error probability is bounded and\n"
"                   reported, and every output says so)\n"
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
"  --threads n      threads for prime classification and independent Rips components (default: all cores)\n"
"  --collapse c     Rips/flag: auto | on | full | off   iterated edge collapse (Boissonnat-Pritam,\n"
"                   Glisse-Pritam); auto = on for maxdim >= 2 when it removes enough edges\n"
"  --engine e       Rips/flag: auto | implicit | sparse (explicit flag complex from the radius graph;\n"
"                   auto picks sparse when the flag complex is small next to n^2: n >= 2048,\n"
"                   average degree <= 64, at most n^2/50 simplices)\n"
"                   cubical/npy: auto | planar (2D, no periodic axes: union-finds and duality)\n"
"                   | morse (discrete gradient; dimensions 2-4, no periodic axes) | classic\n"
"                   | explicit (build the complex); auto = planar in 2D, else morse where it applies\n"
"  --direction d    explicit complexes: auto | cohomology | homology | forward\n"
"  --integral       explicit complexes: also the integral homology of the whole complex (torsion\n"
"                   by Smith forms) and, in Q mode, the F_p barcode of every bad prime\n"
"  --image s,t      explicit complexes: the image of H_*(K_s; Z) -> H_*(K_t; Z)\n"
"  --torsion        explicit complexes, Q mode: the torsion of H_q(K_t; Z) at every t where it changes,\n"
"                   with its invariant factors and p-primary exponents (Howell sweeps over Z/p^k)\n"
"  --classify-by m  explicit complexes: auto (sweeps when a sample of reruns predicts they are cheaper)\n                   | sweep (always batched Howell sweeps over Z/m) | rerun (one\n"
"                   rerun per candidate prime); big candidates (>= 2^63) are always rerun\n"
"  --no-reduce      skip the algebraic Morse reduction (unit pairs with equal births); it runs for\n"
"                   cells/sset, and for other explicit complexes with --integral, --image or --torsion\n"
"  --lmax l         magnitude: lengths 0..l (default 4)\n"
"  --directed       magnitude: the edges are directed\n"
"  --pairs s,t,...  magnitude: only these (start, end) pairs\n"
"  --ring p         pcgroup, sylow: also the mod-p cohomology ring through degree maxdim (default maxdim 4)\n"
"  --no-homology    pcgroup, sylow: skip the integral homology\n"
"  --all-products   pcgroup, sylow --ring: every product of basis classes (default: products of the generators with every class)\n""  --keep-pcgs      pcgroup, sylow, pcseries: run on the given pcgs (default: the cheapest of it, its prime refinement and\n"
"                   its coarsening, by the work of a short trial resolution)\n"
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
    std::string mode = "proved"; std::string seed;
    bool classify = true, zero = false, stats = false;
    bool integral = false, image = false, reduce = true, directed = false, torsion = false;
    std::string classify_by = "auto";
    double img_s = 0, img_t = 0; int lmax = 4; int ring = 0; bool pc_homology = true; bool all_products = false; bool keep_pcgs = false;
    std::vector<std::pair<int, int>> pairs;
};

// integral outputs for an explicit complex (computed on the given, possibly reduced, complex)
struct IntegralSpec { const FilteredComplex* C = nullptr; };
static std::string json_integral(const IntegralHomology& H) {
    std::string s = "[";
    for (size_t k = 0; k < H.free.size(); ++k) {
        s += (k ? "," : ""); s += "{\"degree\": " + std::to_string(k) + ", \"free\": " + std::to_string(H.free[k]) + ", \"torsion\": [";
        for (size_t i = 0; i < H.torsion[k].size(); ++i) s += (i ? "," : "") + H.torsion[k][i].str();
        s += "]}";
    }
    return s + "]";
}
static std::string text_integral(const IntegralHomology& H) {
    std::string s;
    for (size_t k = 0; k < H.free.size(); ++k) {
        s += "  H_" + std::to_string(k) + " = ";
        bool any = false;
        if (H.free[k]) { s += "Z" + (H.free[k] > 1 ? "^" + std::to_string(H.free[k]) : std::string()); any = true; }
        for (auto& t : H.torsion[k]) { s += std::string(any ? " + " : "") + "Z/" + t.str(); any = true; }
        if (!any) s += "0";
        s += "\n";
    }
    return s;
}

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
static bool qmode_of(const Opts& o) { std::string c = o.coeff.empty() ? "Q" : o.coeff; return c == "Q" || c == "q" || c == "0"; }
static std::string json_cert(const Barcode& q, const Certificate& c, int coeff_q, const std::string& extra = "") {
    std::string s = "{\n  \"coefficients\": \"" + std::string(coeff_q ? "Q" : "F_p") + "\",\n";
    s += "  \"barcode\": " + json_bars(q) + ",\n";
    if (coeff_q) {
        s += "  \"mode\": \"" + c.mode + "\",\n";
        s += "  \"proof\": " + std::string(c.mode != "fast" && c.certified() ? "true" : "false") + ",\n";
        if (c.mode == "fast") s += "  \"warning\": \"--mode fast: probabilistic, not a proof\",\n";
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
    char tb[200]; snprintf(tb, sizeof tb, "  \"seconds\": {\"Q\": %.4f, \"reconstruction\": %.4f, \"classification\": %.4f}", c.t_q, c.t_recon, c.t_classify);
    return s + tb + extra + "\n}\n";
}

static int choose_direction(const Prepared& P, const std::string& d0, int maxdim) {
    std::string d = d0;
    if (d == "auto" && getenv("PERSEUS_DIRECTION")) d = getenv("PERSEUS_DIRECTION");   // testing hook
    if (d == "homology") return ExplicitBackend::HOMOLOGY;
    if (d == "cohomology") return ExplicitBackend::COHOMOLOGY;
    if (d == "forward") return ExplicitBackend::FORWARD;
    return auto_direction(P, maxdim);
}


static std::string json_torsion(const torsion::TorsionResult& T) {
    std::string s = "{\"degrees\": [";
    for (size_t i = 0; i < T.degrees.size(); ++i) {
        auto& tr = T.degrees[i];
        s += std::string(i ? ", " : "") + "{\"degree\": " + std::to_string(tr.degree) + ", \"steps\": [";
        for (size_t j = 0; j < tr.steps.size(); ++j) {
            auto& st = tr.steps[j];
            s += std::string(j ? ", " : "") + "\n      {\"value\": " + fmt_num(st.value) + ", \"invariant_factors\": [";
            auto d = torsion::invariant_factors(T, st);
            for (size_t t = 0; t < d.size(); ++t) s += std::string(t ? ", " : "") + "\"" + d[t].str() + "\"";
            s += "], \"primary\": {";
            for (size_t t = 0; t < st.primary.size(); ++t) {
                s += std::string(t ? ", " : "") + "\"" + T.factors[st.primary[t].first].p.str() + "\": [";
                for (size_t e = 0; e < st.primary[t].second.size(); ++e) s += std::string(e ? ", " : "") + std::to_string(st.primary[t].second[e]);
                s += "]";
            }
            s += "}}";
        }
        s += "]}";
    }
    s += "], \"groups\": [";
    bool f = true; for (auto& F : T.factors) if (F.group) { s += std::string(f ? "" : ", ") + "\"" + F.p.str() + "\""; f = false; }
    char b[200]; snprintf(b, sizeof b, "], \"complete\": %s, \"sweeps\": %d, \"seconds\": %.4f", T.complete ? "true" : "false", T.sweeps, T.seconds);
    s += b;
    if (!T.note.empty()) s += ", \"note\": \"" + T.note + "\"";
    return s + "}";
}
static std::string text_torsion(const torsion::TorsionResult& T) {
    std::string s;
    for (auto& tr : T.degrees) {
        s += "torsion of H_" + std::to_string(tr.degree) + "(K_t; Z), where it changes:" + (tr.steps.empty() ? " always 0\n" : "\n");
        for (auto& st : tr.steps) {
            auto d = torsion::invariant_factors(T, st);
            s += "  t = " + fmt_num(st.value) + ": ";
            if (d.empty()) s += "0";
            for (size_t t = 0; t < d.size(); ++t) s += std::string(t ? " + " : "") + "Z/" + d[t].str();
            s += "\n";
        }
    }
    for (auto& F : T.factors) if (F.group) s += "  (" + F.p.str() + " is a composite whose prime factors all behave alike)\n";
    if (!T.complete) s += "  INCOMPLETE: " + T.note + "\n";
    return s;
}

static void output_barcode(const Opts& o, const Barcode& main_bars, const Certificate& cert, bool qmode, const std::string& ijson) {
    std::string pr = o.print.empty() ? (o.prefix.empty() ? "ripser" : "none") : o.print;
    if (pr == "ripser") print_ripser(stdout, main_bars);
    else if (pr == "perseus") { for (size_t k = 0; k < main_bars.size(); ++k) for (auto& b : main_bars[k]) printf("%zu %s %s\n", k, fmt_num(b.birth).c_str(), std::isinf(b.death) ? "-1" : fmt_num(b.death).c_str()); }
    else if (pr == "json") fputs(json_cert(main_bars, cert, qmode, ijson).c_str(), stdout);
    if (!o.prefix.empty()) {
        write_files(o.prefix, main_bars);
        if (qmode) { FILE* f = fopen((o.prefix + "_cert.json").c_str(), "w"); fputs(json_cert(main_bars, cert, 1, ijson).c_str(), f); fclose(f); }
    }
}

template <class B>
static int solve(const B& be, const Opts& o, const IntegralSpec& ispec = IntegralSpec()) {
    auto T0 = std::chrono::steady_clock::now();
    std::string coeff = o.coeff.empty() ? "Q" : o.coeff;
    Barcode main_bars; Certificate cert; bool qmode = (coeff == "Q" || coeff == "q" || coeff == "0");
    const Prepared* prep = prepared_of(be);
    torsion::Cache tcache;
    if (qmode) {
        CertifyOptions co; co.threads = o.threads; co.classify = o.classify;
        co.mode = o.mode; co.seed = o.seed; co.target_log2 = -o.target;
        if (prep && o.classify_by != "rerun") {
            int nt = o.threads > 0 ? o.threads : (int)std::max(1u, std::thread::hardware_concurrency());
            co.sweep_policy = o.classify_by == "sweep" ? "sweep" : "auto";
            co.sweep = [&, nt](const Certificate& c, const Barcode& qb, int kmax, double budget, const std::set<uint64_t>& skip, std::vector<BadPrime>& bad, std::string& diag) {
                return torsion::classify(*prep, c, qb, kmax, nt, be.keep_zero_length, bad, tcache, &diag, budget, skip); };
        }
        auto res = certify(be, co);
        main_bars = res.first.bars; cert = res.second;
    } else {
        uint64_t p = std::stoull(coeff);
        if (!is_prime64(p)) { fprintf(stderr, "--coeff: %llu is not a prime\n", (unsigned long long)p); return 2; }
        main_bars = run_field(be, p, 0, be.maxdim(), nullptr).bars;
    }
    auto T1 = std::chrono::steady_clock::now();
    // integral outputs: homology of the whole complex, the image of H(K_s) -> H(K_t), F_p barcodes of bad primes
    std::string ijson, itext;
    if (ispec.C && (o.integral || o.image)) {
        int md = o.maxdim >= 0 ? o.maxdim : std::max(0, ispec.C->max_dim());
        try {
            if (o.integral) {
                IntegralHomology H = integral_homology(*ispec.C, md);
                ijson += ",\n  \"integral_homology\": " + json_integral(H);
                itext += "integral homology of the whole complex:\n" + text_integral(H);
            }
            if (o.image) {
                IntegralHomology H = integral_image(*ispec.C, o.img_s, o.img_t, md);
                char b[160]; snprintf(b, sizeof b, "{\"s\": %s, \"t\": %s, \"groups\": ", fmt_num(o.img_s).c_str(), fmt_num(o.img_t).c_str());
                ijson += ",\n  \"integral_image\": " + std::string(b) + json_integral(H) + "}";
                itext += "image of H_*(K_" + fmt_num(o.img_s) + "; Z) -> H_*(K_" + fmt_num(o.img_t) + "; Z):\n" + text_integral(H);
            }
            if (o.integral && qmode_of(o) && cert.classified) {
                std::string fj;
                for (auto& bp : cert.bad) {
                    if (bp.is_big) continue;
                    Barcode b = run_field(be, bp.p, 0, be.maxdim(), nullptr).bars;
                    fj += std::string(fj.empty() ? "" : ", ") + "\"" + std::to_string(bp.p) + "\": " + json_bars(b);
                    itext += "barcode over F_" + std::to_string(bp.p) + " (a bad prime):\n";
                    for (size_t k = 0; k < b.size(); ++k) { itext += "  degree " + std::to_string(k) + ":"; for (auto& x : b[k]) itext += " [" + fmt_num(x.birth) + "," + (std::isinf(x.death) ? std::string("inf") : fmt_num(x.death)) + ")"; itext += "\n"; }
                }
                ijson += ",\n  \"bad_prime_barcodes\": {" + fj + "}";
            }
        } catch (std::exception& e) { itext += std::string("integral output: ") + e.what() + "\n"; ijson += std::string(",\n  \"integral_error\": \"") + e.what() + "\""; }
    }
    if (o.torsion) {
        if (!qmode) itext += "--torsion: needs Q coefficients\n";
        else if (!prep) itext += "--torsion: explicit complexes only (simplices, cells, sset, explicit cubical, sparse Rips)\n";
        else if (!cert.barcode_ok || !cert.candidates_complete()) itext += "--torsion: the candidate primes are not completely known\n";
        else {
            std::vector<std::pair<BigInt, bool>> pr = torsion::primes_to_follow(cert);
            int nt = o.threads > 0 ? o.threads : (int)std::max(1u, std::thread::hardware_concurrency());
            int md = std::min(o.maxdim >= 0 ? o.maxdim : prep->D, be.maxdim());
            torsion::TorsionResult T = torsion::trajectory(*prep, cert, md, nt, tcache, pr);
            std::string tj = json_torsion(T);
            tj.insert(1, std::string("\"proof\": ") + (cert.mode != "fast" && cert.certified() ? "true" : "false") + ", ");
            ijson += ",\n  \"torsion\": " + tj;
            if (cert.mode == "fast") itext += "  (--mode fast: the primes followed come from a probabilistic certificate; not a proof)\n";
            itext += text_torsion(T);
        }
    }
    output_barcode(o, main_bars, cert, qmode, ijson);
    if (qmode) {
        if (cert.mode == "fast")
            fprintf(stderr, "WARNING: --mode fast: the Q barcode and the bad primes below are probabilistic, NOT proved.\n"
                            "         Run without --mode fast (the default, --mode proved) for a proof.\n");
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
    if (!itext.empty()) fputs(itext.c_str(), stderr);
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
        else if (a == "--integral") o.integral = true;
        else if (a == "--torsion") o.torsion = true;
        else if (a == "--classify-by") { o.classify_by = val(); if (o.classify_by != "auto" && o.classify_by != "sweep" && o.classify_by != "rerun") { fprintf(stderr, "--classify-by: auto, sweep or rerun\n"); return 1; } }
        else if (a == "--image") { std::string v = val(); for (char& ch : v) if (ch == ',') ch = ' '; auto t = split_ws(v); if (t.size() != 2) { fprintf(stderr, "--image s,t\n"); return 1; } o.image = true; o.img_s = std::stod(t[0]); o.img_t = std::stod(t[1]); }
        else if (a == "--no-reduce") o.reduce = false;
        else if (a == "--lmax") o.lmax = std::stoi(val());
        else if (a == "--ring") o.ring = std::stoi(val());
        else if (a == "--no-homology") o.pc_homology = false;
        else if (a == "--all-products") o.all_products = true;
        else if (a == "--keep-pcgs") o.keep_pcgs = true;
        else if (a == "--directed") o.directed = true;
        else if (a == "--pairs") { std::string v = val(); for (char& ch : v) if (ch == ',') ch = ' '; auto t = split_ws(v); if (t.size() % 2) { fprintf(stderr, "--pairs s,t,...\n"); return 1; } for (size_t i = 0; i < t.size(); i += 2) o.pairs.push_back({std::stoi(t[i]), std::stoi(t[i + 1])}); }
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
        else if (f == "sset") {
            SimplicialSet S = read_simplicial_set(is);
            std::string e = S.validate(); if (!e.empty()) { fprintf(stderr, "invalid simplicial set: %s\n", e.c_str()); return 3; }
            in.complex = S.chains();
        }
        else if (f == "laurent") {
            // barcode over Q(t), bad irreducibles and bad primes; homology over Q[t, t^-1] and its bad primes
            laurent::LaurentComplex X = laurent::read_laurent(is);
            std::string e = X.validate(); if (!e.empty()) { fprintf(stderr, "invalid complex: %s\n", e.c_str()); return 3; }
            int md = o.maxdim >= 0 ? o.maxdim : std::max(0, X.max_dim());
            laurent::LaurentResult L = laurent::laurent_persistence(X, md, o.classify, true);
            std::string pr = o.print.empty() ? "text" : o.print;
            if (pr == "json") { fputs(laurent::to_json(L).c_str(), stdout); return 0; }
            auto bars = [](const laurent::LBarcode& B) {
                std::string s;
                for (size_t k = 0; k < B.size(); ++k) { s += "    degree " + std::to_string(k) + ":"; for (auto& b : B[k]) s += " [" + fmt_num(b.first) + "," + (std::isinf(b.second) ? std::string("inf") : fmt_num(b.second)) + ")"; s += "\n"; }
                return s;
            };
            auto hom = [](const laurent::LHomology& H, const char* R) {
                std::string s;
                for (size_t k = 0; k < H.free.size(); ++k) {
                    s += "    H_" + std::to_string(k) + " = "; bool any = false;
                    if (H.free[k]) { s += std::string(R) + (H.free[k] > 1 ? "^" + std::to_string(H.free[k]) : ""); any = true; }
                    for (auto& t : H.torsion[k]) { s += std::string(any ? " + " : "") + R + "/(" + laurent::to_string(t) + ")"; any = true; }
                    if (!any) s += "0";
                    s += "\n";
                }
                return s;
            };
            printf("barcode over Q(t):\n%s", bars(L.bars).c_str());
            if (o.classify) {
                printf("bad irreducibles (%zu candidates):%s\n", L.candidate_polys, L.bad_irreducibles.empty() ? " none" : "");
                for (auto& b : L.bad_irreducibles) printf("  %s%s: barcode over Q[t]/(f):\n%s", b.irreducible ? "f = " : "every irreducible factor f of ", laurent::to_string(b.g).c_str(), bars(b.bars).c_str());
                printf("bad primes over F_p(t) (%zu candidates):%s\n", L.candidate_primes.size(), L.bad_primes.empty() ? " none" : "");
                for (auto& kv : L.bad_primes) printf("  p = %llu:\n%s", (unsigned long long)kv.first, bars(kv.second).c_str());
            }
            printf("homology of the whole complex over Q[t, t^-1]:\n%s", hom(L.homology, "L").c_str());
            if (o.classify) {
                printf("primes where H_*( ; F_p[t, t^-1]) is not the reduction mod p (%zu candidates):%s\n", L.homology_candidate_primes.size(), L.homology_bad_primes.empty() ? " none" : "");
                for (auto& kv : L.homology_bad_primes) printf("  p = %llu:\n%s", (unsigned long long)kv.first, hom(kv.second, "L_p").c_str());
                if (!L.complete) { printf("INCOMPLETE: unfactored integers"); for (auto& u : L.unresolved) printf(" %s", u.c_str()); printf("\n"); }
            }
            return 0;
        }
        else if (f == "pcseries") {
            std::string text((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char>());
            pcg::PcSeriesData Dt = pcg::parse_pcseries(text);
            int md = o.maxdim >= 0 ? o.maxdim : 4;
            int p = 0;
            if (!o.coeff.empty() && o.coeff != "Q" && o.coeff != "q") p = (int)std::stoul(o.coeff);
            else {                                         // the prime of the groups, when they are p-groups of one prime
                std::set<int> ps; for (auto& G : Dt.groups) { int n = G.n; for (int q = 2; q <= n; ++q) if (n % q == 0) { ps.insert(q); while (n % q == 0) n /= q; } }
                if (ps.size() != 1) { fprintf(stderr, "pcseries: the groups are not p-groups of a single prime; give the field with --coeff p\n"); return 2; }
                p = *ps.begin();
            }
            pcg::PcSeriesResult R = pcg::pc_series(Dt.groups, Dt.maps, md, p, !o.keep_pcgs, Dt.alternatives);
            size_t m = Dt.groups.size();
            if (o.print == "json") {
                printf("{\n  \"p\": %d, \"groups\": [", p);
                for (size_t g = 0; g < m; ++g) printf("%s%d", g ? ", " : "", Dt.groups[g].n);
                printf("], \"maxdim\": %d, \"seconds\": %.4f,\n  \"dims\": [", md, R.seconds);
                for (size_t g = 0; g < m; ++g) { printf("%s[", g ? ", " : ""); for (int n = 0; n <= md; ++n) printf("%s%zu", n ? ", " : "", R.dims[g][n]); printf("]"); }
                printf("],\n  \"degrees\": [");
                for (int n = 0; n <= md; ++n) {
                    printf("%s\n    {\"degree\": %d, \"rank\": [", n ? "," : "", n);
                    for (size_t i = 0; i < m; ++i) { printf("%s[", i ? ", " : ""); for (size_t j = 0; j < m; ++j) printf("%s%zu", j ? ", " : "", R.rank[n][i][j]); printf("]"); }
                    printf("], \"barcode\": [");
                    for (size_t b = 0; b < R.barcode[n].size(); ++b) printf("%s[%d, %d, %zu]", b ? ", " : "", R.barcode[n][b].birth, R.barcode[n][b].death, R.barcode[n][b].mult);
                    printf("], \"maps\": [");
                    for (size_t i = 0; i + 1 < m; ++i) {
                        printf("%s[", i ? ", " : "");
                        const auto& M = R.maps[i][n];
                        for (size_t t = 0; t < M.size(); ++t) { printf("%s[", t ? ", " : ""); for (size_t s2 = 0; s2 < M[t].size(); ++s2) printf("%s%u", s2 ? ", " : "", M[t][s2]); printf("]"); }
                        printf("]");
                    }
                    printf("]}");
                }
                printf("\n  ]\n}\n");
                return 0;
            }
            printf("chain of %zu groups (orders", m);
            for (auto& G : Dt.groups) printf(" %d", G.n);
            printf("), homology with F_%d coefficients through degree %d; %.3f s\n", p, md, R.seconds);
            for (int n = 0; n <= md; ++n) {
                printf("degree %d: dims", n); for (size_t g = 0; g < m; ++g) printf(" %zu", R.dims[g][n]);
                printf("\n  persistence matrix (rank of H_%d(G_i) -> H_%d(G_j), i <= j):\n", n, n);
                for (size_t i = 0; i < m; ++i) { printf("   "); for (size_t j = 0; j < m; ++j) { if (j < i) printf("   ."); else printf(" %3zu", R.rank[n][i][j]); } printf("\n"); }
                printf("  barcode (groups numbered from 1; [i, j] x k: k classes born in G_i, alive through G_j):");
                for (auto& b : R.barcode[n]) printf(" [%d, %d]x%zu", b.birth + 1, b.death + 1, b.mult);
                printf("\n");
            }
            printf("  (the matrices of the maps: --print json)\n");
            return 0;
        }
        else if (f == "sylow") {
            std::string text((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char>());
            pcg::SylowData Dt = pcg::parse_sylow(text);
            int md = o.maxdim >= 0 ? o.maxdim : 4;
            pcg::SylowResult R = pcg::sylow_compute(Dt, md, o.pc_homology, o.ring, o.all_products, !o.keep_pcgs);
            auto inv_str = [](size_t fr, const std::vector<uint64_t>& tors) {
                std::string s; if (fr) s = fr > 1 ? "Z^" + std::to_string(fr) : "Z";
                for (auto t : tors) { if (!s.empty()) s += " + "; s += "Z/" + std::to_string(t); }
                return s.empty() ? std::string("0") : s;
            };
            auto vec = [](const std::vector<uint32_t>& x) { std::string s = "["; for (size_t i = 0; i < x.size(); ++i) s += (i ? ", " : "") + std::to_string(x[i]); return s + "]"; };
            if (o.print == "json") {
                printf("{\n  \"order\": %lld, \"primes\": [", R.order);
                for (size_t i = 0; i < R.primes.size(); ++i) printf("%s%d", i ? ", " : "", R.primes[i]);
                printf("], \"complete\": %s, \"maxdim\": %d, \"seconds\": %.4f", R.complete ? "true" : "false", md, R.seconds);
                if (R.has_homology) {
                    printf(",\n  \"homology\": [");
                    for (int d = 0; d <= md; ++d) {
                        printf("%s\n    {\"degree\": %d, \"free\": %zu, \"torsion\": [", d ? "," : "", d, R.H.free[d]);
                        for (size_t i = 0; i < R.H.torsion[d].size(); ++i) printf("%s%llu", i ? ", " : "", (unsigned long long)R.H.torsion[d][i]);
                        printf("]}");
                    }
                    printf("\n  ]");
                }
                if (R.has_ring) {
                    const auto& C = R.R;
                    printf(",\n  \"cohomology\": {\"p\": %d, \"dims\": [", C.p);
                    for (size_t i = 0; i < C.dims.size(); ++i) printf("%s%zu", i ? ", " : "", C.dims[i]);
                    printf("], \"generators_by_degree\": [");
                    for (size_t i = 0; i < C.indecomposables.size(); ++i) printf("%s%zu", i ? ", " : "", C.indecomposables[i]);
                    printf("], \"sylow_dims\": [");
                    for (size_t i = 0; i < C.sylow_dims.size(); ++i) printf("%s%zu", i ? ", " : "", C.sylow_dims[i]);
                    printf("],\n    \"generators\": [");
                    for (size_t i = 0; i < C.generators.size(); ++i) printf("%s[%d, %d]", i ? ", " : "", C.generators[i].first, C.generators[i].second);
                    printf("],\n    \"basis\": [");
                    for (size_t n = 0; n < C.basis.size(); ++n) { printf("%s[", n ? ", " : ""); for (size_t i = 0; i < C.basis[n].size(); ++i) printf("%s%s", i ? ", " : "", vec(C.basis[n][i]).c_str()); printf("]"); }
                    printf("],\n    \"products_of\": \"%s\", \"products\": [", C.all_products ? "all basis classes" : "generators");
                    bool first = true;
                    for (auto& kv : C.products) { printf("%s\n      [%d, %d, %d, %d, %s]", first ? "" : ",", std::get<0>(kv.first), std::get<1>(kv.first), std::get<2>(kv.first), std::get<3>(kv.first), vec(kv.second).c_str()); first = false; }
                    printf("\n    ]}");
                }
                printf("\n}\n");
                return 0;
            }
            printf("group of order %lld, through its Sylow subgroups for p =", R.order);
            for (int q : R.primes) printf(" %d", q);
            printf("%s; %.3f s\n", R.complete ? "" : " (not every prime dividing |G|: homology is only the p-parts for these primes)", R.seconds);
            if (R.has_homology) {
                printf("integral homology%s:\n", R.complete ? "" : " (p-parts)");
                for (int d = 0; d <= md; ++d) printf("  H_%d(G; Z) = %s\n", d, inv_str(R.H.free[d], R.H.torsion[d]).c_str());
            }
            if (R.has_ring) {
                printf("mod-%d cohomology ring through degree %d (stable elements in H^*(P; F_%d)):\n  dim H^n(G; F_%d):", R.R.p, md, R.R.p, R.R.p);
                for (auto x : R.R.dims) printf(" %zu", x);
                printf("\n  minimal generators by degree (n = 1..%d):", md);
                for (int d = 1; d <= md; ++d) printf(" %zu", R.R.indecomposables[d]);
                printf("\n  (generators, basis and %s: --print json)\n", R.R.all_products ? "all products of basis classes" : "products of the generators with every class");
            }
            return 0;
        }
        else if (f == "pcgroup") {
            // n k / pcgs element ids / relative orders / identity / n x n multiplication table ('#' lines are comments)
            std::vector<long long> tok; std::string line;
            while (std::getline(is, line)) {
                size_t h = line.find('#'); if (h != std::string::npos) line.resize(h);
                for (char& ch : line) if (ch == ',') ch = ' ';
                for (auto& t : split_ws(line)) tok.push_back(std::stoll(t));
            }
            size_t at = 0;
            auto next = [&]() { if (at >= tok.size()) throw std::runtime_error("pcgroup: the file ends early"); return tok[at++]; };
            int n = (int)next(), k = (int)next();
            if (n < 1 || k < 0 || n > (1 << 16)) throw std::runtime_error("pcgroup: bad order or pcgs length");
            std::vector<int> pcgs(k), rel(k);
            for (auto& x : pcgs) x = (int)next();
            for (auto& x : rel) x = (int)next();
            int e = (int)next();
            std::vector<std::vector<int>> T(n, std::vector<int>(n));
            for (auto& r : T) for (auto& x : r) x = (int)next();
            pcg::PcgsList alts;                       // optional alternative pcgs after the table: k / elements / relative orders
            while (at < tok.size()) {
                int ka = (int)next(); if (ka < 0 || ka > 64) throw std::runtime_error("pcgroup: bad alternative pcgs length");
                std::vector<int> pa(ka), ra(ka); for (auto& x : pa) x = (int)next(); for (auto& x : ra) x = (int)next();
                alts.push_back({pa, ra});
            }
            pcg::PcGroup G(T, pcgs, rel, e);
            int md = o.maxdim >= 0 ? o.maxdim : 4;
            pcg::PcResult R = pcg::pc_group_compute(G, md, o.pc_homology, o.ring, o.all_products, !o.keep_pcgs, alts);
            auto group_str = [](size_t fr, const std::vector<uint64_t>& tors) {
                std::string s; if (fr) s = fr > 1 ? "Z^" + std::to_string(fr) : "Z";
                std::map<uint64_t, size_t> c; for (auto t : tors) ++c[t];
                for (auto& kv : c) { if (!s.empty()) s += " + "; s += "(Z/" + std::to_string(kv.first) + ")" + (kv.second > 1 ? "^" + std::to_string(kv.second) : ""); }
                return s.empty() ? std::string("0") : s;
            };
            if (o.print == "json") {
                printf("{\n  \"order\": %d, \"relative_orders\": [", n);
                for (int i = 0; i < k; ++i) printf("%s%d", i ? ", " : "", G.rel[i]);
                printf("], \"maxdim\": %d, \"seconds\": %.4f, \"states\": %zu,", md, R.seconds, R.states);
                auto used = [](const char* what, const pcg::PcUsed& u) {
                    printf("\n  \"pcgs_%s\": {\"choice\": \"%s\", \"pcgs\": [", what, u.choice.c_str());
                    for (size_t i = 0; i < u.pcgs.size(); ++i) printf("%s%d", i ? ", " : "", u.pcgs[i]);
                    printf("], \"relative_orders\": ["); for (size_t i = 0; i < u.rel.size(); ++i) printf("%s%d", i ? ", " : "", u.rel[i]); printf("]},");
                };
                if (R.has_homology) used("homology", R.used_homology);
                if (R.has_ring) used("ring", R.used_ring);
                printf("\n  \"chains\": [");
                for (size_t i = 0; i < R.chains.size(); ++i) printf("%s%zu", i ? ", " : "", R.chains[i]);
                printf("]");
                if (R.has_homology) {
                    printf(",\n  \"homology\": [");
                    for (int d = 0; d <= md; ++d) {
                        printf("%s\n    {\"degree\": %d, \"free\": %zu, \"torsion\": [", d ? "," : "", d, R.H.free[d]);
                        for (size_t i = 0; i < R.H.torsion[d].size(); ++i) printf("%s%llu", i ? ", " : "", (unsigned long long)R.H.torsion[d][i]);
                        printf("]}");
                    }
                    printf("\n  ]");
                }
                if (R.has_ring) {
                    printf(",\n  \"cohomology\": {\"p\": %d, \"dims\": [", R.R.p);
                    for (size_t i = 0; i < R.R.dims.size(); ++i) printf("%s%zu", i ? ", " : "", R.R.dims[i]);
                    printf("], \"generators_by_degree\": [");
                    for (size_t i = 0; i < R.R.indecomposables.size(); ++i) printf("%s%zu", i ? ", " : "", R.R.indecomposables[i]);
                    printf("], \"generators\": [");
                    for (size_t i = 0; i < R.R.generators.size(); ++i) printf("%s[%d, %d]", i ? ", " : "", R.R.generators[i].first, R.R.generators[i].second);
                    printf("],\n    \"products_of\": \"%s\", \"products\": [", R.R.all_products ? "all basis classes" : "generators");
                    bool first = true;
                    for (auto& kv : R.R.products) {
                        printf("%s\n      [%d, %d, %d, %d, [", first ? "" : ",", std::get<0>(kv.first), std::get<1>(kv.first), std::get<2>(kv.first), std::get<3>(kv.first));
                        for (size_t i = 0; i < kv.second.size(); ++i) printf("%s%u", i ? ", " : "", kv.second[i]);
                        printf("]]"); first = false;
                    }
                    printf("\n    ]}");
                }
                printf("\n}\n");
                return 0;
            }
            printf("group of order %d, pcgs of length %d (relative orders", n, k);
            for (int i = 0; i < k; ++i) printf(" %d", G.rel[i]);
            printf("); %.3f s\n", R.seconds);
            if (R.has_homology) {
                printf("integral homology:\n");
                for (int d = 0; d <= md; ++d) printf("  H_%d(G; Z) = %s\n", d, group_str(R.H.free[d], R.H.torsion[d]).c_str());
            }
            if (R.has_ring) {
                printf("mod-%d cohomology ring through degree %d:\n  dim H^n(G; F_%d):", R.R.p, md, R.R.p);
                for (auto x : R.R.dims) printf(" %zu", x);
                printf("\n  minimal generators by degree (n = 1..%d):", md);
                for (int d = 1; d <= md; ++d) printf(" %zu", R.R.indecomposables[d]);
                printf("\n  (%s: --print json; they are cup products)\n", R.R.all_products ? "all products of basis classes" : "products of the generators with every basis class");
            }
            if (o.stats) { printf("  [Anick chains per degree:"); for (auto c : R.chains) printf(" %zu", c); printf("; %zu homotopy states]\n", R.states); }
            return 0;
        }
        else if (f == "magnitude") {
            // edge list: "u v" per line (0-based vertices); MH_{k,l} for l <= lmax with torsion
            std::vector<std::pair<int, int>> E; std::string line; int n = 0;
            while (std::getline(is, line)) { if (line.empty() || line[0] == '#') continue; for (char& ch : line) if (ch == ',') ch = ' '; auto t = split_ws(line); if (t.size() < 2) continue; int a = std::stoi(t[0]), b = std::stoi(t[1]); E.push_back({a, b}); n = std::max(n, std::max(a, b) + 1); }
            MagGraph G(n, E, o.directed);
            MagResult R = magnitude_homology(G, o.lmax, o.threads > 0 ? o.threads : (int)std::max(1u, std::thread::hardware_concurrency()), o.pairs);
            if (o.print == "json") {
                printf("{\n  \"vertices\": %d, \"edges\": %zu, \"directed\": %s, \"pairs\": %zu, \"seconds\": %.4f,\n  \"groups\": [", n, E.size(), o.directed ? "true" : "false", R.pairs, R.seconds);
                bool first = true;
                for (int l = 0; l <= o.lmax; ++l) for (int k = 0; k <= l; ++k) {
                    const MagEntry& e = R.table[l][k];
                    if (!e.generators) continue;
                    printf("%s\n    {\"k\": %d, \"l\": %d, \"free\": %zu, \"torsion\": {", first ? "" : ",", k, l, e.free);
                    bool f2 = true; for (auto& kv : e.torsion) { printf("%s\"%s\": %zu", f2 ? "" : ", ", kv.first.c_str(), kv.second); f2 = false; }
                    printf("}, \"generators\": %zu, \"after_reduction\": %zu}", e.generators, e.residual);
                    first = false;
                }
                printf("\n  ]\n}\n");
                return 0;
            }
            printf("magnitude homology MH_{k,l} over Z (%d vertices, %zu edges, %s; %zu (s, t) pairs; %.3f s)\n", n, E.size(), o.directed ? "directed" : "undirected", R.pairs, R.seconds);
            for (int l = 0; l <= o.lmax; ++l) for (int k = 0; k <= l; ++k) {
                const MagEntry& e = R.table[l][k];
                if (!e.free && e.torsion.empty() && !o.stats) continue;
                printf("  MH_{%d,%d} = ", k, l);
                bool any = false;
                if (e.free) { printf("Z%s", e.free > 1 ? ("^" + std::to_string(e.free)).c_str() : ""); any = true; }
                for (auto& kv : e.torsion) { printf("%s(Z/%s)%s", any ? " + " : "", kv.first.c_str(), kv.second > 1 ? ("^" + std::to_string(kv.second)).c_str() : ""); any = true; }
                if (!any) printf("0");
                printf("    [%zu generators, %zu after reduction]\n", e.generators, e.residual);
            }
            return 0;
        }
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
            if (use_rips_components(E, o.coeff, o.mode, md, o.fp, o.zero, o.torsion)) {
                if (auto parts = split_rips(R, E, md, o.threads, o.collapse, o.stats)) {
                    auto options = o;
                    options.threads = 1;  // the worker budget is already assigned to components
                    return solve(*parts, options);
                }
            }
            std::vector<WEdge>().swap(E);
            if (sparse_src) R.densify();
            // edge collapse pays for itself in degrees >= 2; for H0/H1 the active H1 engine is cheaper
            collapse_implicit_input(R, md, o.collapse, o.stats);
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
        bool chain_input = (f == "cells" || f == "sset");
        std::string err = chain_input ? validate_chain_complex(in.complex) : in.complex.validate(false);
        if (!err.empty()) { fprintf(stderr, "invalid complex: %s\n", err.c_str()); return 3; }
        lap("validate");
        // Integral output needs the residual complex small enough for dense Smith forms, so the
        // reduction also runs for any explicit complex when integral, image or torsion output is asked for.
        if (o.reduce && (chain_input || o.integral || o.image || o.torsion)) {
            // algebraic Morse reduction (unit pairs with equal births): a filtered chain homotopy equivalence over Z
            ReducedComplex R = algebraic_morse(in.complex);
            if (o.stats) fprintf(stderr, "  [algebraic Morse reduction: %zu -> %zu cells (%zu pairs)]\n", in.complex.size(), R.C.size(), R.pairs);
            int md = o.maxdim >= 0 ? o.maxdim : in.complex.max_dim();
            Prepared P(R.C);
            lap("reduce and prepare");
            ReducedBackend be(R, P, md, choose_direction(P, o.direction, md)); be.keep_zero_length = o.zero;
            IntegralSpec ispec; ispec.C = &R.C;
            return solve(be, o, ispec);
        }
        Prepared P(in.complex);
        lap("prepare");
        ExplicitBackend be(P, o.maxdim >= 0 ? o.maxdim : P.D, choose_direction(P, o.direction, o.maxdim >= 0 ? o.maxdim : P.D)); be.keep_zero_length = o.zero;
        IntegralSpec ispec; ispec.C = &in.complex;
        return solve(be, o, ispec);
    } catch (std::exception& e) {
        fprintf(stderr, "error: %s\n", e.what()); return 2;
    }
}
