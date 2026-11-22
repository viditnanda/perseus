// Independent components of the final radius graph give a direct sum at every filtration value.
#pragma once
#include "certify.hpp"
#include "collapse.hpp"
#include "rips_engine.hpp"
#include "rips_sparse.hpp"
#include <atomic>
#include <future>
#include <optional>

namespace perseus {

inline bool collapse_implicit_input(RipsInput& R, int md, const std::string& mode, bool stats = false,
                                    std::string* log = nullptr) {
    const bool automatic = mode == "auto";
    if (mode != "full" && mode != "on" && !(automatic && md >= 2 && R.n >= 32)) return false;
    CollapseStats cs;
    RipsInput collapsed = edge_collapse(R, mode == "full" ? 1000 : 32, cs,
                                        mode == "full" ? 0.0 : 0.02, automatic ? 0.10 : 0.0, log);
    if (stats) {
        char line[192];
        if (cs.used) snprintf(line, sizeof(line), "  [edge collapse: %zu -> %zu edges (%zu shifted), %d passes, %.3f s]\n",
                             cs.edges_before, cs.edges_after, cs.shifted, cs.passes, cs.seconds);
        else snprintf(line, sizeof(line), "  [edge collapse not worth it here, abandoned after %.3f s]\n", cs.seconds);
        if (log) log->append(line);
        else fputs(line, stderr);
    }
    if (cs.used) R.lower.swap(collapsed.lower);
    return cs.used;
}

// The automatic sparse-engine probe has already built the radius graph.
inline bool use_rips_components(const std::vector<WEdge>& edges, const std::string& coeff,
                                const std::string& mode, int md, const std::vector<uint64_t>& extra_fields,
                                bool zero, bool torsion) {
    if (edges.empty() || md < 0 || zero || torsion || !extra_fields.empty()) return false;
    if (coeff.empty() || coeff == "Q" || coeff == "q" || coeff == "0") return mode != "fast";
    return is_prime64(std::stoull(coeff));
}

struct RipsComponents {
    mutable std::vector<RipsInput> parts;
    int isolated = 0;
    int md, threads;
    std::string collapse;
    bool stats;
    bool keep_zero_length = false;
    mutable bool prepared = false;
    int maxdim() const { return md; }

    // Each run owns the worker pool. Between runs retain only the sparse graph,
    // with the collapsed weights, so local clearing indices remain valid.
    using Clear = std::vector<RipsBackend::Clear>;
    template<class Ring>
    BRun<Clear> run(const Ring& ring, int kmin, int kmax, const Clear* in) const;
};

inline std::optional<RipsComponents> split_rips(RipsInput& input, std::vector<WEdge>& edges,
                                                int md, int threads, std::string collapse, bool stats = false) {
    const int n = input.n;
    std::vector<int> parent(n, -1);
    auto root = [&](int v) {
        while (parent[v] >= 0) {
            if (parent[parent[v]] >= 0) parent[v] = parent[parent[v]];
            v = parent[v];
        }
        return v;
    };
    int count = n;
    for (const auto& edge : edges) {
        if (!std::isfinite(edge.w)) continue;
        int a = root(edge.a), b = root(edge.b);
        if (a == b) continue;
        if (parent[a] > parent[b]) std::swap(a, b);
        parent[a] += parent[b];
        parent[b] = a;
        if (--count == 1) return std::nullopt;
    }

    // All distances needed below are in edges. The connected fallback above still owns its input.
    const double threshold = input.threshold;
    input = RipsInput();
    std::vector<RipsInput> parts;
    std::vector<int> part_id(n, -1), local(n);
    int isolated = 0;
    for (int v = 0; v < n; ++v) {
        int r = root(v);
        if (parent[r] == -1) { ++isolated; continue; }
        if (part_id[r] == -1) {
            part_id[r] = static_cast<int>(parts.size());
            parts.emplace_back();
            parts.back().threshold = threshold;
            parts.back().from_edges = true;
        }
        local[v] = parts[part_id[r]].n++;
    }
    for (const auto& edge : edges) {
        if (edge.a == edge.b || !std::isfinite(edge.w)) continue;
        int a = local[edge.a], b = local[edge.b];
        parts[part_id[root(edge.a)]].edges.emplace_back(a, b, edge.w);
    }

    std::vector<WEdge>().swap(edges);
    std::stable_sort(parts.begin(), parts.end(), [](const RipsInput& a, const RipsInput& b) { return a.n > b.n; });
    if (threads <= 0) threads = std::max(1u, std::thread::hardware_concurrency());
    int workers = std::min<int>(threads, parts.size());
    if (stats) fprintf(stderr, "  [Rips: %zu nontrivial components, %d isolated vertices, %d workers]\n", parts.size(), isolated, workers);
    return RipsComponents{std::move(parts), isolated, md, workers, std::move(collapse), stats};
}

// Rings with a recoverable arithmetic failure report it to the caller after the
// independent runs. The existing certifier retries with big rationals or a factor.
template<class R> inline void collect_ring_status(const R&, const R&) {}
inline void collect_ring_status(const Q64& dst, const Q64& src) { dst.overflow |= src.overflow; }
template<int L> void collect_ring_status(const MontN<L>& dst, const MontN<L>& src) {
    if (src.split_found && !dst.split_found) { dst.split = src.split; dst.split_found = true; }
}

template<class Ring>
BRun<RipsComponents::Clear> RipsComponents::run(const Ring& ring, int kmin, int kmax, const Clear* in) const {
    std::vector<std::string> logs(parts.size());
    std::vector<BRun<RipsBackend::Clear>> results(parts.size());
    std::vector<Ring> rings(parts.size(), ring);
    std::atomic<size_t> next{0};
    auto worker = [&] {
        for (size_t i; (i = next.fetch_add(1)) < parts.size();) {
            auto& part = parts[i];
            part.densify();
            if (!prepared && collapse_implicit_input(part, md, collapse, stats, &logs[i])) {
                // Collapse only removes or shifts existing edges. Save their new
                // weights for candidate reruns without retaining a distance table.
                for (auto& [a, b, weight] : part.edges) weight = part.dist(a, b);
                auto removed = [](const auto& edge) { return !std::isfinite(std::get<2>(edge)); };
                part.edges.erase(std::remove_if(part.edges.begin(), part.edges.end(), removed), part.edges.end());
            }
            RipsBackend backend(part, md);
            std::vector<double>().swap(part.lower);
            results[i] = backend.run(rings[i], kmin, kmax, in ? &(*in)[i] : nullptr);
        }
    };
    std::vector<std::future<void>> jobs;
    for (int i = 1; i < threads; ++i) jobs.push_back(std::async(std::launch::async, worker));
    worker();
    for (auto& job : jobs) job.get();
    prepared = true;

    for (const auto& log : logs) if (!log.empty()) fputs(log.c_str(), stderr);

    BRun<Clear> out;
    out.kmin = kmin;
    out.bars.resize(kmax + 1);
    out.clear_out.resize(kmax + 1);
    if (kmin == 0) out.bars[0].assign(isolated, {0.0, std::numeric_limits<double>::infinity()});
    for (size_t i = 0; i < results.size(); ++i) {
        auto& r = results[i];
        collect_ring_status(ring, rings[i]);
        for (int k = kmin; k <= kmax; ++k) {
            out.bars[k].insert(out.bars[k].end(), r.bars[k].begin(), r.bars[k].end());
            out.clear_out[k].push_back(std::move(r.clear_out[k]));
        }
        out.nonunit.insert(out.nonunit.end(), std::make_move_iterator(r.nonunit.begin()),
                           std::make_move_iterator(r.nonunit.end()));
    }
    for (auto& degree : out.bars) std::sort(degree.begin(), degree.end());
    return out;
}
} // namespace perseus
