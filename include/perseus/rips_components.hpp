// Independent components of the final radius graph give a direct sum at every filtration value.
#pragma once
#include "certify.hpp"
#include "collapse.hpp"
#include "rips_engine.hpp"
#include "rips_sparse.hpp"
#include <atomic>
#include <future>
#include <numeric>
#include <optional>

namespace perseus {

inline void collapse_implicit_input(RipsInput& R, int md, const std::string& mode, bool stats = false) {
    const bool automatic = mode == "auto";
    if (mode != "full" && mode != "on" && !(automatic && md >= 2 && R.n >= 32)) return;
    CollapseStats cs;
    RipsInput collapsed = edge_collapse(R, mode == "full" ? 1000 : 32, cs,
                                        mode == "full" ? 0.0 : 0.02, automatic ? 0.10 : 0.0);
    if (stats) {
        if (cs.used) fprintf(stderr, "  [edge collapse: %zu -> %zu edges (%zu shifted), %d passes, %.3f s]\n",
                             cs.edges_before, cs.edges_after, cs.shifted, cs.passes, cs.seconds);
        else fprintf(stderr, "  [edge collapse not worth it here, abandoned after %.3f s]\n", cs.seconds);
    }
    if (cs.used) R.lower.swap(collapsed.lower);
}

// Called only when the automatic sparse-engine probe has already built the radius graph.
inline std::optional<Barcode> component_field(int n, const std::vector<WEdge>& edges, double threshold,
                                               int md, uint64_t prime, const std::string& collapse,
                                               int threads, bool stats = false) {
    if (edges.empty()) return std::nullopt;
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
        }
        local[v] = parts[part_id[r]].n++;
    }
    const double inf = std::numeric_limits<double>::infinity();
    for (auto& part : parts) part.lower.assign(size_t(part.n) * (part.n - 1) / 2, inf);
    for (const auto& edge : edges) {
        if (edge.a == edge.b || !std::isfinite(edge.w)) continue;
        int a = local[edge.a], b = local[edge.b];
        if (a < b) std::swap(a, b);
        parts[part_id[root(edge.a)]].lower[size_t(a) * (a - 1) / 2 + b] = edge.w;
    }

    std::vector<Barcode> results(parts.size());
    std::atomic<size_t> next{0};
    auto worker = [&] {
        for (size_t i; (i = next.fetch_add(1)) < parts.size();) {
            auto& part = parts[i];
            collapse_implicit_input(part, md, collapse, stats);
            RipsBackend backend(part, md);
            std::vector<double>().swap(part.lower);
            results[i] = run_field(backend, prime, 0, md, nullptr).bars;
        }
    };
    if (threads <= 0) threads = std::max(1u, std::thread::hardware_concurrency());
    int workers = std::min<int>(threads, parts.size());
    if (stats) fprintf(stderr, "  [Rips: %d components, %d isolated vertices, %d workers]\n", count, isolated, workers);
    std::vector<std::future<void>> jobs;
    for (int i = 1; i < workers; ++i) jobs.push_back(std::async(std::launch::async, worker));
    worker();
    for (auto& job : jobs) job.get();

    Barcode bars(md + 1);
    bars[0].assign(isolated, {0.0, inf});
    for (const auto& result : results)
        for (size_t k = 0; k < result.size(); ++k)
            bars[k].insert(bars[k].end(), result[k].begin(), result[k].end());
    for (auto& degree : bars) std::sort(degree.begin(), degree.end());
    return bars;
}
} // namespace perseus
