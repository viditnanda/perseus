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

inline void collapse_implicit_input(RipsInput& R, int md, const std::string& mode, bool stats = false,
                                    std::string* log = nullptr) {
    const bool automatic = mode == "auto";
    if (mode != "full" && mode != "on" && !(automatic && md >= 2 && R.n >= 32)) return;
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
}

inline std::optional<uint64_t> component_prime(const std::vector<WEdge>& edges, const std::string& coeff,
                                               int md, const std::vector<uint64_t>& extra_fields,
                                               bool zero, bool torsion) {
    if (edges.empty() || md < 0 || zero || torsion || !extra_fields.empty() || coeff.empty() ||
        coeff == "Q" || coeff == "q" || coeff == "0") return std::nullopt;
    uint64_t prime = std::stoull(coeff);
    if (!is_prime64(prime)) return std::nullopt;  // normal field path reports the invalid coefficient
    return prime;
}

// Called only when the automatic sparse-engine probe has already built the radius graph.
inline std::optional<Barcode> component_field(RipsInput& input, std::vector<WEdge>& edges,
                                               int md, uint64_t prime, const std::string& collapse,
                                               int threads, bool stats = false) {
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
    const double inf = std::numeric_limits<double>::infinity();
    for (const auto& edge : edges) {
        if (edge.a == edge.b || !std::isfinite(edge.w)) continue;
        int a = local[edge.a], b = local[edge.b];
        parts[part_id[root(edge.a)]].edges.emplace_back(a, b, edge.w);
    }

    std::vector<WEdge>().swap(edges);
    std::stable_sort(parts.begin(), parts.end(), [](const RipsInput& a, const RipsInput& b) { return a.n > b.n; });
    std::vector<std::string> logs(parts.size());
    std::vector<Barcode> results(parts.size());
    std::atomic<size_t> next{0};
    auto worker = [&] {
        for (size_t i; (i = next.fetch_add(1)) < parts.size();) {
            auto& part = parts[i];
            part.densify();
            decltype(part.edges)().swap(part.edges);
            collapse_implicit_input(part, md, collapse, stats, &logs[i]);
            RipsBackend backend(part, md);
            std::vector<double>().swap(part.lower);
            results[i] = run_field(backend, prime, 0, md, nullptr).bars;
        }
    };
    if (threads <= 0) threads = std::max(1u, std::thread::hardware_concurrency());
    int workers = std::min<int>(threads, parts.size());
    if (stats) fprintf(stderr, "  [Rips: %zu nontrivial components, %d isolated vertices, %d workers]\n", parts.size(), isolated, workers);
    std::vector<std::future<void>> jobs;
    for (int i = 1; i < workers; ++i) jobs.push_back(std::async(std::launch::async, worker));
    worker();
    for (auto& job : jobs) job.get();

    for (const auto& log : logs) if (!log.empty()) fputs(log.c_str(), stderr);

    Barcode bars(md + 1);
    bars[0].assign(isolated, {0.0, inf});
    for (const auto& result : results)
        for (size_t k = 0; k < result.size(); ++k)
            bars[k].insert(bars[k].end(), result[k].begin(), result[k].end());
    for (auto& degree : bars) std::sort(degree.begin(), degree.end());
    return bars;
}
} // namespace perseus
