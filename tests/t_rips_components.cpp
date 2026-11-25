// Force rare arithmetic failures in just one worker, using the production rings.
#include "perseus/rips_components.hpp"
#include <cassert>

using namespace perseus;

struct ArithmeticProbe : RipsBackend {
    bool inject;
    ArithmeticProbe(const RipsInput& input, int md) : RipsBackend(input, md), inject(input.n == 4) {}

    template<class Ring>
    BRun<Clear> run(const Ring& ring, int lo, int hi, const Clear* in) const {
        auto result = RipsBackend::run(ring, lo, hi, in);
        if (inject && hi == 1) {
            if constexpr (std::is_same_v<Ring, Q64>) {
                ring.mul(ring.from_int(Q64::LIM), ring.from_int(2));
                result.bars[1].clear();  // the caller must discard this failed reduction
            } else if constexpr (std::is_same_v<Ring, MontN<2>>) {
                ring.inv(ring.from_int(3));
            }
        }
        return result;
    }
};

struct ProbedComponents : RipsComponents {
    template<class Ring>
    BRun<Clear> run(const Ring& ring, int lo, int hi, const Clear* in) const {
        return RipsComponents::run<Ring, ArithmeticProbe>(ring, lo, hi, in);
    }
};

int main() {
    RipsInput input;
    input.n = 8;
    std::vector<WEdge> edges{{0,1,1}, {1,2,1}, {2,3,1}, {3,0,1}, {4,5,1}, {5,6,1}};
    auto split = split_rips(input, edges, 1, 8, "off");
    assert(split && split->workers == 2 && split->candidate_threads == 4);
    const ProbedComponents be{std::move(*split)};
    std::string arithmetic;
    const auto exact = run_exact(be, 1, arithmetic);
    assert(arithmetic == "int64, big rationals in degree 1");
    assert(exact.bars[1].size() == 1 && exact.bars[1][0].birth == 1 && std::isinf(exact.bars[1][0].death));

    const BigInt modulus = BigInt::from_u64((1ULL << 62) + 1) * BigInt(3);
    BigInt factor;
    BRun<RipsComponents::Clear> result;
    const bool ok = run_mod_big(be, modulus, 1, 1, &exact.clear_out[0], result, factor);
    assert(!ok && factor == BigInt(3));

    std::vector<std::future<Barcode>> jobs;
    for (uint64_t p : {2, 3, 5, 7}) {
        jobs.push_back(std::async(std::launch::async, [&, p] {
            return run_field(be, p, 1, 1, &exact.clear_out[0]).bars;
        }));
    }
    for (auto& job : jobs) assert(job.get()[1] == exact.bars[1]);

    try { run_lanes(be, {3, 5}, 0, 1, nullptr); }
    catch (const std::logic_error&) { return 0; }
    return 1;
}
