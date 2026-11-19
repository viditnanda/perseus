# Performance roadmap: status, decisions and protocol

This file is for anyone picking up Perseus 6 performance work, whether a person or a Claude instance. It records what is in flight, what was agreed with Vidit Nanda (the maintainer), how to measure, and what has already been tried. Update it in the same pull request as the work it describes.

Last updated: 19 Nov 2026 (merge order and plan agreed with Vidit).

## Ground rules (from Vidit)

* **Do not merge, tag, release or upload to PyPI.** Vidit reviews every pull request. Work on a branch and open a pull request against `main`.
* Commit messages end with the `Co-Authored-By` / `Claude-Session` lines of the session doing the work. Pull-request descriptions end with the "Generated with Claude Code" line and the session link.
* **The Python build does not notice header changes.** Delete `python/build` before every `pip install`. Install into a private target such as `pip install --no-deps --no-build-isolation --target /tmp/pyX ./python`, never into site-packages, so that parallel sessions do not overwrite each other's modules.
* Every speed claim needs identical output. Compare barcodes, and for Q runs the certificate (bad primes and their intervals), against the previous binary or another engine.
* State limitations plainly in the pull request: what was not measured, and what was tried and did not help.

## Measurement protocol

* The test machine is a 2-core VM with about 6 GB usable (the command shell's cgroup limit). Timings vary by about ±10% between runs, so take the best of three where runs are short, and do not claim gains below about 10% from single runs.
* While benchmarking, create `/tmp/perseus_bench.lock`. Other sessions must not start heavy jobs while it exists; pause theirs with SIGSTOP if needed, and resume with SIGCONT. Remove the lock afterwards.
* Python tools are timed in-process around the call, after their imports. GUDHI's sklearn wrapper alone takes about a second to import, and an older script timed it by mistake.
* Memory is peak RSS read right after the call, before any barcode comparison. Building Python `Counter`s of millions of intervals used to dominate the measured memory.
* Command-line runs are timed separately: they include process start and reading the file.
* Scripts: `bench/cubical2d.py` (2D images), `bench/cubical.py` (3D/4D grids; `--others-only` re-measures only GUDHI and CubicalRipser), and `bench/vs_*.py` for the rest. Inputs are regenerated deterministically into `/tmp/cubical_bench`.
* For phase timings, set `PERSEUS_DEBUG=1`. The engines then print per-degree reduction times, additions, merged entries and the longest working column to stderr.

## Open pull requests and merge order

Vidit's merge order (19 Nov 2026): **#4, then #2, then #3**. He verifies each one himself before merging, and #3, the largest, gets the most scrutiny.

| PR | branch | content | state |
|---|---|---|---|
| #4 | `cubical-bittree` | bit-tree working columns (step 1 below); this file | merged 19 Nov 2026 (merge commit b706d25) |
| #2 | `cubical-2d` | planar engine for 2D images (union-finds, Alexander duality, unimodular certificate); the README cubical tables, re-measured in one sitting on top of #4 | merged 19 Nov 2026 (merge commit af0e491) |
| #3 | `simplicial-sets` | simplicial sets, chain complexes over Z, integral outputs, magnitude homology, Laurent coefficients | merged 19 Nov 2026 (merge commit 6237171) |
| #5 | `certify-warmstart` | classification profiled on Linial–Meshulam; big candidates in parallel (step 2) | merged 19 Nov 2026 (merge commit 00a8ec0) |
| #6 | `certify-proved-default` | the default certificate is a proof everywhere; fast mode only on request, labelled (step 3) | merged 19 Nov 2026 (merge commit d089aec) |
| #7 | `certificate-per-prime-theorem` | docs only: Theorem 4 (bad primes = primes of torsion in relative homology) in `CERTIFICATE.md` §7; step 4 plan | open; Vidit asked for this as a separate docs PR before step 4 |
| #8 | `step4-howell` | Howell sweeps: batched classification of small candidates (default on explicit complexes), `--torsion` trajectories; based on #7 | open |

The README's cubical tables (2D, 3D and 4D) are the single source of truth for cubical timings. They were re-measured on the code of #4 and #2 together, all programs in one sitting. Any later change that affects cubical speed re-measures the affected rows, or the whole table, and dates it.

## The plan, in order

The first items came from an external review of `main` at `5716401`; Vidit asked for an assessment, and the steps below are the agreed version. Prototype figures from the review are re-measured here before they are quoted. Steps 1 onward are ordered by what serves the research uses (torsion in random complexes, Milnor fibres, Novikov homology), not by benchmark value.

0. **Land #4, #2, #3 in that order** and re-measure every table once (above).

1. **Bit-tree working columns: done (#4).**
   * Long working columns in the cubical gradient and classic engines were binary heaps that keep duplicate rows until they surface, then sorted the rest to store the column. They are now PHAT's bit tree with a coefficient array (`include/perseus/bitcol.hpp`), serving F_2, F_p, multi-modular Q and Q64.
   * It is ring-generic, unlike the review's F_2-only bitset, because the default run is Q certified, which cost twice F_2 on white noise 48⁴.
   * With identical output: white noise 48⁴ F_2 went from 66.1 to 52.0 s (degree 1 from 26.3 to 7.9 s), and Q certified from 134.5 to 76.6 s. White noise 40⁴, degree 1, went from 6.5 to 2.4 s. The classic engine on white noise 32⁴ went from 20.3 to 12.4 s.
   * The explicit and Rips engines have no such heap (checked with grep for `push_heap`).
   * Afterwards, on white noise 48⁴ (F_2, 37.5 s), the largest phase is building the Morse complex: 14.8 s, 101 million path cells. The gradient takes 6.8 s and the reduction 12.8 s. This is the next lever for noisy 4D inputs; it is not scheduled.

2. **Make the classification of candidate primes cheap: profiled (19 Nov 2026); big candidates now run in parallel.**
   * **Instance.** `bench/certify_lm.py` regenerates it deterministically: all edges of K_90 at time 1, then the first 4,000 triangles of a random order. This is near the threshold where H_1 dies and torsion appears. Its F_2 time matches the README's n = 90 row, but it is not the same random complex: it has 496 candidates, not 463. There is also a smaller n = 70 instance.
   * **Profile before the change** (2 cores, total 147 s):

     | phase | time |
     |---|---|
     | Q run (8 lanes; 274 non-unit pivots, up to 170 bits) | 1.9 s |
     | factoring the pivot values into 386 small and 110 big candidates | 8.5 s |
     | 110 big candidates (≥ 2^63, runs modulo multi-word N), **one after another** | 80 s |
     | 386 small candidates, 2 threads (115 s CPU, 0.30 s per rerun) | 57 s |

     Six candidates are bad: 2, 5, 479, 22303, 1157339 and a 117-bit prime. Each of them changes exactly one H_1 interval and one H_2 interval near the end, where transient torsion in H_1 lives.
   * **Done.** The big candidates now go through a thread pool sharing the work queue; a split puts its two parts back on the queue. Total 147 → 105 s on 2 cores (big candidates 80 → 39 s). This scales with cores, and the output is identical for 1 and 2 threads.
   * **Warm starts, measured and not adopted.**
     * `PERSEUS_CLASSIFY_DUMP=file` (with `PERSEUS_DEBUG`) records, per small candidate, the columns of the pivots it divides and the cost of its rerun.
     * The explicit engine reduces in cohomology here, processing columns from high to low index, and its cost is back-loaded: 1% of the time has passed at column 803 of 4,005, 15% at 400, 56% at 200.
     * Each candidate's first divisible pivot comes at columns 89–803 (median 251). Restarting every rerun at that column would save 49% of the small-candidate CPU: 54 of 112 s.
     * But a restart needs the reduction's state modulo p at that column. In fast mode (int64 overflowed here) the Q run exists only modulo random 62-bit primes, so the state is not available. Exact checkpoints mean the proved run, which took 259 s on the README's instance.
     * The alternative, a run modulo a product of candidates that forks its state at each split, shares only the prefix among the 2–3 primes fitting in one word, so about 30% at best. It also needs resumable state in every engine.
     * So this is not worth building now.
   * **The bigger lever: now a theorem.** `docs/CERTIFICATE.md` §7 (Theorem 4) proves that p is bad in degree q iff some H_q(K_j, K_i; Z) or H_{q−1}(K_i; Z) has p-torsion (per-prime refinement of Obayashi–Yoshiwaki; checked independently against FLINT Smith forms on LM n = 28, 30). On LM only tors H_1(K_t; Z) matters, and the six bad primes on n = 90 all come from t ∈ [3564, 3568). The cheap necessary condition (p divides the product of pivot values over a value quadrant) excludes none of the 205 / 496 candidates, so classification through torsion needs actual torsion orders. That is step 4.
   * **Smaller remaining items:** factoring takes 8.5 s on one thread (parallelise across pivot values), and the small-candidate pool starts only after the big one finishes (one pool would avoid idle threads at the end).

3. **Make the default certificate exact all the way: done (19 Nov 2026).**
   * The default is now `proved`, and `auto` is a synonym. Random primes run only with `--mode fast`, and are labelled loudly: a WARNING on stderr, `"proof": false` and a `"warning"` entry in the JSON and in the Python certificate, and a Python `UserWarning`. Every certificate carries `"proof"`.
   * To make this affordable, exact arithmetic got:
     * inline limbs in `BigInt` (up to 256 bits);
     * a binary gcd on 128-bit words, with one-word shortcuts;
     * `QHyb`, rationals that are int64 pairs while they fit;
     * Henrici's formulas for sums and products.
   * The exact Q run on n = 90 went from 202 s to 24 s (fast mode: 2 s). Overall, proved costs 1.2× fast on n = 90 (126.7 s against 103.7 s) and 1.35× on n = 70 (6.9 s against 5.1 s), with the same bad primes.
   * Further: the exact run still spends most of its time in big-rational arithmetic. Fraction-free integer columns, which scale columns instead of keeping rationals, would remove most gcds; that is not done.
   * Tests: `t_exact` (BigInt identities, QHyb against plain big rationals) in `make test-chain`. The Python suite checks a torsion-heavy instance (Linial–Meshulam n = 60, int64 overflows): the default is proved with bad primes 2, 3, 922457715037; `auto` is the same; and fast agrees and is labelled in the module, the JSON and on stderr.

4. **Scale the integral machinery of #3, and classify candidates through torsion where it pays (CERTIFICATE.md §7). Howell sweeps and torsion trajectories: done (19 Nov 2026, branch `step4-howell`, on top of #7).**
   * **Classification by torsion orders: measured first (19 Nov 2026).** By Theorem 4 the bad primes in degree q are the primes of torsion in H_q(K_j, K_i; Z) and H_{q−1}(K_i; Z). For a value quadrant (rows born after s, columns born by t) of ∂_r, the torsion order of H_{r−1}(K_t, K_s) divides Δ(s, t) = ∏ π over the non-unit pairs in it. `PERSEUS_DUMP` now writes each non-unit pivot's birth and death; the counts come from `bench/quadrants.py` on the step-4 branch.

     | instance | values | quadrants with Δ ≠ ±1 (share) | non-unit pairs | candidates | bad |
     |---|---|---|---|---|---|
     | LM n = 70, tied (edges all at 1) | 2,401 | 727 (0.03%), one row threshold | 133 | 205 | 7 |
     | LM n = 90, tied | | one row threshold | 274 | 496 | 6 |
     | LM n = 40, distinct interleaved values | 1,480 | 1.2·10^5 (11%) | 30 | 6 | 6 |
     | LM n = 70, distinct interleaved values | 4,815 | 4.9·10^5 (4.3%) | 210 | 347 | 347 |
     | random metric, 300 points, H≤1 | 43,483 | 2.5·10^7 (2.7%) | 61 | 15 | 13 |
     | RP² Rips, 300 points, H≤2 | 18,452 | 8.6·10^7 (6–44% by degree) | 2 | 1 | 1 |
     | cubical white noise 96³ | 884,736 | 7.9·10^10 (20%) | 20 | 1 | 1 |
     | cubical white noise 24⁴ | 331,776 | 2.4·10^10 (14–29%) | 26 | 1 | 1 |

     The non-tied LM n = 40 result was checked independently by plain reductions over Q and F_p (`tests/fp_reference.py`): 2–13 bad, 17, 19, 23 harmless.
   * **Conclusion.** Torsion classification pays only where most candidates are harmless, and that is the tied case. On non-tied inputs almost every candidate is bad, and each bad prime needs its F_p barcode anyway (the certificate lists the changed intervals), so there is at most 2 of 15 reruns to save. Even so, the quadrant family there is two-parameter and huge. Do not build general two-parameter torsion classification.
   * **Vidit's decision (19 Nov 2026).** No general two-parameter classification. Build the sweep, starting with the tied case, validated against reruns and independent Smith forms. Expose the torsion trajectory itself (tors H_{d−1}(K_t; Z) with exponents at every t; for random-complex studies, Kahle–Lutz–Newman–Parsons). Handle d-dimensional Linial–Meshulam processes, and any filtration whose q-cells take few values.
   * **What was built (`howell.hpp`, `torsion.hpp`; CERTIFICATE.md §3 Lemma 5 and §7).**
     * *A Howell sweep* of a boundary matrix over Z/m, m any modulus (one word below 2^63, or BigInt): columns in filtration order, one basis vector per leading row with pivot ideal (d_r), extended gcds, annihilator multiples. For squarefree m, the rows with p ∤ d_r are the F_p lows (Lemma 5), so one sweep gives the F_p pairs of every prime in m. For m = ∏ p^k, Σ_r (k − min(k, v_p(d_r))) = log_p |span mod p^k|.
     * *Why not a plain run modulo a composite:* every candidate divides some pivot value, so such a run splits at every candidate (Lemma 3 style); the Howell form carries on through non-unit pivots instead.
     * *A correction to the plan above.* "One sweep over t replaces all 496 reruns" was wrong: it assumed one modulus holding every candidate, which is the modulus-Δ route that Vidit's second caution rules out. With 62-bit batches the cost is one sweep per batch and degree: 124 sweeps for 386 small candidates at n = 90.
     * *The tied-case restriction turned out unnecessary.* The sweep reads off F_p pairs directly, not relative torsion of quadrants, so it works for every filtration of an explicit complex. Absolute torsion needs only the column prefixes (s = −1), for any filtration too.
     * *Classification (`--classify-by auto|sweep|rerun`, default sweep on explicit complexes).* Small candidates are packed per degree into 62-bit moduli. A reference prime ℓ dividing no pivot value gives the Q pairs; its barcode must equal the Q barcode, or the code falls back to reruns. Big candidates (≥ 2^63) are still one run each modulo N.
     * *Torsion trajectory (`--torsion`, `torsion=True`).* For each bad prime p (Theorem 4: only bad primes have torsion), sweeps modulo p, then p^2, p^3, … until no exponent reaches k, give the p-primary exponents of tors H_q(K_t; Z) at every t. Output: per degree, the steps where the group changes, with invariant factors and primary parts. Composites that could not be factored are followed as such modulo N^k and flagged as groups. A split marks the trajectory incomplete.
   * **Validation.**
     * `tests/t_howell.cpp`: Lemma 5 against plain F_p reductions, and the length formula against brute-force spans; 6,000 random cases in both arithmetics; part of `make test-chain`.
     * Sweeps and reruns give identical bad primes and intervals on LM n = 40, 60, 70, 90 (tied and distinct values) and the Python fixtures.
     * The trajectory was checked against FLINT Smith forms at every t (`tests/torsion_flint.py`, needs python-flint, installed privately in `/tmp/pyflint`):
       * LM processes in dimension 2 (n = 18, 24) and 3 (n = 12, 13); only Z/2 at these sizes;
       * 20 random filtered chain complexes over Z, tied and distinct, invariant factors up to 1.7·10^7, primes up to 701, higher exponents.
     * On the benchmark processes, a full Smith form is too slow (over 20 min for one at n = 70). `tests/torsion_flint_large.py` instead eliminates unit entries exactly and runs FLINT on the residual (about 80 × 33). It agreed at 123 values at n = 70 and 71 at n = 90, all t in the torsion windows plus samples.
     * Python tests: a hand example, random chain complexes against a pure-Python Smith form (itself checked against FLINT on 3,000 matrices), the 98-bit-prime and composite-group fixtures, and command line = module.
   * **Results (`bench/certify_lm.py`, 2 cores, under the bench lock).**

     | instance | default (sweeps) | reruns | small candidates: sweeps / reruns | `--torsion` adds |
     |---|---|---|---|---|
     | LM n = 70, tied | 4.3 s | 7.0 s | 2.1 s (61) / 4.8 s (202) | 0.3 s |
     | LM n = 90, tied | 93.8 s | 125.5 s | 24.4 s (124) / 56.4 s (386) | 19.9 s |
     | LM n = 70, distinct values | 18.7 s | 31.4 s | 7.3 s (114) / 20.1 s (326) | 21.7 s |

     Trajectories found:
     * n = 70: H_1 passes through Z/2 ⊕ Z/11618403313260425604101352.
     * n = 90: a cyclic group of order 1.09·10^50, including the 117-bit prime, for three steps.
     * 3-dimensional LM, n = 25, 2,700 tetrahedra: H_2 passes through Z/8705060926506448518465431062162, with an 81-bit composite followed as a group.
   * **Next, in order of payoff.**
     1. The 110 big candidates at n = 90 (39 s) are now the largest classification cost. Packing them into multi-word moduli does not pay, since schoolbook arithmetic is quadratic in the words. Options: Montgomery arithmetic for `ModBig` (it uses BigInt division now), and running them in the Howell sweep so that a split no longer restarts a run.
     2. `ModBig` sweeps are slow (most of the trajectory's 20 s at n = 90 is the 117-bit prime). Montgomery multi-word arithmetic, as in `run_mod_big`, would fix this.
     3. Rips and cubical inputs still rerun every candidate. They have few candidates (1–15 in the measurements), so this is low priority.
   * The integral outputs use a dense Smith form, capped at about 40M entries, and run only on explicit complexes. Replace this with sparse elimination down to a small residual, followed by a Smith form; the fill-in-free first phase of #3 is the starting point.
   * Then extend the integral outputs to cubical grids, which have torsion from dimension 3.

5. **Laurent (Novikov) coefficients, driven by applications.** Use SnapPy knot-complement triangulations and a Salvetti-complex model of a Milnor fibre as tests, and add the missing pieces in the order those tests need them:
   * irreducibility proofs beyond degree 2 (Berlekamp–Zassenhaus factoring over Z);
   * images over Q[t, t⁻¹];
   * the Z[t, t⁻¹]-module structure (residue fields F_p[t]/(h)).

6. **Parallel column reduction** (Morozov–Nigmetov lock-free reduction, as in giotto-ph) for Rips and cubical. It is the most invasive change, and it multiplies whatever the earlier steps gain, so it comes after them.

7. **Small items, opportunistic or dropped:**
   * *H0-only Rips by a minimum spanning forest.* Prim's algorithm takes O(n²) time and O(n) memory, and the certificate is unimodular (graph incidence). It only helps `--maxdim 0`, and must match the enclosing-radius threshold and zero-length bars. Re-measure the review's 2.4–10×.
   * *Merging the cubical setup passes.* Births and gradient metadata could be written in one pass. The review's gain was in 2D, which the planar engine made moot. Drop it unless smoothed noise 192³ or 48⁴ shows a real gain; the review saw none on white noise 24⁴.
   * *Two threads in the planar engine:* dropped (see below).

8. **Release 6.1** when Vidit decides: PyPI, a new Zenodo version, the web page.

## Tried and measured, not adopted (do not repeat without a new idea)

* **Planar engine (PR #2):**
  * compact node records indexed by filtration rank, and remapping through rank tables: slower;
  * a full-grid pointer-compression pass: slower on every input. Compression from the candidates' neighbours only is kept, and switched on automatically for large sparse inputs;
  * prefetch distance tuning and second-stage prefetch: no measurable change at 4096²;
  * pointing at the middle of the neighbour arc: no measurable change.
  * Running H0 and H1 on two threads conflicts with the H1-first skip; it would save about 0.35 s at 4096². Dropped: step 6 supersedes it.
* **Degree-1 reduction of white noise 48⁴ (before the bit tree):**
  * storing the reductions V instead of the reduced columns, as Ripser does: 5–9× more merged entries, 10× slower;
  * reducing degree 1 by homology: 83 s against 29 s.
* **Bit-tree threshold:** 128–512 entries measured alike; 2048 and 8192 were slower; 512 is the default.

## Test commands

* `make perseus`: the command-line binary.
* `make test-cubical`:
  * `t_cubical_morse` (gradient against classic over F_2, F_3, F_5 and Q, with certificates), also run with `PERSEUS_LONG_COLUMN=0`, which sends every column through the bit tree;
  * `t_bitcol` (the bit tree against a `std::map`, over F_2, F_3, F_7, multi-modular Q and Q64);
  * on PR #2, also `t_cubical_planar`.
* `python3 tests/run_tests.py B E F G`: differential tests. Sections A, C and D need Ripser, p2ref and Perseus 5, which are not on the test machine.
* `PYTHONPATH=<private target> PERSEUS=./perseus python3 tests/test_python.py`: Python tests.
* `PYTHONPATH=<private target> python3 tests/fuzz_cubical.py N [seed]`: fuzzing against CubicalRipser and between engines.
