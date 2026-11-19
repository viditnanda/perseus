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
| #4 | `cubical-bittree` | bit-tree working columns (step 1 below); this file | verified by Vidit, to be merged first |
| #2 | `cubical-2d` | planar engine for 2D images (union-finds, Alexander duality, unimodular certificate); the README cubical tables, re-measured in one sitting on top of #4 | rebased onto #4; force-pushed once #4 is on `main` |
| #3 | `simplicial-sets` | simplicial sets, chain complexes over Z, integral outputs, magnitude homology, Laurent coefficients | to be rebased onto `main` after #2, then the same independent check |

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

2. **Make the classification of candidate primes cheap.** The bad-prime certificate is the core of the enterprise. On Linial–Meshulam n = 90, the Q reduction takes 1.3 s of about 113 s; the rest reruns the reduction for each of 463 candidate primes (`certify.hpp`, classification loop).
   * **Profile first:** for each candidate, the first column whose pivot it divides, and how many of the 463 are actually bad.
   * **Then warm-start each prime** from that column, using checkpoints of the Q reduction. Up to that column, the reduction modulo p is identical to the one over Q. The Q run already records its non-unit pivots (`out.nonunit`).
   * **Composite-modulus batching is the weaker idea.** Every candidate divides some pivot by construction, so every batch splits eventually, and only the prefix before the first split is shared. Large candidates are already batched this way.
   * Random 2-complexes have large torsion. Each bad prime still needs its full barcode modulo p, and only the warm start helps there.

3. **Make the default certificate exact all the way.** Vidit decided that certified Q should end in exact arithmetic. Today, when int64 arithmetic overflows, `--mode auto` falls back to `fast`, which is probabilistic. Make the fallback `proved` (big-integer arithmetic), or at minimum make the fallback impossible to miss in the output and in `cert.json`.

4. **Scale the integral machinery of #3.**
   * The integral outputs use a dense Smith form, capped at about 40M entries, and run only on explicit complexes. Replace this with sparse elimination down to a small residual, followed by a Smith form; the fill-in-free first phase of #3 is the starting point.
   * Then extend the integral outputs to cubical grids, which have torsion from dimension 3.
   * Then test whether integral persistence can read off bad primes directly, without reruns. That would make step 2 partly redundant on explicit complexes.

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
