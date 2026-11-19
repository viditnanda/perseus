# Performance roadmap: status, decisions and protocol

This file is for anyone picking up Perseus 6 performance work, whether a person or a Claude instance. It records what is in flight, what was agreed with Vidit Nanda (the maintainer), how to measure, and what has already been tried. Update it in the same pull request as the work it describes.

Last updated: 19 Nov 2026.

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

## Open pull requests

| PR | branch | content | state |
|---|---|---|---|
| #2 | `cubical-2d` | planar engine for 2D images (union-finds, Alexander duality, unimodular certificate); new 2D table; re-measured GUDHI and CubicalRipser columns of the 3D/4D table | awaiting review |
| #3 | `simplicial-sets` | simplicial sets, chain complexes over Z, integral outputs, magnitude homology, Laurent coefficients | awaiting review |
| (this) | `cubical-bittree` | bit-tree working columns, item 1 below | awaiting review |

`cubical-bittree` branches from `main` at `5716401` and does not depend on #2 or #3. It touches `cubical_morse.hpp` and `cubical_engine.hpp`, which #2 also edits. Expect a small merge conflict in `cubical_engine.hpp` (#2 adds the planar dispatch to the backend; this branch changes the classic engine's reduction loop) and in the README's 3D/4D section.

## Agreed plan, in order

This plan came from an external review of `main` at `5716401`. Vidit asked for an assessment, and the items below are the agreed version, with the review's prototype figures re-measured here.

1. **Bit-tree working columns: done on `cubical-bittree`.**
   * Long working columns in the cubical gradient and classic engines were binary heaps that keep duplicate rows until they surface, and then sorted the rest to store the column. They are now PHAT's bit tree with a coefficient array (`include/perseus/bitcol.hpp`). The same structure serves F_2, F_p, multi-modular Q and Q64.
   * The review proposed an F_2-only bitset. We made it ring-generic because the default run is Q certified, which cost twice F_2 on white noise 48⁴.
   * Results, with identical output:

     | case | before | after |
     |---|---|---|
     | white noise 48⁴, F_2 (total) | 66.1 s | 52.0 s |
     | white noise 48⁴, F_2 (degree 1) | 26.3 s | 7.9 s |
     | white noise 48⁴, Q certified | 134.5 s | 76.6 s |
     | white noise 40⁴, degree 1 | 6.5 s | 2.4 s |
     | classic engine, 32⁴ | 20.3 s | 12.4 s |
   * Not yet done: the explicit and Rips engines do not have this heap pattern (checked with grep for `push_heap`), so nothing was changed there.

2. **Certificate classification: next.** On Linial–Meshulam n = 90 (README, explicit complexes), the Q reduction takes 1.3 s of about 113 s; the rest classifies 463 candidate primes by reruns (`certify.hpp`, classification loop).
   * Large candidates are already batched: a run modulo their product splits it on a zero divisor.
   * The review proposed batching small primes by a composite modulus. The catch is that every candidate divides some pivot value by construction, so every batch splits eventually, and the only saving is the shared prefix before the first split.
   * The better lever is a warm start. For each prime, restart the reduction at the first column whose pivot it divides, using checkpoints of the Q reduction. The Q run already records the columns with non-unit pivots (`out.nonunit`), so factoring those pivot values gives each candidate its first such column.
   * **Measure first:** the distribution of those first columns, and how many of the 463 candidates are actually bad. Random 2-complexes have large torsion, and each bad prime needs its full barcode modulo p.

3. **H0-only Rips by a minimum spanning forest.** The H0 barcode of a Rips filtration is the minimum spanning forest (Prim, O(n²) time, O(n) memory, distances computed on demand).
   * It applies only to `--maxdim 0`; the default is 1.
   * The certificate is the unimodular one, since this is a graph incidence matrix (see `docs/CERTIFICATE.md` §5 on PR #2).
   * It must match the existing default threshold (the enclosing radius) and zero-length bars.
   * The review measured 2.4–10× on 3,000 and 6,000 points; re-measure here.

4. **Merging the cubical setup passes: low priority.** Births and gradient metadata could be written in one pass. The review's gain (smoothed noise 4096², 2.93 to 2.29 s) is in 2D, where the planar engine (PR #2) now takes 0.24 s, so it matters only for smooth 3D and 4D inputs. Test on smoothed noise 192³ and 48⁴ before investing; the review saw no change on white noise 24⁴.

## Tried and measured, not adopted (do not repeat without a new idea)

* **Planar engine (PR #2):**
  * compact node records indexed by filtration rank, and remapping through rank tables: slower;
  * a full-grid pointer-compression pass: slower on every input. Compression from the candidates' neighbours only is kept, and switched on automatically for large sparse inputs;
  * prefetch distance tuning and second-stage prefetch: no measurable change at 4096²;
  * pointing at the middle of the neighbour arc: no measurable change.
  * Running H0 and H1 on two threads conflicts with the H1-first skip; it would save about 0.35 s at 4096². Not done.
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
