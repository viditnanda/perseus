# Changes

## 6.1.0 (unreleased)

Measured against 6.0.0 on the same 2-core machine, with identical barcodes and bad primes in every case
(`tests/` and the head-to-head in the 6.1 pull request).

### Certificates

* **The default certificate is a proof everywhere.** In 6.0, `--mode auto` used exact int64 arithmetic and
  fell back to random primes where it overflowed, as on torsion-heavy inputs. Now `auto` is the same as
  `proved`: rationals stay in int64 pairs while they fit and become big integers only where they must, so the
  barcode and the bad primes are always proved. Random primes run only with `--mode fast`, which says on
  stderr, in the JSON (`"proof": false`, `"warning"`) and in Python (a `UserWarning`) that it is not a proof.
  Every certificate now carries `"proof"`.
* Exact arithmetic is much faster (inline 256-bit integers, binary gcds on 128-bit words, Henrici's
  cancellation), so the proof costs little: a Linial–Meshulam complex on 60 vertices is proved in 1.0 s,
  where 6.0 took 1.4 s for a probabilistic answer.
* Classifying candidate primes: candidates of at least 2^63 run in parallel; on explicit complexes the small
  candidates can be settled in batches by Howell sweeps over Z/m, chosen automatically when a timed sample
  predicts a saving (`--classify-by auto|sweep|rerun`).
* `docs/CERTIFICATE.md` §7 proves that p is bad in degree q exactly when some H_q(K_j, K_i; Z) or
  H_{q−1}(K_i; Z) has p-torsion.

### New: integral outputs and new kinds of input

* `--torsion`: the torsion of H_q(K_t; Z) at every t where it changes, for any explicit filtration.
* `cells` is now a filtered chain complex over Z, reduced first by algebraic Morse theory; `--integral` gives
  H_*(K; Z) and the F_p barcode of every bad prime, `--image s,t` the image of H_*(K_s; Z) → H_*(K_t; Z).
* `sset`: finite filtered simplicial sets. In Python also nerves of finite categories (`perseus.nerve`) and
  bar constructions (`perseus.bar_construction`).
* `magnitude`: magnitude homology of graphs over Z, torsion included.
* `laurent`: chain complexes over Z[t, t⁻¹] (infinite cyclic covers, Alexander modules, Novikov homology).
* Finite groups: `pcgroup` (solvable groups: H_*(G; Z) and mod-p cohomology rings through the Anick
  resolution), `sylow` (any finite group through its Sylow subgroups), `pcseries` (maps on H_*(−; F_p) along a
  chain of homomorphisms, and their persistence). GAP exporters are in `tools/`.

### Speed

* Cubical, 2D: a new planar engine (union-finds and Alexander duality; no arithmetic, since the boundary
  matrices are totally unimodular). On 1024² and 2048² images it is 6.5–26× faster than 6.0, in less than half the memory, and
  faster than GUDHI's image algorithm on every image in the README.
* Cubical, 3D and 4D: a filtered discrete-gradient (Morse) engine. 2.2–3.9× faster than 6.0, with 25–63% less
  memory, and faster than CubicalRipser and GUDHI on every grid in the README.
* Bit-tree working columns for long reductions (after PHAT), for every coefficient ring.
* `cubes`: a packed builder for sparse cubical complexes.
* Rips: disconnected filtrations are split into components for single-field computations.
* Writing output files is about 8× faster (the shortest round-trip decimals are found with `<charconv>`);
  the files are byte-for-byte the same as 6.0's. On 1024² white noise, 6.5 s → 0.75 s, almost all of which
  was writing.

### Documentation

* `examples/`: nine getting-started examples, each as a shell script and a Python script, with their inputs;
  `make examples` runs them all, and CI runs them on every push.
* `docs/TECHNIQUES.md` catalogues every technique with its origin; `docs/ROADMAP.md` records the plan and
  what was tried.

### Compatibility

* `--mode auto` now proves (see above); scripts that relied on the fallback to random primes for speed should
  pass `--mode fast`.
* The `examples/` input files moved to `examples/data/`.
* Command-line options, output files and existing Python functions are otherwise unchanged.

## 6.0.0 (17 November 2026)

First release of the rewrite. [doi:10.5281/zenodo.24305983](https://doi.org/10.5281/zenodo.24305983)
