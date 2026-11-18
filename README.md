# Perseus 6

Persistent homology **over Q, with a certificate of every prime p for which the F_p barcode is different**, or over any prime field.

One run gives you:

* the rational barcode,
* the complete list of *bad primes* — the p for which F_p coefficients give a different barcode — and, for each, exactly which intervals change,
* the exact rational values of all non-unit pivots, which explain where the torsion lives.

By default the Q reduction runs in exact int64 arithmetic. Then the certificate is a proof, at about the cost of an F_p run. When the numbers outgrow int64, as on torsion-heavy inputs, Perseus switches to a run modulo random 62-bit primes and reports a rigorous bound on the probability of error. `--proved` insists on exact arithmetic throughout. See [docs/CERTIFICATE.md](docs/CERTIFICATE.md).

Perseus 6 is a rewrite from scratch. It does not read the old Perseus 4/5 input formats directly; `tools/convert_legacy.py` turns old files into new ones (see [Old Perseus files](#old-perseus-files)).

## Install

Executables for Linux, macOS and Windows are attached to each [GitHub release](https://github.com/viditnanda/perseus/releases); the Python module is `pip install perseus-tda` (Linux and macOS).

## Build

```
make                 # builds ./perseus (C++17 with 128-bit integers: g++ or clang; MinGW-w64 on Windows)
make test            # differential tests (needs Ripser / Perseus 5 / p2ref; skips what is missing)
cd python && pip install --no-build-isolation .   # Python module `perseus` (needs pybind11, numpy)
```

## Quick start

```
perseus points cloud.txt --maxdim 2                  # Rips over Q + bad-prime certificate (stdout: Ripser style)
perseus points cloud.txt --maxdim 2 --coeff 2        # just F_2
perseus bin32 dist.bin --maxdim 3 --threshold 1.8 out   # writes out_k.txt, out_betti.txt, out_cert.json
perseus simplices complex.txt --fp 2,3               # also print the F_2 and F_3 barcodes
perseus npy image.npy --periodic 0                  # cubical complex of an image, periodic in axis 0
```

On stderr you get a summary line:

```
Q barcode: proved (exact arithmetic: int64)
bad primes: 2  (2 non-unit pivots)
```

For example, 1024 random points of SO(3) ≅ RP³ (`bench/o3_1024`, H≤3, threshold 1.8) give two pivots, −2 and 2. Each sits on one H1 interval:

* over Q: H1 has the finite bars [0.657, 1.288) and [0.690, 1.307);
* over F_2: those two H1 bars are infinite, and H2 gains [1.288, ∞) and [1.307, ∞).

## What is certified

Run the reduction over Q and look at the pivots. If a prime p divides no numerator and no denominator of a pivot value, then the same run reduced mod p is a valid F_p reduction with the same pivots. So the F_p barcode equals the Q barcode, and the bad primes are among the prime factors of the pivot values. Pivot values are invariants of the filtered matrix (ratios of minors), whichever reduction computes them. Perseus 6:

1. **Reduces over Q.** By default (`--mode auto`) it tries exact int64 arithmetic first. If that suffices, the barcode and pivots are exact. Otherwise it reduces modulo K random primes from [2^61, 2^62) and reports a bound on the probability of error that holds for every input (default below 2^-64; `--target`). The bound is over a fresh 256-bit seed; `--seed` reproduces a run but does not certify it. `--mode proved` uses exact arithmetic throughout, with bigints where int64 overflows. `--mode fast` always uses random primes.
2. **Reconstructs** non-unit pivot values in fast mode, by Wang reconstruction checked against two independent primes.
3. **Factors** the pivot values: trial division, division by primes already found, Pollard–Brent rho, and gcds between leftovers. Large primes and unsplit composites N are not factored further.
4. **Settles every candidate.** Primes below 2^63 are settled by reruns over F_p, in parallel. Each remaining N is settled by one run modulo N. Either the run splits N, or every prime factor of N gives the same barcode. Then the certificate states the answer for all of them at once, e.g. "every prime factor of N is bad", without knowing the factors.

The output states what is guaranteed about the barcode (`proved`, or `probabilistic` with the bound), and whether the description of the bad primes is complete. Field dependence of persistence is governed by torsion in the relative homology H_q(X_n, X_m; Z) (Obayashi–Yoshiwaki). The pivots above are where that torsion is detected.

## Inputs

| format | contents |
|---|---|
| `points` | one point per line, whitespace or comma separated (Euclidean Rips) |
| `distance` | full distance matrix |
| `lower-distance` | strictly lower triangle of a distance matrix, row by row (Ripser's format); `inf` = no edge |
| `bin32`, `bin64` | binary lower triangle, float32 / float64 (Ripser's `binary` format) |
| `flag` | clique filtration of a weighted graph: lines `u v birth` (edges) and `v birth` (vertex births, default 0); labels are arbitrary tokens |
| `cubical` | dense grid of top-cell values: a line of sizes `n_1 ... n_d`, then the values in C order (last axis fastest); `nan`/`inf` = absent; `--periodic 0,2` makes axes periodic |
| `npy` | the same from a numpy `.npy` file (any real or integer dtype, C or Fortran order) |
| `cubes` | sparse cubical complex: one top cube per line, its integer anchor coordinates then its birth |
| `simplices` | one simplex per line: vertex labels, then its birth; faces are added with the minimum birth of their cofaces |
| `cells` | general cell complex with integer coefficients: one cell per line `dim birth [face coef]...`, faces are 0-based line numbers of earlier cells |

`#` starts a comment. Births are arbitrary real numbers. Coefficients default to Q (with the bad-prime certificate); `--coeff p` gives F_p.

## Outputs

* `prefix_k.txt`: one line `birth death` per interval, with death `-1` if infinite.
* `prefix_betti.txt`: Betti numbers at every time they change.
* `prefix_cert.json` (Q mode). Fields:
  * `barcode`
  * `mode` (`proved` or `fast`) and `guarantee` (`barcode_certified`, and in fast mode `log2_error_bound`, `seed`, `seed_supplied`, `moduli`)
  * `bad_primes_status` (`complete`, `candidates only`, `incomplete`) and `complete`
  * `nonunit_pivots` (degree, birth, death, exact value)
  * `candidates`
  * `bad_primes`: for each p, the intervals present only over F_p and only over Q. An entry with `"every_prime_factor_of": true` applies to every prime factor of the number given.
  * `unresolved_factors`
  * timings

`--print ripser|perseus|json|none` controls stdout.

## Python

```python
import numpy as np, perseus
r = perseus.rips(points=X, maxdim=2)            # or distances=D (square or condensed)
r["barcode"][1]                                 # (m, 2) array, death = inf if essential
r["certificate"]["mode"]                        # "proved", or "fast" with r["certificate"]["log2_error_bound"]
perseus.bad_primes(r)                           # e.g. [2]
r["certificate"]["bad_primes"][2]["only_over_Fp"]
perseus.cubical(image, periodic=[True, False])  # NaN/inf = absent top cell
perseus.simplices([((0, 1, 2), 1.0), ...])
perseus.flag([(0, 1, 0.5), (1, 2, 0.7)], vertex_births=[0, 0, 0.2])   # clique filtration of a graph
perseus.cells([(0, 0.0, []), (1, 1.0, []), (2, 2.0, [(1, 12)])])   # Moore space M(Z/12, 1): bad primes [2, 3]
perseus.rips(points=X, coeff=3)                 # plain F_3
perseus.rips(points=X, mode="proved")           # exact arithmetic throughout (also "fast", seed=..., target=...)
```

Every function takes `coeff`, `mode`, `seed`, `target`, `classify`, `threads` and `fp`, as on the command line; `rips` and `flag` also take `engine` (`auto`, `implicit`, `sparse`) and `collapse`, and `cubical` takes `engine` (`auto`, `planar`, `morse`, `classic`, `explicit`). For the planar engine (the 2D default) the array is read in place, without a copy. A key of `certificate["bad_primes"]` whose entry has `every_prime_factor_of` set is a number all of whose prime factors are bad; `probable_prime` marks a prime above 2^64 known only by Miller–Rabin. `tests/test_python.py` checks the module against the command line.

## Engines

**Explicit complexes** (simplicial, cubical, general cells). Three reductions, all certified by the same layer, with `--direction auto` choosing between the first two:

* **Backward cohomology with clearing.** This is the default.
* **Forward cohomology.** Cells are inserted in filtration order. Each cell carries its coordinates in the basis of currently alive classes, with an inverted index. It is chosen when many top cells sit over a degree whose cells arrive in a few large batches, as in Linial–Meshulam complexes.
* **Homology with clearing from the top.** Available on request.

H0 uses union-find whenever every 1-cell has boundary ±(v − w).

**Cubical grids** (`cubical`, `npy`): three implicit engines; none builds the complex.

* **Planar** (the default for 2D grids without periodic axes; `--engine planar`). H0 is a union-find on pixels: two closed squares meet exactly when their pixels are 8-adjacent. H1 is, by Alexander duality in the sphere, a union-find on the complement swept downwards, with 4-adjacency and one node for the outside; absent pixels give essential classes. One byte per pixel, its earlier neighbours, and two 256-entry tables tell which pixels can change anything. Only those are sorted. Degree 1 runs first, and the Euler characteristic of each pixel's neighbourhood lets the degree-0 pass skip most of them. Both union-finds reduce graph incidence matrices, which are totally unimodular. So the barcode is the same over Q and every F_p, and the certificate needs no arithmetic (mode `unimodular`). This follows the dual-then-primal structure of GUDHI's image algorithm (Glisse). See `docs/TECHNIQUES.md` §4.7 and `docs/CERTIFICATE.md` §5.
* **Discrete gradient** (the default in dimensions 3 and 4 without periodic axes; `--engine morse`, also in 2D). A cell enters with its first top coface τ. The cells entering with τ are matched among themselves by a per-pattern discrete gradient, which depends only on which of τ's 3^d − 1 neighbours came earlier. This is the filtered Morse reduction of Mischaikow–Nanda. Every pattern in 2D and 3D, and every pattern we have met in 4D, is matched perfectly: it leaves exactly as many critical cells as bars start or end at τ. H0 and degree d−1 are union-finds on critical cells along gradient paths. The degrees in between are cohomology with clearing on the integer Morse complex, which is filtered chain-homotopy equivalent to the cubical complex. So barcodes, Q certificates and bad primes are the same as with the classic engine. See `docs/TECHNIQUES.md` §4.6.
* **Classic** (periodic axes, dimension 1 or above 4; `--engine classic`). It stores one birth and one rank per cell of the doubled grid and generates coboundaries on the fly. H0 is union-find. Degree d−1 is union-find on the dual graph: every (d−1)-cell has at most two cofaces, so δ^(d−1) is an incidence matrix and all its pivots are ±1. The degrees in between are cohomology with clearing.

`--engine explicit` builds the complex instead.

**Sparse Rips / flag** (`--engine sparse`; automatic when the flag complex is small next to n²: n ≥ 2048, average degree at most 64 at the threshold, and at most n²/50 simplices, counted before anything is built): the radius graph is built directly, with a grid for low-dimensional point clouds, and the flag complex is enumerated from adjacency lists. Memory is linear in the number of simplices. Example: 8,000 points in the unit cube with 13,056 edges at r = 0.047 take 0.04 s and 13 MB through H2. The dense engine needs 0.7 s and 750 MB for the same input.

**Vietoris–Rips** (dense): the complex is never built. For maxdim ≥ 2 it is first edge-collapsed (`--collapse auto|on|full|off`), then one of two engines runs:

* **Cone-level cohomology.**
  * The filtration is refined by *levels*, i.e. the rank of the longest edge.
  * A level is a *cone level* when its link has an apex adjacent to everything else in the link. This is the edge-domination condition of Boissonnat–Pritam.
  * In the apex-first order, the cells of a cone level pair up as apparent pairs. They are never enumerated, and their rows never become pivots of base cells.
  * Apparent pairs in the remaining levels are recognised implicitly.
  * The design follows Ripser (lazy coboundaries, clearing, apparent pairs).
* **Active H1 enumeration** in the homology direction, for `--maxdim 1`, in the style of RedZeD.

## Benchmarks

### What certification costs (measured 16 Nov 2026, 2-core Xeon, wall clock)

| case | F_2 | Q auto (default) | Q fast | Q proved | original 6.0 Q (unsound) |
|---|---|---|---|---|---|
| o3_1024 H≤3, r = 1.8 | 0.06 s | 0.07 s, proved | 0.09 s | 0.07 s | 0.07 s |
| o3_4096 H≤3, r = 1.4 | 0.85 s | 0.99 s, proved | 1.03 s | 0.99 s | 0.92 s |
| sphere H≤2 | 0.42 s | 0.45 s, proved | 0.44 s | 0.45 s | 0.47 s |
| rp2_600 H≤2 | 11.5 s | 19.2 s, proved | 27.5 s | 19.7 s | 19.1 s |
| random metric, 300 points, H≤1 | 0.15 s | 2.1 s, proved | 2.0 s | 2.1 s | 1.8 s |
| Linial–Meshulam n = 70, H≤2 | — | 4.5 s, fast (int64 overflows) | 4.4 s | 12.1 s | — |
| Linial–Meshulam n = 90, H≤2 | 0.16 s | ≈ 113 s, fast | 113 s | 259 s | does not finish |

On geometric data a proof costs nothing extra; time beyond F_2 is classification of candidate primes. Exact arithmetic gets expensive only when the integers involved grow, which happens on torsion-heavy inputs, and then the default switches to random primes. On Linial–Meshulam complexes the Q reduction itself is fast (1.3 s at n = 90, with 261 non-unit pivots of up to 170 bits); classifying the 463 candidate primes, of which exactly two are bad, takes the rest.

### Cubical, 2D images (F_2, measured 18 Nov 2026, 2-core Xeon, 7.8 GB)

`bench/cubical2d.py` regenerates the inputs: the 512² camera photograph from scikit-image (8-bit grey levels, so heavy ties), uniform noise, and uniform noise smoothed by a Gaussian of width 4, all as float64. Every program runs in a fresh process and is timed around its in-process Python call, after its imports, with the input array already in memory. Memory is the process's peak RSS right after the call, so it includes the interpreter, the modules and the input array. Each time is the best of three runs (two for white noise 8192²). GUDHI's `CubicalPersistence` is asked for degrees 0 and 1, which selects its dedicated image algorithm. The Perseus command line is timed separately, including process start and reading the `.npy` file, with output suppressed. In every case the positive-length intervals are identical to GUDHI's, and, up to 4096², to CubicalRipser's and the gradient engine's. At 8192² those two need more memory than this machine has.

| image | Perseus 6, planar engine | planar, Q certified | Perseus command line | GUDHI 3.13 | CubicalRipser | Perseus 6, gradient engine | intervals |
|---|---|---|---|---|---|---|---|
| photograph 512² | **9.5 ms / 39 MB** | 11.7 ms / 39 MB | 0.02 s / 14 MB | 13.0 ms / 147 MB | 66.1 ms / 101 MB | 38.6 ms / 47 MB | 36,849 |
| smoothed noise 1024² | **12.3 ms / 52 MB** | 13.9 ms / 52 MB | 0.03 s / 33 MB | 21.2 ms / 162 MB | 0.40 s / 265 MB | 0.17 s / 111 MB | 7,022 |
| white noise 1024² | **65.4 ms / 66 MB** | 65.8 ms / 66 MB | 0.08 s / 40 MB | 87.6 ms / 175 MB | 0.54 s / 302 MB | 0.29 s / 121 MB | 325,501 |
| smoothed noise 4096² | **0.24 s / 252 MB** | 0.24 s / 251 MB | 0.43 s / 260 MB | 0.53 s / 499 MB | 9.91 s / 3.0 GB | 3.42 s / 1.0 GB | 110,329 |
| white noise 4096² | **1.40 s / 450 MB** | 1.48 s / 450 MB | 1.71 s / 424 MB | 1.82 s / 682 MB | 10.13 s / 3.0 GB | 7.25 s / 1.4 GB | 5,216,867 |
| smoothed noise 8192² | **0.99 s / 892 MB** | 0.96 s / 893 MB | 1.80 s / 1.0 GB | 2.48 s / 1.5 GB | out of memory | out of memory | 441,508 |
| white noise 8192² | **7.22 s / 1.6 GB** | 7.25 s / 1.6 GB | 7.57 s / 1.6 GB | 9.33 s / 2.3 GB | out of memory | out of memory | 20,876,216 |

* The planar engine is 1.3–2.5× faster than GUDHI's image algorithm and uses 27–73% less memory: 1.3× on white noise, 1.7–2.5× on smoothed noise, 1.4× on the photograph. It is 7–41× faster than CubicalRipser.
* The certificate costs nothing: the boundary matrices it reduces are totally unimodular, so the Q barcode is proved equal to every F_p barcode without arithmetic (`docs/CERTIFICATE.md` §5).
* Where the time goes (white noise 4096², 1.25 s in C++): classification 0.13 s, sorting 6.4 million births and candidates 0.22 s, degree 1 0.59 s, degree 0 0.28 s. Only 38% of the pixels are sorted; on smoothed noise, under 1%.

### Cubical, 3D and 4D (F_2, measured 18 Nov 2026, 2-core Xeon, 7.8 GB)

`bench/cubical.py` regenerates the inputs (uniform noise, and noise smoothed by a Gaussian of width 2 in 3D and 1.5 in 4D; float64 values) and runs every program in a fresh process. Time is wall clock. Perseus is timed from the command line, including process start and reading the `.npy` file, and its memory is the peak RSS. CubicalRipser and GUDHI are timed around the in-process call, after their imports, and their memory is the peak RSS right after the call, including Python and the input array. Their columns were re-measured on 18 Nov 2026, after an earlier version of the script was found to time GUDHI's import (about a second) and to count the barcode comparison in their memory. In every case Perseus (both engines) and CubicalRipser give identical barcodes; on white noise 48⁴, where CubicalRipser did not finish, the two Perseus engines agree. All three compute the T-construction.

| case | Perseus 6, gradient engine | Perseus 6, classic engine | CubicalRipser (`tcripser`) | GUDHI 3.13 (`CubicalPersistence`) | critical cells | gradient engine, Q certified |
|---|---|---|---|---|---|---|
| smoothed noise 96³ | **0.36 s / 96 MB** | 1.04 s / 143 MB | 0.88 s / 132 MB | 4.70 s / 469 MB | 0.22% | 0.46 s / 97 MB |
| white noise 96³ | **1.14 s / 134 MB** | 2.95 s / 180 MB | 2.50 s / 276 MB | 8.16 s / 553 MB | 8.7% | 1.34 s / 147 MB |
| smoothed noise 24⁴ | **0.36 s / 63 MB** | 0.89 s / 126 MB | 0.68 s / 121 MB | 4.27 s / 396 MB | 0.13% | 0.40 s / 64 MB |
| white noise 24⁴ | **1.63 s / 108 MB** | 3.66 s / 182 MB | 3.02 s / 274 MB | 7.58 s / 445 MB | 4.7% | 1.99 s / 123 MB |
| smoothed noise 192³ | **3.88 s / 764 MB** | 13.34 s / 1.1 GB | 11.51 s / 800 MB | 56.99 s / 2.7 GB | 0.23% | 3.45 s / 764 MB |
| white noise 192³ | **11.25 s / 1.0 GB** | 27.19 s / 1.4 GB | 25.30 s / 1.7 GB | 92.47 s / 3.4 GB | 9.0% | 17.60 s / 1.3 GB |
| smoothed noise 32⁴ | **1.63 s / 183 MB** | 4.35 s / 386 MB | 3.24 s / 276 MB | 15.67 s / 923 MB | 0.13% | 1.69 s / 183 MB |
| white noise 32⁴ | **6.19 s / 342 MB** | 19.94 s / 626 MB | 14.40 s / 903 MB | 37.99 s / 1.1 GB | 5.0% | 9.56 s / 427 MB |
| smoothed noise 48⁴ | **7.87 s / 923 MB** | 24.70 s / 1.7 GB | 21.33 s / 1.2 GB | 90.13 s / 4.0 GB | 0.14% | 7.67 s / 923 MB |
| white noise 48⁴ | **64.76 s / 1.7 GB** | 253.92 s / 3.8 GB | killed (out of memory above 6 GB) | 440.84 s / 5.2 GB | 5.4% | 130.04 s / 3.1 GB |

"Critical cells" is the share of cells the gradient leaves unmatched. In every case above it equals the minimum any filtered gradient can reach, 2·(finite bars) + (essential bars). The engine therefore reduces 0.1–0.2% of the cells on smooth data and 5–9% on white noise. In 2D the same holds: 0.33% on smoothed noise and 15.5% on white noise. Apparent pairs alone leave 1.7–50 times more (`docs/TECHNIQUES.md` §4.6).

* In 3D and 4D the gradient engine is 1.8–3.0× faster than CubicalRipser, 2.2–3.9× faster than the classic engine and 4.7–15× faster than GUDHI. It uses less memory than any of them, even though its time includes process start and reading the file. On white noise 48⁴ it is the only one of the three implicit engines to finish within this machine's memory.
* Q certified (default `--mode auto`: exact int64 arithmetic, a proof) costs 1.0–1.6× F_2 on these inputs, and 2× on white noise 48⁴.
* Where the time goes (white noise 48⁴, F_2, 60 s): gradient 6.4 s (1.2 million distinct block patterns, 259 exact searches), Morse complex 15.5 s (101 million path cells visited), reduction 35 s. Degree 1 alone took 29 s: its long working columns were binary heaps, which keep duplicate rows until they surface (the longest held 2.4 million entries). Two remedies for degree 1 were tried, and both are slower. Storing the reductions (V) instead of the reduced columns, as Ripser does, merged 5–9× more entries and took 10× longer. Reducing degree 1 by homology took 83 s. On smooth inputs the reduction is negligible (under 0.1 s), and the gradient and the path sums dominate.
* **Bit-tree working columns (19 Nov 2026).** Long working columns now live in a 64-ary bit tree with a coefficient array (`bitcol.hpp`, after PHAT), where duplicates cancel on arrival. This applies to the gradient and classic engines and to every coefficient ring. Measured on 19 Nov with identical output:
  * white noise 48⁴, F_2: degree-1 reduction from 26.3 to 7.9 s, total from 66.1 to 52.0 s;
  * white noise 48⁴, Q certified: total from 134.5 to 76.6 s, with degree 1 from 48.6 to 17.2 s and the reruns mod 2 and 3 from 27 and 34 s to 7.5 and 10.7 s;
  * white noise 40⁴, F_2: degree 1 from 6.5 to 2.4 s;
  * white noise 32⁴, classic engine, F_2: total from 20.3 to 12.4 s;
  * smoothed noise 48⁴: total from 8.4 to 7.2 s.

  The table above predates this change.

The table this replaces measured the classic engine on other inputs: smoothed noise 4096² 8.0 s, white noise 4096² 11.6 s, smoothed noise 192³ 10.4 s, white noise 192³ 26.0 s, white noise 32⁴ 19.3 s, smoothed noise 48⁴ 20.3 s, against CubicalRipser's 7.4, 8.5, 9.3, 24.9, 15.2 and 16.0 s.

The tables below were measured on an earlier machine with the original 6.0 engines and are kept for reference.

For degree 2 and up, Perseus 6 first runs its own iterated edge collapse (Boissonnat–Pritam, in the shift-and-trim form of Glisse–Pritam), implemented with bitsets. It keeps the collapse only when it actually removes edges. Against GUDHI's edge collapse followed by Ripser, the strongest Rips pipeline we know of (`bench/vs_collapse.py`, identical barcodes):

| case | edges before → after collapse | collapse + Ripser | Perseus 6 |
|---|---|---|---|
| o3_1024 H≤3 | 33,189 → 7,179 | 0.20 s | 0.06 s |
| o3_4096 H≤3 | 234,145 → 29,724 | 1.43 s | 0.72 s |
| random16 H≤7 | 942 → 331 | 0.02 s | 0.01 s |
| fractal H≤2 | 122,835 → 27,871 | 6.22 s | 1.21 s |
| sphere H≤2 | 17,857 → 17,668 | 1.84 s | 0.50 s |
| dragon H≤1 | 1,388,645 → 12,034 | 22.2 s | 0.64 s |

Against plain Ripser and RedZeD:

| case | Perseus 6, Q certified | Perseus 6, F_2 | Ripser 1.2 (F_2) | RedZeD (F_2) |
|---|---|---|---|---|
| sphere H≤2 | 0.52 s / 109 MB | 0.50 s / 83 MB | 0.79 s / 103 MB | 0.26 s / 16 MB |
| o3_1024 H≤3 | 0.10 s / 19 MB | 0.09 s / 18 MB | 2.57 s / 22 MB | 0.90 s / 110 MB |
| o3_4096 H≤1 | 0.59 s / 117 MB | 0.50 s / 116 MB | 1.32 s / 68 MB | 0.32 s / 51 MB |
| o3_4096 H≤3 | 0.86 s / 208 MB | 0.74 s / 208 MB | 50.32 s / 244 MB | 29.22 s / 1337 MB |
| random16 H≤1 | 0.03 s / 13 MB | 0.03 s / 13 MB | 0.03 s / 13 MB | 0.03 s / 13 MB |
| random16 H≤7 | 0.03 s / 13 MB | 0.03 s / 13 MB | 4.45 s / 124 MB | 5.91 s / 746 MB |
| fractal H≤1 | 0.10 s / 13 MB | 0.10 s / 13 MB | 0.20 s / 21 MB | 0.07 s / 13 MB |
| fractal H≤2 | 1.21 s / 13 MB | 1.22 s / 13 MB | 5.75 s / 24 MB | 2.46 s / 36 MB |
| dragon H≤1 | 0.63 s / 100 MB | 0.63 s / 99 MB | 1.70 s / 138 MB | 0.21 s / 29 MB |
| randdist180 H≤1 | 0.06 s / 13 MB | 0.04 s / 13 MB | 0.66 s / 56 MB | 0.04 s / 13 MB |

Against GUDHI on simplicial complexes (`bench/vs_gudhi.py`, F_2, identical intervals in every case; for cubical see the table above):

| case | GUDHI | Perseus 6 |
|---|---|---|
| Linial–Meshulam, n=70 | 0.35 s | 0.27 s |
| Linial–Meshulam, n=100 | 4.7 s | 2.9 s |

On a 1000×1000 image in the legacy `cubtop` text format, Perseus 6 took 2.2 s over F_2 (2.8 s certified over Q) and 470 MB, against 21.9 s and 1.8 GB for Perseus 5, with identical intervals. That was measured with the classic engine, before the planar engine existed, and most of it is reading text. For the planar engine, see the 2D table above.

Random metrics have many bad primes. For 450 uniformly random distances (H≤1):

* the Q barcode takes 3.2 s;
* the complete certificate covers 179 non-unit pivots and 225 candidate primes, of which 161 are bad, the largest ≈ 2.6·10¹³;
* the certificate takes 288 s on 2 cores. Classification is embarrassingly parallel, so `--threads` makes it scale with the core count.

## Old Perseus files

The Perseus 4/5 input formats (`simtop`, `nmfsimtop`, `cubtop`, `scubtop`, `cellcomp`, `distmat`, `corrmat`, `rips`, `brips`) are not read directly. Convert them:

```
python3 tools/convert_legacy.py cubtop old.txt new.txt
# prints e.g.:  perseus cubical new.txt --coeff 2 --periodic 1
```

The converter reproduces the old semantics exactly, and the test suite checks converted files against the Perseus 5 executable. Those semantics include faces inheriting the oldest coface, the +1 birth shift for cubical grids containing a 0, the first occurrence of a repeated cube winning, step-quantised distances, and growing balls. The printed command uses `--coeff 2`, because Perseus 4/5 computed over F_2; drop it to get the rational barcode and its bad primes.

Two things differ from the old programs:

* For Rips-type files, the old programs built simplices up to a dimension cap C and reported spurious infinite intervals in degree C. The printed command uses `--maxdim C−1`, the degrees the old output got right.
* Betti files list Betti numbers at every value where they change.

## Tests

`tests/run_tests.py` runs differential checks against four oracles:

* **A:** Ripser, exact interval multisets over F_2, F_3, F_5, F_7, F_127 and Q.
* **B:** the Rips engines against the explicit engine on random metrics with ties: barcodes and bad primes.
* **C:** an independent Python reference (Z-elimination then field reductions) on Linial–Meshulam, random simplicial and torsion-rich chain complexes.
* **D:** old Perseus files of every kind, converted by `tools/convert_legacy.py`, against the Perseus 5 binary.
* **E:** the new formats against independent constructions (`flag` against explicit cliques and against `lower-distance`, `cubes` against `cubical`).
* **F:** certificate regressions with exact expected answers in auto, fast and proved mode. These include inputs built so that a pivot vanishes or looks like a unit modulo fixed primes, a 98-bit bad prime, and products of large primes that are only classified as groups.
* **Python:** `tests/test_python.py` runs every entry point of the module against the command line, in every certificate mode.
* **G:** the implicit cubical engines against the explicit complex and against each other (1–4 dimensions, periodic axes, absent cells, zero-length bars, 4D white noise and quantised fields), and the sparse Rips route against the dense engine.
* **Cubical engines:** `make test-cubical` (`tests/t_cubical_morse.cpp`) runs the discrete-gradient engine against the classic one on random grids with heavy ties, absent cells and nearly distinct values. It compares barcodes over F_2, F_3, F_5 and Q (int64 and multi-modular), with and without zero-length bars, and from every degree with clearing data. It also compares certificates in auto and fast mode. It checks the Morse complex itself: every column is recomputed by recursion over gradient paths, along with the filtration order, δδ = 0, and the apparent pairs. It also covers degenerate shapes (axes of length 0, 1 and 2, very elongated grids, a single cube, grids entirely absent or constant) and the int64 guard: with the bound on path sums lowered, the engine must refuse, and the default must fall back to the classic engine.
* **Fuzzing:** `tests/fuzz_cubical.py` drives the Python module against the explicit complex (F_2, F_3, and Q with bad primes, with absent cells) and against CubicalRipser (F_2). It uses value distributions built for ties: few integer levels, quantised smooth fields, step functions, ±0, extreme magnitudes, plus degenerate shapes.

The current tree passes all of them.

## Limitations and roadmap

* Dense Rips inputs (distance matrices) are still n × n; point clouds and graphs go through the sparse route when it pays.
* Classifying many candidate primes dominates on torsion-heavy inputs (Linial–Meshulam at the threshold: about 100 s for 359 candidates at n = 90). Settling groups of primes in one run modulo their product is the obvious next step.
* Twisted coefficients (Z[t]/(t^n − 1) via primes p ≡ 1 mod n) for Milnor fibres and Novikov homology are planned.
* Checkpointed reruns for candidate primes.

## Techniques and credits

`docs/TECHNIQUES.md` catalogues every technique used, with its origin, and says plainly what is new here and what is not.

## Credits

Vidit Nanda (Perseus, 2010–). Perseus 6 builds on ideas from:

* U. Bauer, *Ripser*;
* K. Mischaikow and V. Nanda, *Morse theory for filtrations and efficient computation of persistent homology* (the filtered discrete gradient of the cubical engine);
* V. Robins, P. J. Wood and A. P. Sheppard, *ProcessLowerStars* (discrete gradients of grayscale images);
* P. Guillou, J. Vidal and J. Tierny, *Discrete Morse Sandwich* (union-find on critical cells for the extreme degrees);
* S. Kaji, T. Sudo and K. Ahara, *Cubical Ripser* (the implicit cubical design of the classic engine);
* J.-D. Boissonnat and S. Pritam, edge collapse and edge domination;
* *RedZeD*;
* I. Obayashi and M. Yoshiwaki, field choice and torsion;
* J.-D. Boissonnat and C. Maria, multi-field persistence and modular reconstruction.

## Citing

Vidit Nanda, *Perseus, the persistent homology software*, version 6.0.0, [doi:10.5281/zenodo.24305983](https://doi.org/10.5281/zenodo.24305983). GitHub's "Cite this repository" button gives the same from `CITATION.cff`.
