# Perseus 6

Persistent homology **over Q, with a certificate of every prime p for which the F_p barcode is different**, or over any prime field.

One run gives you:

* the rational barcode,
* the complete list of *bad primes* — the p for which F_p coefficients give a different barcode — and, for each, exactly which intervals change,
* the exact rational values of all non-unit pivots, which explain where the torsion lives.

Beyond simplicial and cubical complexes, Perseus 6 takes filtered simplicial sets, nerves of finite categories and arbitrary filtered chain complexes over Z, and returns integral information: the torsion of the homology, the image of H_*(K_s; Z) → H_*(K_t; Z), the F_p barcode of every bad prime, and the torsion of H_q(K_t; Z) at every step of a filtration (for instance along a Linial–Meshulam process; see [Torsion along a filtration](#torsion-along-a-filtration)). It computes magnitude homology of graphs over Z, torsion included, and handles chain complexes over the Laurent ring Z[t, t^-1] (infinite cyclic covers, Novikov homology). See [Chain complexes over Z and beyond](#chain-complexes-over-z-and-beyond).

The Q reduction is exact, so the certificate is a proof: int64 arithmetic while it fits, at about the cost of an F_p run, and exact big rationals in degrees where it overflows, as on torsion-heavy inputs. `--mode fast` runs modulo random 62-bit primes instead and reports a rigorous bound on the probability of error; it is not a proof, and says so. See [docs/CERTIFICATE.md](docs/CERTIFICATE.md).

Perseus 6 is a rewrite from scratch. It does not read the old Perseus 4/5 input formats directly; `tools/convert_legacy.py` turns old files into new ones (see [Old Perseus files](#old-perseus-files)).

## Install

Executables for Linux, macOS and Windows are attached to each [GitHub release](https://github.com/viditnanda/perseus/releases); the Python module is `pip install perseus-tda` (Linux and macOS).

## Build

```
make                 # builds ./perseus (C++17 with 128-bit integers: g++ or clang; MinGW-w64 on Windows)
make test            # differential tests (needs Ripser / Perseus 5 / p2ref; skips what is missing)
make test-chain      # chain complexes over Z, simplicial sets, magnitude homology, Laurent coefficients
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

1. **Reduces over Q, exactly.** By default (`--mode proved`; `--mode auto` is the same) the reduction runs in exact rational arithmetic: int64 numerators and denominators with every operation overflow-checked, and, in a degree where they overflow, rationals that stay in int64 while they fit and become big integers only where they must. The barcode and the pivot values are exact, so the certificate is a proof. `--mode fast` reduces modulo K random primes from [2^61, 2^62) instead and reports a bound on the probability of error that holds for every input (default below 2^-64; `--target`). It is not a proof, and every output says so: a warning on stderr, `"proof": false` and a `"warning"` entry in the JSON, and a Python warning. The bound is over a fresh 256-bit seed; `--seed` reproduces a run but does not certify it.
2. **Reconstructs** non-unit pivot values in fast mode, by Wang reconstruction checked against two independent primes.
3. **Factors** the pivot values: trial division, division by primes already found, Pollard–Brent rho, and gcds between leftovers. Large primes and unsplit composites N are not factored further.
4. **Settles every candidate.** Primes below 2^63 are rerun over F_p, in parallel. On explicit complexes (simplices, cells, simplicial sets) they can instead be packed into 62-bit moduli and settled by one sweep per modulus over Z/m in Howell form, which gives every packed prime's F_p barcode at once. By default a timed sample of reruns decides which is cheaper for the input at hand (`--classify-by auto|sweep|rerun`). Each remaining N is settled by one run modulo N. Either the run splits N, or every prime factor of N gives the same barcode. Then the certificate states the answer for all of them at once, e.g. "every prime factor of N is bad", without knowing the factors.

The output states what is guaranteed about the barcode (`proved`, or `probabilistic` with the bound), and whether the description of the bad primes is complete. Field dependence of persistence is governed by torsion in the relative homology H_q(X_n, X_m; Z) (Obayashi–Yoshiwaki): p is bad in degree q exactly when some H_q(K_j, K_i; Z) or H_{q−1}(K_i; Z) has p-torsion (`docs/CERTIFICATE.md` §7, Theorem 4). The pivots above are where that torsion is detected.

## Torsion along a filtration

For explicit complexes, `--torsion` (Python: `torsion=True`) adds the torsion of H_q(K_t; Z) at every value t where it changes: its invariant factors and, for each prime, the exponents of its p-primary part. Nothing is sampled. The trajectory comes from sweeps of the boundary matrices modulo p, p^2, … for the bad primes, batched into 62-bit words, with big integers for larger primes and for composites that could not be factored (`docs/CERTIFICATE.md` §7). This is what studies of torsion in random complexes look at (Kahle–Lutz–Newman–Parsons), at every step of the process. On a Linial–Meshulam process near the threshold (`bench/certify_lm.py`, n = 90; all edges at time 1, then triangles one by one):

```
$ perseus simplices lm90_4000.txt --torsion --print none
...
torsion of H_1(K_t; Z), where it changes:
  t = 3564: Z/2
  t = 3565: Z/109384485698954982065478566306368771057679577324100
  t = 3568: 0
```

The big cyclic group is 2^2 · 5^2 · 479 · 22303 · 1157339 · p with p a 117-bit prime. These are exactly the bad primes, as Theorem 4 predicts when the whole 1-skeleton comes first. d-dimensional processes work the same way (complete (d−1)-skeleton, then d-simplices; torsion in H_{d−1}), as does any filtration of a simplicial complex, simplicial set or chain complex over Z, with ties or without. The results were checked against FLINT Smith forms at every t on small processes and random chain complexes, and at 194 values in the torsion windows of the two benchmark processes.

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
| `cells` | filtered based chain complex over Z: one cell per line `dim birth [face coef]...`, faces are 0-based line numbers of earlier cells; d∘d = 0 is checked exactly |
| `sset` | finite filtered simplicial set: one nondegenerate simplex per line `dim birth d_0 ... d_dim`; a face is an earlier line `y`, or `y@k1.k2...km` for the degenerate simplex s_k1 s_k2 ... s_km y |
| `laurent` | filtered based chain complex over Z[t, t^-1]: one cell per line `dim birth face:poly ...`, e.g. `2 0 1:t^2-t+1 2:-1+t^-1` |
| `magnitude` | a graph, one edge `u v` per line (vertices 0, 1, ...), for magnitude homology |

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

For `cells`, `sset` and `simplices`, `--integral` adds `integral_homology` (free rank and torsion of H_k of the whole complex, by degree) and, in Q mode, `bad_prime_barcodes` (the full F_p barcode of each bad prime below 2^64); `--image s,t` adds `integral_image`, the image of H_k(K_s; Z) → H_k(K_t; Z); `--torsion` adds `torsion`: per degree, the steps `{"value", "invariant_factors", "primary": {p: exponents}}` where tors H_q(K_t; Z) changes, `groups` (composites standing for all their prime factors) and `proof`. The text forms go to stderr.

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

Chain complexes over Z, simplicial sets, nerves, magnitude homology and Laurent coefficients:

```python
perseus.cells(cells, integral=True, image=(0, 1))     # r["integral"]: homology, image, bad_prime_barcodes
perseus.simplices(S, torsion=True)                    # r["torsion"]["degrees"][q]: steps of tors H_q(K_t; Z)
perseus.simplicial_set([(0, 0, []), (1, 0, [0, 0]), (2, 1, [1, (0, [0]), 1])], integral=True, image=(0, 1))   # RP^2
perseus.nerve(src, tgt, identities, composition, weights=w, maxdim=3, integral=True)
perseus.bar_construction(mult, identity=0, maxdim=5)   # H_k(G; Z) for a finite group G
perseus.magnitude(edges, lmax=6)                       # r["table"][(k, l)] = {"free", "torsion": {order: multiplicity}, ...}
perseus.laurent([(0, 0, []), (1, 0, [(0, "t - 1")])])  # barcode over Q(t), bad irreducibles, homology over Q[t^±1]
```

Every function takes `coeff`, `mode`, `seed`, `target`, `classify`, `threads` and `fp`, as on the command line; `simplices`, `cells` and `simplicial_set` also take `torsion` and `classify_by` (`auto`, `sweep`, `rerun`); `rips` and `flag` also take `engine` (`auto`, `implicit`, `sparse`) and `collapse`, and `cubical` takes `engine` (`auto`, `planar`, `morse`, `classic`, `explicit`). For the planar engine (the 2D default) the array is read in place, without a copy. A key of `certificate["bad_primes"]` whose entry has `every_prime_factor_of` set is a number all of whose prime factors are bad; `probable_prime` marks a prime above 2^64 known only by Miller–Rabin. `tests/test_python.py` checks the module against the command line.

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

## Chain complexes over Z and beyond

**Filtered chain complexes over Z** (`cells`). Any based chain complex with integer boundary coefficients and a filtration by births. Perseus checks d∘d = 0 exactly, then cancels every pair σ < τ with [τ : σ] = ±1 and equal births (algebraic Morse reduction). That is a filtered chain homotopy equivalence over Z, so every barcode, the bad primes and every integral group below are unchanged; `--no-reduce` skips it. The Q barcode and its bad-prime certificate are computed as for any complex. With `--integral` and `--image s,t`:

* the integral homology of the whole complex, H_k(K_∞; Z) = Z^b ⊕ ⊕ Z/d_i, from Smith forms of the reduced boundary matrices;
* the image of H_k(K_s; Z) → H_k(K_t; Z), computed as (Z_k(K_s) + B_k(K_t)) / B_k(K_t). Two-stage filtrations are the main use. For RP^2 built on a projective line (K_0 = RP^1, K_1 = RP^2) the image in degree 1 is Z/2: the generator of H_1(RP^1) = Z maps onto H_1(RP^2) = Z/2;
* the F_p barcode of every bad prime.

All of it is exact integer arithmetic (int64 with overflow checks, then bigints). The Smith forms are dense, so they run on the residual complex after reduction.

**Simplicial sets** (`sset`, `perseus.simplicial_set`). Nondegenerate simplices with their faces, which may be degenerate; degeneracies are normalised away (normalised chains). The simplicial identities are checked. Constructions in the library: nerves of finite categories (`perseus.nerve`), the bar construction BG of a finite group or monoid (`perseus.bar_construction`), quotients X/A, cones, suspensions, and minimal spheres. In a **weighted nerve** a chain (f_1, …, f_n) of composable morphisms is born at the largest weight of a composite f_j ∘ … ∘ f_i of consecutive morphisms, so faces are never born after their simplices.

**Magnitude homology** (`magnitude`, `perseus.magnitude`). MH_{k,l}(G; Z) of a graph (Hepworth–Willerton), directed or not, for lengths l ≤ `--lmax`, with torsion. The complex splits by endpoints (s, t) and length l. Each piece is enumerated in lexicographic order, stored compactly, reduced by the fill-in-free part of the algebraic Morse reduction in place, then by the general reduction and Smith forms. `--pairs s,t,...` restricts to chosen endpoints. The answer is exact; nothing is sampled or reduced mod p. [gramag](https://gramag.readthedocs.io/) computes ranks over F_2 only. On every graph of `bench/vs_gramag.py` where it finishes, they equal Perseus's integral answer reduced by universal coefficients (see [Benchmarks](#magnitude-homology-measured-18-nov-2026-2-core-xeon)).

**Laurent coefficients** (`laurent`, `perseus.laurent`). Chain complexes over Z[t, t^-1]: infinite cyclic covers, Alexander modules, Novikov homology. Z[t, t^-1] is not a principal ideal domain, so its homology has no barcode and no Smith form. Perseus computes what the principal ideal domains around it see, and says which:

* the barcode over the field Q(t). It equals the barcode over the Novikov field Q((t)), and its infinite bars count the Novikov Betti numbers;
* the **bad irreducibles**: the irreducible f ∈ Q[t] (f ≠ t) whose residue field Q[t]/(f) gives a different barcode, each with that barcode. When a polynomial is not proved irreducible, the answer holds for every irreducible factor of it;
* the **bad primes** p for the barcode over F_p(t);
* the homology of the whole complex over Q[t, t^-1]: free rank and invariant factors (the Alexander polynomials, up to units ±t^k);
* the primes p for which the homology over F_p[t, t^-1] is not the reduction mod p of the rational answer, each with its homology over F_p[t, t^-1]. This is the integral content Perseus reports. For Z[t, t^-1]/(2, t + 1), which is zero over Q[t, t^-1], it gives F_2[t, t^-1]/(t + 1) at p = 2.

Not computed: the Z[t, t^-1]-module itself, residue fields F_p[t]/(h), and torsion over the Novikov ring Z((t)). Pairs with unit incidence ±t^k are cancelled first. `docs/CERTIFICATE.md` §6 states what is certified.

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

On geometric data a proof costs nothing extra; time beyond F_2 is classification of candidate primes. Exact arithmetic gets expensive only when the integers involved grow, which happens on torsion-heavy inputs, and then the default switches to random primes. On Linial–Meshulam complexes the Q reduction itself is fast (1.3 s at n = 90, with 261 non-unit pivots of up to 170 bits); classifying the 463 candidate primes, of which exactly two are bad, takes the rest. Since 19 Nov 2026 the candidates of at least 2^63 are classified in parallel too. On the reproducible instance of `bench/certify_lm.py` (n = 90, the first 4,000 triangles; 496 candidates, 6 bad), the default run went from 147 to 105 s on 2 cores, and the reruns scale with cores.

**Since 19 Nov 2026 the default is proved everywhere.** The table above is from 16 Nov, when the default (`auto`) fell back to random primes where int64 overflowed, as in its Linial–Meshulam rows. Now the default never falls back. Its exact arithmetic was made fast for this:
* rationals stay in int64 pairs while they fit;
* big integers keep up to 256 bits inline;
* gcds finish as binary gcds on 128-bit words;
* sums and products cancel common factors before multiplying (Henrici).

On the reproducible instances of `bench/certify_lm.py` (measured 19 Nov 2026, 2 cores, raw results in `bench/certify_lm.json`). Since 19 Nov the small candidates of explicit complexes can be classified by batched Howell sweeps. The table compares forced sweeps with `--classify-by rerun`, which gives identical results:

| instance | default, proved | `--classify-by rerun` | small candidates: sweeps | small candidates: reruns | `--torsion` adds | `--mode fast` |
|---|---|---|---|---|---|---|
| n = 70 | 4.3 s | 7.0 s | 2.1 s (61 sweeps) | 4.8 s (202 reruns) | 0.3 s | 2.9 s |
| n = 90 | 93.8 s | 125.5 s | 24.4 s (124 sweeps) | 56.4 s (386 reruns) | 19.9 s | 74.8 s |
| n = 70, distinct values | 18.7 s | 31.4 s | 7.3 s (114 sweeps) | 20.1 s (326 reruns) | 21.7 s | 11.2 s |

The exact Q run alone takes 1.8 s and 23.9 s on the two tied instances (it took 202 s at n = 90 with plain big rationals). What remains at n = 90 is mostly the 110 candidates of at least 2^63 (39 s), still one run each. In the third instance edges and triangles have distinct interleaved values; there all 347 candidates are bad. Both modes and both classifications find the same bad primes.

The default, `--classify-by auto`, chooses between the two per input. It reruns a timed sample of candidates first, and sweeps only if reference sweeps run under that budget predict a saving. On `bench/classify_auto.py` (20 Nov 2026, 2 cores; seconds spent classifying the small candidates; identical bad primes in all three columns):

| instance | candidates | reruns | sweeps | auto |
|---|---|---|---|---|
| LM n = 90, tied | 496 | 54.2 | 24.2 | 24.0 (sweeps) |
| LM n = 70, distinct values | 347 | 18.8 | 6.8 | 6.9 (sweeps) |
| LM d = 3, n = 25 | 333 | 74.5 | 29.8 | 31.7 (sweeps) |
| RP² Rips, 300 points, sparse route | 1 | 0.5 | 44.9 | 0.3 (reruns) |
| Rips with 0.9M tetrahedra + LM n = 40 | 6 | 0.17 | 1.69 | 0.34 (reruns, sweeps cut off) |
| Rips with 0.9M tetrahedra + LM n = 70 | 347 | 111.0 | 18.4 | 18.3 (sweeps) |

### Cubical (F_2, measured 19 Nov 2026 in one sitting, 2-core Xeon, 7.8 GB)

These two tables are the single source of truth for cubical timings. Both were measured on the same code, with the planar engine and the bit-tree working columns, in one sitting. `bench/cubical2d.py` and `bench/cubical.py` regenerate the inputs and run every program in a fresh process. The raw results are in `bench/cubical2d.json` and `bench/cubical.json`.
* Python tools (Perseus's module, GUDHI, CubicalRipser) are timed around the in-process call, after their imports, with the input array already in memory. Their memory is the process's peak RSS read right after the call, so it includes the interpreter, the modules and the input array.
* The Perseus command line is timed from outside, including process start and reading the `.npy` file, with output suppressed; its memory is the peak RSS.
* Each time is the best of three runs, except GUDHI and CubicalRipser in 3D and 4D, which ran once.
* The positive-length intervals agree in every case between every pair of programs that finished.

**2D images.** The inputs are the 512² camera photograph from scikit-image (8-bit grey levels, so heavy ties), uniform noise, and uniform noise smoothed by a Gaussian of width 4, all float64. GUDHI's `CubicalPersistence` is asked for degrees 0 and 1, which selects its dedicated image algorithm. At 8192² the gradient engine and CubicalRipser need more memory than this machine has.

| image | Perseus 6, planar engine | planar, Q certified | Perseus command line | GUDHI 3.13 | CubicalRipser | Perseus 6, gradient engine | intervals |
|---|---|---|---|---|---|---|---|
| photograph 512² | **9.8 ms / 40 MB** | 11.9 ms / 41 MB | 14.3 ms / 14 MB | 13.1 ms / 147 MB | 64.1 ms / 101 MB | 38.9 ms / 48 MB | 36,849 |
| smoothed noise 1024² | **11.9 ms / 53 MB** | 14.9 ms / 54 MB | 27.5 ms / 33 MB | 21.8 ms / 162 MB | 0.38 s / 264 MB | 0.19 s / 112 MB | 7,022 |
| white noise 1024² | **68.4 ms / 67 MB** | 69.0 ms / 68 MB | 77.7 ms / 39 MB | 85.8 ms / 175 MB | 0.49 s / 302 MB | 0.31 s / 122 MB | 325,501 |
| smoothed noise 4096² | **0.23 s / 253 MB** | 0.26 s / 253 MB | 0.41 s / 260 MB | 0.55 s / 499 MB | 8.54 s / 3.0 GB | 3.62 s / 1.0 GB | 110,329 |
| white noise 4096² | **1.60 s / 451 MB** | 1.65 s / 451 MB | 1.79 s / 424 MB | 1.88 s / 682 MB | 10.12 s / 3.0 GB | 7.48 s / 1.4 GB | 5,216,867 |
| smoothed noise 8192² | **1.06 s / 893 MB** | 0.92 s / 893 MB | 1.98 s / 1.0 GB | 2.56 s / 1.5 GB | out of memory | out of memory | 441,508 |
| white noise 8192² | **7.51 s / 1.6 GB** | 7.88 s / 1.6 GB | 8.08 s / 1.6 GB | 9.85 s / 2.3 GB | out of memory | out of memory | 20,876,216 |

* The planar engine is 1.2–2.4× faster than GUDHI's image algorithm and uses 27–73% less memory: 1.2–1.3× on white noise, 1.8–2.4× on smoothed noise, 1.3× on the photograph. It is 6–37× faster than CubicalRipser.
* The certificate costs nothing: the boundary matrices it reduces are totally unimodular, so the Q barcode is proved equal to every F_p barcode without arithmetic (`docs/CERTIFICATE.md` §5).
* Where the time goes (white noise 4096², 1.25 s in C++): classification 0.13 s, sorting 6.4 million births and candidates 0.22 s, degree 1 0.59 s, degree 0 0.28 s. Only 38% of the pixels are sorted; on smoothed noise, under 1%.

**3D and 4D grids.** The inputs are uniform noise, and noise smoothed by a Gaussian of width 2 in 3D and 1.5 in 4D, all float64. All three programs compute the T-construction. On white noise 48⁴, where CubicalRipser did not finish, the two Perseus engines agree.

| case | Perseus 6, gradient engine | Perseus 6, classic engine | CubicalRipser (`tcripser`) | GUDHI 3.13 (`CubicalPersistence`) | critical cells | gradient engine, Q certified |
|---|---|---|---|---|---|---|
| smoothed noise 96³ | **0.37 s / 96 MB** | 0.98 s / 142 MB | 1.12 s / 132 MB | 4.80 s / 470 MB | 0.22% | 0.36 s / 96 MB |
| white noise 96³ | **1.02 s / 134 MB** | 2.30 s / 179 MB | 2.64 s / 272 MB | 8.93 s / 553 MB | 8.7% | 1.27 s / 150 MB |
| smoothed noise 24⁴ | **0.38 s / 63 MB** | 0.86 s / 126 MB | 0.79 s / 121 MB | 4.44 s / 396 MB | 0.13% | 0.39 s / 63 MB |
| white noise 24⁴ | **1.40 s / 107 MB** | 2.59 s / 166 MB | 3.24 s / 274 MB | 7.79 s / 445 MB | 4.7% | 1.68 s / 120 MB |
| smoothed noise 192³ | **3.41 s / 764 MB** | 12.40 s / 1.1 GB | 10.68 s / 800 MB | 52.81 s / 2.7 GB | 0.23% | 3.37 s / 764 MB |
| white noise 192³ | **11.06 s / 1.0 GB** | 26.38 s / 1.4 GB | 26.40 s / 1.7 GB | 98.42 s / 3.4 GB | 9.0% | 15.99 s / 1.3 GB |
| smoothed noise 32⁴ | **1.27 s / 183 MB** | 3.49 s / 386 MB | 3.02 s / 276 MB | 15.06 s / 923 MB | 0.13% | 1.29 s / 183 MB |
| white noise 32⁴ | **5.48 s / 340 MB** | 11.19 s / 562 MB | 14.54 s / 903 MB | 36.23 s / 1.1 GB | 5.0% | 7.11 s / 404 MB |
| smoothed noise 48⁴ | **7.35 s / 923 MB** | 24.26 s / 1.7 GB | 18.83 s / 1.2 GB | 88.03 s / 4.0 GB | 0.14% | 7.45 s / 923 MB |
| white noise 48⁴ | **37.52 s / 1.6 GB** | 89.99 s / 3.5 GB | killed (out of memory above 6 GB) | 435.84 s / 5.2 GB | 5.4% | 65.53 s / 2.8 GB |

"Critical cells" is the share of cells the gradient leaves unmatched. In every case above it equals the minimum any filtered gradient can reach, 2·(finite bars) + (essential bars). The engine therefore reduces 0.1–0.2% of the cells on smooth data and 5–9% on white noise; in 2D, 0.33% on smoothed noise and 15.5% on white noise. Apparent pairs alone leave 1.7–50 times more (`docs/TECHNIQUES.md` §4.6).

* In 3D and 4D the gradient engine is 2.1–3.1× faster than CubicalRipser, 1.8–3.6× faster than the classic engine and 5.6–15× faster than GUDHI. It uses less memory than any of them, even though its time includes process start and reading the file. On white noise 48⁴ it is the only one of the three implicit engines to finish within this machine's memory.
* Q certified (the default, a proof: exact int64 arithmetic suffices on these inputs) costs 1.0–1.5× F_2, and 1.75× on white noise 48⁴.
* Where the time goes (white noise 48⁴, F_2, 37.5 s):
  * setup 24.4 s: births 1.2 s; gradient 6.8 s (1.2 million distinct block patterns, 259 exact searches); critical cells 1.0 s; the Morse complex 14.8 s (101 million path cells visited);
  * reduction 12.8 s: degree 1 7.9 s, degree 2 4.9 s.

  The reduction used to dominate. Degree 1 alone took 26–29 s while its long working columns were binary heaps, which keep duplicate rows until they surface. They are now bit trees with a coefficient array (`bitcol.hpp`, after PHAT; `docs/TECHNIQUES.md` §1.13), where duplicates cancel on arrival, for every coefficient ring. Before that, two other remedies were tried, and both were slower. Storing the reductions (V) instead of the reduced columns, as Ripser does, merged 5–9× more entries and took 10× longer. Reducing degree 1 by homology took 83 s. On smooth inputs the reduction is negligible, and the gradient and the path sums dominate.

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

### Magnitude homology (measured 18 Nov 2026, 2-core Xeon)

**Torsion families** (`bench/mag_families.cpp`, results in `bench/mag_families.json`). The Kaneta–Yoshinaga graph of a regular CW complex M is the Hasse diagram of its face poset with a bottom 0̂ and a top 1̂ added. H̃_{k−2}(M) embeds in MH_{k,d}(0̂, 1̂), d = d(0̂, 1̂) = dim M + 2. Sazdanovic and Summers use these graphs, for projective spaces and lens spaces, to get torsion of every prime order. Cases: `rpN` is RP^N as the boundary of the (N+1)-cube modulo the antipodal map; `rp2_6` is the 6-vertex RP^2; `lensP_M` is L(P, 1), of dimension 2M − 1, as a join of M cycles C_2P modulo a diagonal rotation. Each row is the largest l reached with one thread and at most 3.5 GB of address space. The next even l ran out of memory in every (0̂, 1̂) row; each step of 2 multiplies the number of tuples by about 30 to 75. The all-pairs rows stop at the 900 s limit. Times are CPU seconds. Every group is computed over Z: the torsion below is exact, and so are the free ranks, which are omitted here.

| graph | pairs | vertices | d | l ≤ | generators | largest complex | time | memory | torsion (k, l) |
|---|---|---|---|---|---|---|---|---|---|
| `rp2_6` | (0̂, 1̂) | 33 | 4 | 10 | 20,437,390 | 19,988,734 | 20.4 s | 1619 MB | (3,4): Z/2 |
| `rp2` | (0̂, 1̂) | 15 | 4 | 10 | 1,987,552 | 1,920,676 | 1.7 s | 209 MB | (3,4): Z/2; (5,6): Z/2; (7,8): Z/2; (8,10): Z/2^4; (9,10): Z/2 |
| `rp3` | (0̂, 1̂) | 42 | 5 | 9 | 3,647,549 | 3,587,888 | 5.5 s | 306 MB | (3,5): Z/2; (5,7): Z/2^8; (7,9): Z/2^17 |
| `rp4` | (0̂, 1̂) | 123 | 6 | 8 | 1,687,264 | 1,675,022 | 2.9 s | 141 MB | (3,6): Z/2; (5,6): Z/2; (5,8): Z/2^16 |
| `rp5` | (0̂, 1̂) | 366 | 7 | 7 | 211,929 | 211,929 | 0.7 s | 41 MB | (3,7): Z/2; (5,7): Z/2 |
| `lens3_2` | (0̂, 1̂) | 58 | 5 | 9 | 8,524,653 | 8,411,880 | 13.0 s | 729 MB | (3,5): Z/3; (5,7): Z/3 |
| `lens5_2` | (0̂, 1̂) | 90 | 5 | 9 | 27,389,309 | 27,123,536 | 33.5 s | 2020 MB | (3,5): Z/5; (5,7): Z/5^3 |
| `lens7_2` | (0̂, 1̂) | 122 | 5 | 7 | 481,141 | 478,212 | 0.4 s | 40 MB | (3,5): Z/7; (5,7): Z/7^5 |
| `lens2_3` | (0̂, 1̂) | 366 | 7 | 7 | 211,929 | 211,929 | 0.5 s | 41 MB | (3,7): Z/2; (5,7): Z/2 |
| `lens3_3` | (0̂, 1̂) | 734 | 7 | 7 | 474,937 | 474,937 | 2.1 s | 97 MB | (3,7): Z/3; (5,7): Z/3 |
| `rp2_6` | all | 33 | 4 | 9 | 807,438,559 | 2,512,138 | 651.7 s | 260 MB | (3,4): Z/2^2 |
| `rp2` | all | 15 | 4 | 9 | 62,032,711 | 507,431 | 41.5 s | 54 MB | (3,4): Z/2^2; (4,5): Z/2^8; (5,6): Z/2^14; (6,7): Z/2^14; (6,8): Z/2^2; (7,8): Z/2^15; (7,9): Z/2^22; (8,9): Z/2^14 |
| `rp3` | all | 42 | 5 | 8 | 435,632,452 | 677,916 | 320.9 s | 65 MB | (3,5): Z/2^2; (4,6): Z/2^16; (5,7): Z/2^48; (6,8): Z/2^91 |

"Generators" counts all tuples up to length l; "largest complex" is the largest single piece MC_{*,l}(s, t). In the (0̂, 1̂) rows, the torsion at l = d is the Kaneta–Yoshinaga embedding of H_odd(RP^N) = Z/2 and H_odd(L(P, 1)) = Z/P. The torsion at larger l, such as Z/2^17 in MH_{7,9} for RP^3, goes beyond it.

**Against gramag** (`bench/vs_gramag.py`, results in `bench/vs_gramag.json`). gramag 0.4.1 computes ranks over F_2: its documentation describes representatives as Z_2-linear combinations of paths, and its reduction is lophat's. On the Kaneta–Yoshinaga graph of RP^2 it reports rank 1 in both MH_{3,4} and MH_{4,4}, where the integral groups are Z/2 and 0. It reports no torsion, so the comparison is of F_2 ranks, from Perseus's integral groups by universal coefficients. Both programs run in a fresh process (Perseus with 2 threads, gramag on all cores) under a 3.5 GB limit:

| graph | l ≤ | generators | Perseus | gramag | F_2 ranks |
|---|---|---|---|---|---|
| C_7 | 7 | 13,797 | 0.13 s, 9 MB | 0.06 s, 17 MB | agree |
| C_12 | 10 | 702,432 | 0.30 s, 9 MB | 2.23 s, 140 MB | agree |
| directed C_6 | 8 | 1,464 | 0.00 s, 9 MB | 0.03 s, 14 MB | agree |
| grid 4x4 | 6 | 111,468 | 0.04 s, 9 MB | 0.31 s, 34 MB | agree |
| G(20, 0.2) | 5 | 294,562 | 0.09 s, 9 MB | 0.70 s, 64 MB | agree |
| directed G(16, 0.18) | 5 | 25,446 | 0.01 s, 9 MB | 0.08 s, 20 MB | agree |
| KY rp2 (all pairs) | 6 | 383,455 | 0.14 s, 9 MB | 1.01 s, 83 MB | agree |
| KY rp2_6 (all pairs) | 5 | 393,623 | 0.15 s, 9 MB | 1.29 s, 87 MB | agree |
| KY rp2 (pair (0,1)) | 8 | 66,876 | 0.04 s, 10 MB | 19.15 s, 1325 MB | agree |
| KY rp3 (pair (0,1)) | 7 | 59,661 | 0.04 s, 9 MB | out of memory | — |
| KY lens3_2 (pair (0,1)) | 7 | 112,773 | 0.09 s, 13 MB | out of memory | — |

gramag enumerates the paths between all pairs of vertices before computing anything, even when asked for one pair. On the single pair (0̂, 1̂) of the RP^3 and L(3, 1) graphs at l ≤ 7 it ran out of memory after 1–2 minutes; Perseus takes 0.1 s. Where both finish, the ranks agree everywhere. On the inputs with more than 100,000 generators Perseus is 7 to 9 times faster and uses a fraction of the memory.

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
* **Chain complexes, simplicial sets, magnitude, Laurent:** `make test-chain`. `tests/t_chain.cpp`: Smith forms (with the bigint fallback), integral homology and images on RP^2 and two-stage filtrations, minimal spheres, torus, Klein bottle, RP^2 as simplicial sets, BZ/2, BZ/3, BS_3 (H_3 = Z/6), nerves with non-identity endomorphisms, suspensions, weighted poset nerves against order complexes, and on random complexes that the algebraic Morse reduction preserves barcodes over F_2, F_3, F_5 and Q, bad primes, integral homology and every image. `tests/t_magnitude.cpp`: the integral answer reduced by universal coefficients against field reductions of the unreduced complexes (trees, cycles, random graphs, directed and not), and the Z/2 of the Kaneta–Yoshinaga graph of RP^2. `tests/t_laurent.cpp`: Alexander polynomials of knots by Fox calculus, Z[t, t^-1]/(2, t + 1), the cylinder covering the torus, and on random complexes the Q(t) and F_p(t) barcodes against the rank function, the bad irreducibles against evaluation at t = t0 through the explicit engine, and the homology over Q[t, t^-1] and F_p[t, t^-1] against evaluation. `tests/t_howell.cpp`: the Howell sweeps against plain F_p reductions (pairs) and brute-force enumeration of spans modulo p^k (lengths), in one-word and BigInt arithmetic. The Python suite checks torsion trajectories against a Smith form at every t on random chain complexes, and that sweeps and reruns classify alike.
* **Fuzzing:** `tests/fuzz_cubical.py` drives the Python module against the explicit complex (F_2, F_3, and Q with bad primes, with absent cells) and against CubicalRipser (F_2). It uses value distributions built for ties: few integer levels, quantised smooth fields, step functions, ±0, extreme magnitudes, plus degenerate shapes.

The current tree passes all of them.

## Limitations and roadmap

* Dense Rips inputs (distance matrices) are still n × n; point clouds and graphs go through the sparse route when it pays.
* Classifying many candidate primes dominates on torsion-heavy inputs. On explicit complexes the small candidates are settled in batches by Howell sweeps when that is cheaper (2.2–6× faster on the benchmarks above). Candidates of at least 2^63 are still one run each, and implicit Rips and cubical inputs still rerun every candidate.
* Laurent coefficients cover Q(t), Q[t, t^-1], its residue fields and F_p(t), F_p[t, t^-1], not Z[t, t^-1] itself; images over Q[t, t^-1] are not computed yet. Finite twisted coefficients (Z[t]/(t^n − 1), for Milnor fibres) are planned.
* Integral homology and images use dense Smith forms on the residual complex; inputs that do not reduce well are limited to residual blocks of about 40 million entries.

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
* J. Luo and G. Henselman-Petrusek, interval decomposition over a PID and field choice;
* J.-D. Boissonnat and C. Maria, multi-field persistence and modular reconstruction.

## Citing

Vidit Nanda, *Perseus, the persistent homology software*, version 6.0.0, [doi:10.5281/zenodo.24305983](https://doi.org/10.5281/zenodo.24305983). GitHub's "Cite this repository" button gives the same from `CITATION.cff`.
