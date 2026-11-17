# Techniques in Perseus 6, and where they come from

This document lists every technique Perseus 6 uses, with its origin. For each one it says what it is, where it lives in the code, whom it is due to, and what (if anything) we changed. The last section is our own assessment of what is new.

File references are to `include/perseus/`.

---

## 1. Persistence algorithms

**1.1 Matrix reduction and the pairing.** The standard column reduction of the filtered boundary matrix, and the elder rule.
Edelsbrunner, Letscher, Zomorodian, *Topological persistence and simplification*, DCG 2002; Zomorodian, Carlsson, *Computing persistent homology*, DCG 2005.

**1.2 Clearing (twist).** A column whose index is the pivot row of a reduced column one degree up reduces to zero, so it is skipped.
Chen, Kerber, *Persistent homology computation with a twist*, EuroCG 2011. Also implemented in PHAT: Bauer, Kerber, Reininghaus, Wagner, *PHAT*, J. Symb. Comput. 2017.
Used in `explicit_engine.hpp` (cohomology and homology directions) and `rips_engine.hpp`.

**1.3 Cohomology instead of homology.** Persistent cohomology yields the same barcode. Combined with clearing it is usually much cheaper.
de Silva, Morozov, Vejdemo-Johansson, *Dualities in persistent (co)homology*, Inverse Problems 2011.
This is the default direction for explicit complexes and for all Rips degrees ≥ 2.

**1.4 H0 by union-find.** Degree 0 is a minimum spanning forest computation. Tarjan's union-find with path compression is enough.
Used whenever every 1-cell has boundary ±(v − w) (`Prepared::graph_like`), and always for Rips.

**1.5 Forward cohomology with annotations.** Simplices are inserted in filtration order. Each k-cell carries its coordinates in the basis of the currently alive k-cocycles, and an inverted index records which cells carry each coordinate. A new (k+1)-cell either creates a class, when its boundary has zero coordinates, or kills the youngest coordinate.
Dey, Fan, Wang, *Computing topological persistence for simplicial maps*, SoCG 2014; Boissonnat, Dey, Maria, *The compressed annotation matrix*, Algorithmica 2015 (ESA 2013). It is the algorithm behind GUDHI's persistent cohomology.
*Ours:* we do not compress identical annotations. In `explicit_engine.hpp` (`run_forward`) it is generalised from our H1 Rips engine (1.10). A size-and-level heuristic (`auto_direction`) picks forward or backward, calibrated in `bench/direction.py`.

**1.6 Homology direction with clearing from the top.** Available as `--direction homology`. It never won in our benchmarks.

### Vietoris–Rips specific

**1.7 Implicit Rips coboundaries.** The coboundary matrix is never stored. Simplices are indexed by the combinatorial number system (colexicographic rank), the filtration order is computed on the fly, and reduction uses a heap and lazy column operations. Only the reduction matrix is stored, not the reduced matrix.
Bauer, *Ripser: efficient computation of Vietoris–Rips persistence barcodes*, JACT 2021. The combinatorial number system itself is classical.
*Ours:* the same design, with coefficients generalised to multi-modular lanes (§2).

**1.8 Apparent and emergent pairs.** If τ is the earliest cofacet of σ and σ is the latest facet of τ, then (σ, τ) is a persistence pair and needs no reduction. *Emergent* pairs are recognised when a column is first examined.
Bauer, *Ripser* (as in 1.7); Bauer, Roll, *Gromov hyperbolicity, geodesic defect, and apparent pairs in Vietoris–Rips filtrations*, SoCG 2022.
*Ours:* apparent pairs of non-cone levels are recognised implicitly. In each level the youngest facet is the cell minus its smallest vertex outside the level edge, and the earliest cofacet is found by an early-exit scan in colex order. Nothing is stored.

**1.9 Dominated edges (edge collapse).** An edge e = ab is dominated by a vertex v when its closed neighbourhood lies in v's. Removing dominated edges, or shifting them later in the filtration, preserves flag persistence.
Barmak, Minian, *Strong homotopy types, nerves and collapses*, DCG 2012 (dominated vertices); Boissonnat, Pritam, *Edge collapse and persistence of flag complexes*, SoCG 2020; Glisse, Pritam, *Swap, shift and trim to edge collapse a filtration*, SoCG 2022.
*Ours:* `collapse.hpp` implements the backward shift-and-trim pass with bitsets: the graph at the current time and the final graph are kept as bitsets. Common neighbourhoods and dominator tests are word-parallel. An edge whose dominator is already adjacent to all of its future common neighbours is trimmed with no time-ordered scan. On dragon (1.39M edges) the first pass takes 0.9 s, against 21.7 s for GUDHI's collapser. Our "cone level" condition (§4.1) is exactly domination of the edge in the subgraph present at its insertion time.

**1.10 Active enumeration for H1.** Only triangles with an active face are enumerated, and zero-length classes take a fast path that never stores them.
Kapulkin, Kershaw, *RedZeD: computing persistent homology by reduction to zero differentials*, arXiv:2606.06310 (2026).
*Ours:* re-implemented with multi-modular coefficients in `rips_engine.hpp` (`h1_degree`).

### Explicit complexes

**1.11 Discrete Morse reduction before persistence.** This was Perseus 4/5's coreduction engine.
Forman, *Morse theory for cell complexes*, Adv. Math. 1998; Mischaikow, Nanda, *Morse theory for filtrations and efficient computation of persistent homology*, DCG 2013.
Perseus 6 does **not** currently pre-reduce (see §5). Lower-star versions for images: Robins, Wood, Sheppard, *Theory and algorithms for constructing discrete Morse complexes from grayscale digital images*, IEEE TPAMI 2011.

**1.12 Cubical complexes in doubled ("Khalimsky") coordinates.** Separable min-propagation of top-cell values, with periodic axes. This is standard; see GUDHI's cubical complexes (Dłotko) and Cubical Ripser (Kaji, Sudo, Ahara, arXiv:2005.12692).

**1.13 Implicit cubical persistence (`cubical_engine.hpp`).** The design is Cubical Ripser's:
* cells are never stored, and coboundaries are generated from grid positions;
* H0 is union-find;
* degree d−1 is union-find on the dual graph, because δ^(d−1) is an incidence matrix (their "joint pairs", an instance of Alexander duality);
* the degrees between are cohomology with clearing and apparent pairs (1.8).

Engineering choices of ours:
* births are stored as 32-bit indices into the sorted distinct values;
* only dimensions 1..d−1 are sorted, with a stable 8-bit radix sort;
* the union-finds are indexed by grid position, for locality;
* the working column switches from a sorted vector to a heap when it grows long, after the Ripser heap; without this, smooth 4D inputs were quadratic.

**1.14 Sparse flag complexes (`rips_sparse.hpp`).** When the graph at the threshold is sparse, the radius graph is built directly (a uniform grid for low-dimensional point clouds), and the flag complex is enumerated from sorted adjacency lists. This is what GUDHI's simplex tree does; memory is linear in the number of simplices.

## 2. Coefficients and arithmetic

**2.1 Computation over Q: exact, or multi-modular with random primes.** By default the rational run uses exact int64 arithmetic with overflow checks (bigints in proved mode). When int64 overflows, the run is carried through K random primes from [2^61, 2^62) at once, and repeated with fresh primes if a pivot vanishes modulo some of them but not all. This is the "unlucky primes" method of computer algebra (von zur Gathen, Gerhard, *Modern Computer Algebra*), with the primes drawn at random. With a fixed set of primes the method can be fooled by an input built for them; Perseus 6.0's first release could be. The error bound uses Hadamard's inequality and a rigorous lower bound for the number of primes in the window (Dusart, *Estimates of some functions over primes*, 2010); see `docs/CERTIFICATE.md`. For persistence over several fields in one pass: Boissonnat, Maria, *Computing persistent homology with various coefficient fields in a single pass*, JACT 2019 (ESA 2014).

**2.2 Montgomery multiplication.** Montgomery, *Modular multiplication without trial division*, Math. Comp. 1985. Implemented in `ring.hpp`, with 32-bit lanes for moduli below 2^31 and 64-bit lanes below 2^63.

**2.3 CRT and rational reconstruction.** Wang, *A p-adic algorithm for univariate partial fractions*, SYMSAC 1981; Wang, Guy, Davenport, *P-adic reconstruction of rational numbers*, SIGSAM Bull. 1982. We use the Wang bound √(M/2) on all lanes but two, and accept a value only if those two independent random primes confirm it (`certify.hpp`, `numtheory.hpp`).

**2.4 Big integers.** Schoolbook multiplication and Knuth's Algorithm D for division (TAOCP vol. 2, §4.3.1), in `bigint.hpp`.

**2.5 Primality and factoring.** Miller–Rabin with the fixed bases {2, …, 37}, which are deterministic below 2^64. Trial division up to 2^16, then division by primes already found in other pivots, then Pollard's rho with Brent's cycle detection and batched gcds (Pollard 1975; Brent 1980), then gcds between leftover cofactors.

Numbers that still resist are not factored further. A run modulo such an N either splits it (a pivot shares a factor with N) or proves that all its prime factors give the same barcode. This is the idea of computing modulo a number as if it were prime and learning from the failure, as in Lenstra's elliptic-curve method; here it settles every prime factor at once. Multi-limb Montgomery arithmetic is in `numtheory.hpp`.

## 3. Field dependence: the mathematics we rely on

**3.1 Bad primes and relative torsion.** The F_p barcode differs from the Q barcode exactly when some relative homology group H_q(X_n, X_m; Z) has p-torsion.
Obayashi, Yoshiwaki, *Field choice problem in persistent homology*, DCG 2023. Their algorithm decides field independence for a given complex.

**3.2 Torsion in random complexes.** Context for the experiments, not used by the code:
Kahle, Lutz, Newman, Parsons, *Cohen–Lenstra heuristics for torsion in homology of random complexes*, Exp. Math. 2020; Linial, Peled, *On the phase transition in random simplicial complexes*, Ann. Math. 2016; Aronshtam, Linial, Łuczak, Meshulam, *Collapsibility and vanishing of top homology in random simplicial complexes*, DCG 2013; Kahle, *Topology of random clique complexes*, Discrete Math. 2009.

## 4. What we added, and how it relates to the above

**4.1 Cone levels in place (`rips_engine.hpp`).**

* Levels are the ranks of the longest edge.
* A level e = ab is a *cone level* when the common earlier neighbourhood L_e has an apex v adjacent to all of L_e − v. Equivalently, e is dominated at its insertion time (1.9).
* In the refined order (level, [apex ∉ cell], colex), every pair (β, β ∪ v) of a cone level is an apparent pair (1.8). So these cells are never enumerated as columns.
* *Projection lemma:* rows of non-apex ("base") cells of cone levels are never pivots, so they are deleted from coboundaries.

*Relation to prior work.* This is a lazy, single-pass realisation of Boissonnat–Pritam edge collapse inside a Ripser-style reduction. It needs no rewriting of the graph and no iteration. The freedom to reorder within a level is Glisse–Pritam's "swap". Iterated collapse, using "shift" and "trim", removes strictly more edges.

*Interaction with collapse.* Shifting a dominated edge moves it to a time where, by construction, it is no longer dominated, so shifting removes cone levels. When a collapse shifts many edges but removes few (sphere: 13,177 shifted, 189 removed), the reduction gets *slower* (0.60 s against 0.34 s). Perseus 6 therefore keeps a collapse only if its first pass removes at least 10% of the edges, and abandons it early (after 10% or a third of the first pass) if the removal rate is far below that. Collapse is used only for maxdim ≥ 2; for H1 the active engine (1.10) is cheaper than collapsing.

*Measured against GUDHI edge collapse followed by Ripser* (`bench/vs_collapse.py`, F_2, identical barcodes):

| case | edges before → after collapse | collapse + Ripser | Perseus 6 |
|---|---|---|---|
| o3_1024 H≤3 | 33,189 → 7,179 | 0.20 s | 0.06 s |
| o3_4096 H≤3 | 234,145 → 29,724 | 1.43 s | 0.72 s |
| random16 H≤7 | 942 → 331 | 0.02 s | 0.01 s |
| fractal H≤2 | 122,835 → 27,871 | 6.22 s | 1.21 s |
| sphere H≤2 | 17,857 → 17,668 | 1.84 s | 0.50 s |
| dragon H≤1 | 1,388,645 → 12,034 | 22.2 s | 0.64 s |

Before our own collapser was added, the collapse pipeline beat us on o3_4096 and random16. It no longer does.

**4.2 Certificate of all bad primes (`certify.hpp`, `docs/CERTIFICATE.md`).**

1. Run over Q: exact arithmetic if it suffices, otherwise random 62-bit primes with a stated error bound (2.1).
2. The bad primes are among the prime divisors of the numerators and denominators of the pivot values (Lemma 2). Pivot values are invariants of the filtered matrix: ratios of minors of the blocks of nested pairs. So they are the same for every reduction algorithm, which is what makes a Hadamard-type bound possible.
3. Reconstruct the values (fast mode) and factor them (2.5).
4. Settle each candidate prime by a rerun mod p that starts in the first degree where p divides a pivot, reusing the clearing data of the Q run. Settle each unfactored number by one run modulo it (2.5).

The output is the Q barcode, a complete description of the bad primes (explicitly, or as "every prime factor of N"), and for each one the intervals that differ.

*Relation to prior work.* Boissonnat–Maria compute the barcodes for a prescribed finite set of fields in one pass. Obayashi–Yoshiwaki decide independence for a given complex. We know of no tool that returns the full set of primes with the differences, but the underlying principle is standard computer algebra.

**4.3 Reduction direction (`explicit_engine.hpp`).** Explicit complexes support backward cohomology, forward cohomology (1.5) and homology. A cheap rule picks between the first two: forward when a degree with many top cells sits over a degree whose cells arrive in a few large batches. Engineering, not theory.

**4.4 Fast iterated edge collapse.** This is engineering on Glisse–Pritam's algorithm (1.9): bitset common neighbourhoods and the "already covered" trim test. Together with the keep-only-if-it-removes rule above, this is why the combination beats both edge collapse followed by Ripser and our engine without collapse.

**4.5 Multi-modular lanes throughout.** Every engine is templated on a coefficient ring: F_2 as bytes, F_p, or Q as K Montgomery lanes. The certificate layer is generic over backends. Engineering.

## 5. Not used yet, but planned

* Settling many candidate primes in one run modulo their product (2.5). Classification dominates on torsion-heavy complexes such as Linial–Meshulam near c_2.
* Discrete Morse pre-reduction (1.11) with unit incidences in Z[t^±1], for twisted coefficients (Milnor fibres, Novikov).
* Parallel column reduction: Morozov, Nigmetov, *Towards lockfree persistent homology*, SPAA 2020; giotto-ph (Burella Pérez et al., 2021).
