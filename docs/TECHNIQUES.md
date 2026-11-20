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
Perseus 6 uses a filtered discrete gradient for cubical grids in dimensions 2–4 (4.6), and an algebraic Morse reduction of unit pairs for chain complexes over Z and over Z[t, t^-1] (5.2, 5.6). It does not pre-reduce simplicial or Rips inputs (see §6). Lower-star versions for images: Robins, Wood, Sheppard, *Theory and algorithms for constructing discrete Morse complexes from grayscale digital images*, IEEE TPAMI 2011 (ProcessLowerStars). Union-find on critical cells along gradient paths for the extreme degrees: Guillou, Vidal, Tierny, *Discrete Morse Sandwich: fast computation of persistence diagrams for scalar data*, IEEE TVCG 2023.

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
* the working column switches from a sorted vector to a bit tree when it grows long (`bitcol.hpp`; below). Without a switch, smooth 4D inputs were quadratic. Until November 2026 the switch went to a binary heap, as in Ripser, which keeps duplicate rows until they surface.

*Bit-tree working columns (`bitcol.hpp`).* PHAT's "bit tree" pivot column is a 64-ary tree of bit words over the rows: a leaf bit per row, and a bit per nonzero word one level up. Bauer, Kerber, Reininghaus, Wagner, *Phat — Persistent Homology Algorithms Toolbox*, J. Symbolic Comput. 78 (2017).
* Over F_2, adding an entry flips a bit, so a duplicate cancels at once. The pivot (the smallest row) is found by descending log₆₄ n words.
* We add a dense coefficient array, so the same structure serves F_p, multi-modular Q and exact int64 rationals: adding at a set bit adds coefficients and clears the bit on zero.
* The tree is allocated once per degree, on the first long column. Storing a reduced column drains the tree in row order, touching only its own entries, so no sort is needed.
* Columns move into the tree past 512 entries (`PERSEUS_LONG_COLUMN` overrides; 0 sends every column through it, which the tests use). Thresholds of 128–512 measured alike on white noise 40⁴; 2048 and 8192 were slower.

**1.14 Sparse flag complexes (`rips_sparse.hpp`).** When the graph at the threshold is sparse, the radius graph is built directly (a uniform grid for low-dimensional point clouds), and the flag complex is enumerated from sorted adjacency lists. This is what GUDHI's simplex tree does; memory is linear in the number of simplices.

**1.15 2D images by two union-finds and Alexander duality.** For a top-cell filtration of a planar grid, H0 is the connected components of the pixels present and H1 is, by Alexander duality, the components of the complement swept downwards. Both are union-finds. GUDHI's image algorithm does this: a discrete gradient inside each pixel leaves critical edges, which are sorted. A dual pass over them in decreasing order gives degree 1, and the edges it leaves give degree 0 in a primal pass.
Glisse, `Persistence_on_rectangle_from_top_cells.h` (GUDHI, 2023, MIT licence), used by GUDHI's `CubicalPersistence` for 2D inputs. Kaji, Sudo, Ahara (Cubical Ripser, 1.12) use the same duality for the top degree ("joint pairs"). Le Breton, Szustakowski, Piraud, *Fast cubical persistent homology on 2D and 3D images via union-find, pruning, and lookup tables* (Flash Cubical, arXiv:2606.04801, 2026) combine union-find, duality and lookup tables for V-filtrations over F_2. That is a different complex from the T-construction, so we do not benchmark against it.
*Ours* (4.7): we classify pixels, not cells. The local counts come from one byte per pixel and two 256-entry tables. Only births and merge candidates are sorted. Degree 1 runs first, and the Euler characteristic of the pixel's local picture tells the degree 0 pass which candidates to skip. That skip is our analogue of GUDHI's handing the leftover edges from the dual pass to the primal one. The certificate is total unimodularity. Absent pixels give essential classes. We read GUDHI's source for this idea and then wrote the code independently.

## 2. Coefficients and arithmetic

**2.1 Computation over Q: exact by default, multi-modular with random primes on request.** The rational run uses exact int64 arithmetic with overflow checks, and where that overflows, exact rationals of any size (`QHyb`, 2.4). With `--mode fast` the run is instead carried through K random primes from [2^61, 2^62) at once, and repeated with fresh primes if a pivot vanishes modulo some of them but not all. This is the "unlucky primes" method of computer algebra (von zur Gathen, Gerhard, *Modern Computer Algebra*), with the primes drawn at random. With a fixed set of primes the method can be fooled by an input built for them; Perseus 6.0's first release could be. The error bound uses Hadamard's inequality and a rigorous lower bound for the number of primes in the window (Dusart, *Estimates of some functions over primes*, 2010); see `docs/CERTIFICATE.md`. For persistence over several fields in one pass: Boissonnat, Maria, *Computing persistent homology with various coefficient fields in a single pass*, JACT 2019 (ESA 2014).

**2.2 Montgomery multiplication.** Montgomery, *Modular multiplication without trial division*, Math. Comp. 1985. Implemented in `ring.hpp`, with 32-bit lanes for moduli below 2^31 and 64-bit lanes below 2^63.

**2.3 CRT and rational reconstruction.** Wang, *A p-adic algorithm for univariate partial fractions*, SYMSAC 1981; Wang, Guy, Davenport, *P-adic reconstruction of rational numbers*, SIGSAM Bull. 1982. We use the Wang bound √(M/2) on all lanes but two, and accept a value only if those two independent random primes confirm it (`certify.hpp`, `numtheory.hpp`).

**2.4 Big integers and exact rationals.** Schoolbook multiplication and Knuth's Algorithm D for division (TAOCP vol. 2, §4.3.1), in `bigint.hpp`. Limbs are kept inline up to 256 bits, so the numbers a reduction meets do not touch the heap.
* *Gcds.* Euclid's remainders run until both numbers fit in 128 bits, then binary gcd takes over (Stein 1967; Knuth §4.5.2), with one-word shortcuts.
* *Exact rationals (`QHyb` in `ring.hpp`).* Values stay as int64 pairs while they fit, and as big rationals beyond. Sums and products use Henrici's formulas, which cancel common factors before multiplying, so only gcds of the small cross terms are needed (Henrici, *A subroutine for computations with rational numbers*, JACM 1956; Knuth §4.5.1).
* *Effect.* On the n = 90 Linial–Meshulam instance the exact Q run went from 202 s to 24 s.

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

**4.6 A filtered discrete gradient for cubical grids (`cubical_morse.hpp`).** The default engine for top-cell filtrations in dimensions 2–4 without periodic axes. It is Mischaikow–Nanda's filtered Morse reduction (1.11), specialised to grids. Neither Cubical Ripser nor the classic engine (1.13) uses it.

* *Blocks.* Top cubes are totally ordered by (birth, position). A cell enters with its first top coface τ. The cells entering with τ form L(τ) = cl τ − A(τ), where A(τ) is the union of the faces τ shares with the earlier cubes among its 3^d − 1 neighbours. A(τ) depends only on which neighbours are earlier. By excision, rank H_*(cl τ, A(τ)) counts the bars born or dying at τ. Summed over τ, this is 2·#(finite bars) + #(essential bars), the least number of critical cells any filtered gradient can have.
* *Patterns.* Dually, L(τ) is the star of a vertex in the dual grid. Its link is a subcomplex K of the boundary of the d-dimensional cross-polytope, and every such K occurs. The block has a perfect matching exactly when K, with an empty cell added, has a perfect discrete Morse function. H̃(K) is torsion-free, by Alexander duality in S³. So on a perfect block the Morse differential inside the block vanishes over Z.
  * In 2D (14 patterns up to symmetry) and in 3D (522 patterns, 15,935 subcomplexes), every pattern is perfect, checked exhaustively. This agrees with Robins–Wood–Sheppard for lower stars. In 3D the deterministic rule alone is perfect on all 15,935.
  * In 4D we do not know. Every pattern we have met is perfect: the exhaustive enumeration through 495,048 orbits, and about 460,000 distinct patterns from smooth, white-noise and quantised 4D images. No dunce-hat-type obstruction exists: none of the 32.5 million triangle sets without a free edge in the boundary of the cross-polytope is acyclic.
* *Matching per pattern.* A deterministic greedy removal from the top. It is accepted when its count of critical cells equals |χ(L)|, which is a lower bound. Otherwise the count is compared with the F_2 rank, and when they differ the engine tries a greedy that also removes from the bottom, then an exact search over removal orders. Patterns are cached by their 3^d-bit mask. A pattern that is not perfect would only cost extra critical cells, never correctness.
* *Order inside a block.* A removal from the top lists every matched face just before its coface, and every other incidence coface first. So the reverse of the removal sequence orders the block's modified Hasse diagram. The key (birth, block, order) then strictly increases along every coboundary V-path step and strictly decreases along every boundary step. That makes it a filtration of the Morse complex, and a processing order for path sums (a cell's key can even be taken from its partner, so a path step reads only the partner's word).
* *Reductions.* H0 is union-find on critical vertices: the Morse boundary of a critical edge joins the ends of the descending paths from its two vertices, memoised so every vertex is walked once. Degree d−1 is the dual union-find on critical top cells. In both, the Morse complex is again a graph incidence matrix, so all pivots are ±1. This is the strategy of Guillou–Vidal–Tierny's Discrete Morse Sandwich. The degrees between use the explicit Morse coboundary, whose integer path sums are computed once and checked for int64 overflow. They are reduced by cohomology with clearing. Apparent pairs are decided on the integer matrix, using only ±1 entries, so the decision holds in every ring.
* *Coefficients.* Path sums are integers, kept within ±2^62 by checked arithmetic, so every product with an incidence ±1 is exact. A larger sum stops the construction: the default then falls back to the classic engine, and `--engine morse` reports an error. Entries that vanish mod p are dropped. Q mode is exact, and the bad-prime list is unchanged (`docs/CERTIFICATE.md`). Zero-length bars are counted rather than paired: n_k(v) = births_k(v) + deaths_{k−1}(v).
* *Engineering.* Birth index, mate, block offset and order are packed into one 64-bit word per position. Critical cells keep their index there instead. Births are radix-sorted, by the high halves of order-preserving keys, with short runs fixed up afterwards. Earlier neighbours are found as an 80-bit word from a compact array of top-cube births, and A(τ) is assembled from byte tables. Blocks are processed in tiles. Path ends and Morse columns are computed in position order, for locality.

**4.7 The planar engine for 2D images (`cubical_2d.hpp`).** This is the default for top-cell filtrations of 2D grids without periodic axes; `--engine planar` selects it explicitly. Pixels are totally ordered by (value, index), as in the other engines. X_s is the union of the closed squares of the pixels up to s; non-finite values (NaN, ±inf) are absent.

* *H0.* Two closed squares meet exactly when their pixels are 8-adjacent. So H0(X_s) is free on the 8-connected components of the pixels present.
* *H1 by duality.* X_s is a compact polyhedron in S² = R² ∪ {∞}. Alexander duality gives H̃_1(X_s; Z) ≅ H̃^0(S² − X_s; Z), naturally in s. The complement consists of the open squares of the pixels not yet present, the absent pixels and the outside, glued along open edges. Two of these regions touch when they share an edge. Two diagonal regions whose other two neighbours are present meet only at a corner, and that corner lies in X_s. So the complement components are the 4-connected components of {later, absent, outside} pixels, with one node for the outside. As s decreases, pixels join the complement. The H1 barcode of X is therefore the H0 barcode of this decreasing sweep, read backwards. Components that exist before any pixel joins (the outside, and each 4-connected cluster of absent pixels) are born at +∞. When two of them meet at p, the class [v(p), ∞) is essential.
* *Local counts.* The earlier neighbours of p are one byte E: bit j is set for ring neighbour j (E, NE, N, NW, W, SW, S, SE). Ties count as earlier for the neighbours of smaller index (W, SW, S, SE). The complement mask is ~E, because absent and outside cells are never earlier. Two tables of 256 entries give:
  * c0(p), the number of local groups of earlier neighbours. Cells are connected when they are consecutive on the ring, or when they are two edge neighbours at a right angle, since those share a corner of p.
  * c1(p), the number of runs of complement cells, consecutive on the ring, that contain an edge neighbour.
  * a representative cell for each group, and for a single group the cell nearest the arc's circular mean. The pointer then approximates steepest descent, which keeps chains short on smooth data; the first cell in ring order would follow level sets.
* *Who must be processed.* If c0 = 0, p is an H0 birth. If c0 = 1, p joins one component and nothing changes, so it gets a union-find pointer at once and is never processed. If c1 = 0, p fills a hole: it is an H1 death and a root of the downward sweep, and it needs no sorting either. Otherwise p is a *candidate*, and c0 = c1 = c ≥ 2: on the circle ∂p the earlier part and the complement part are alternating arcs. Only births and candidates are sorted.
* *Euler characteristic.* Adding the closed square p along A = p ∩ X_{<p}, which is a union of c arcs of ∂p, changes χ by 1 − c. Neither degree loses a class unless A is all of ∂p. So at a candidate, (H0 merges) + (H1 births) = c − 1. H1 runs first, over the candidates in decreasing order. When the c complement groups lie in r distinct components, r − 1 classes are born at p. If r = c, every event at p was an H1 birth, and the H0 pass treats p as a pointer. On white noise, smoothed noise and the photograph this skips 55–62% of the candidates.
* *Union-finds.* One int32 array over the padded grid (a frame of outside pixels) serves H1 and is then reset for H0. A root is a negative code:
  * H0: −1 − (rank of the birth), so the elder rule compares ranks;
  * H1: −4 − f for the f-th hole fill, whose value is kept in a compact array; −2 for an absent cluster; −3 for the outside.

  Finds use path halving. On large inputs with sparse candidates, a pre-pass in index order compresses every pointer chain that starts next to a candidate. Chains are spatially local, so this is cache friendly, and afterwards finds in sorted order cost one step.
* *Ties and zero-length bars.* Positive-length bars do not depend on how ties are broken, so they agree with the other engines. Zero-length bars, when requested, are counted per value. Every k-cell of the doubled grid is a birth in degree k or a death in degree k − 1. So z_k(v) = n_k(v) − b_k(v) − z_{k−1}(v) − d_{k−1}(v), where n counts cells of value v, b counts positive-length births and d positive-length deaths.
* *Certificate.* Both union-finds reduce graph incidence matrices: the 1-skeleton of pixels, and the dual graph with the outside node. These are totally unimodular, so every pivot is ±1 (`docs/CERTIFICATE.md` §5). The certificate reports mode `unimodular` and does no arithmetic.
* *Engineering.*
  * The neighbour bytes of a row are computed four pixels at a time with vector compares (GCC/Clang vector extensions; SSE2 or NEON width, so no special build flags are needed).
  * The bytes and pointers are written in one streaming pass. Rows with few special pixels branch on them; dense rows append without branches.
  * Keys are order-preserving 64-bit images of the values, radix-sorted by their high halves; runs of equal high halves are fixed up afterwards and skipped when all their keys are equal.
  * The random-access arrays are allocated on transparent huge pages where Linux offers them.
  * The processing loops prefetch two levels ahead: the pixel's neighbourhood, then the pointer targets.
  * The Python module passes the NumPy array to the engine without copying it.

## 5. Chain complexes over Z, simplicial sets, magnitude homology, Laurent coefficients, finite groups

**5.1 Filtered based chain complexes over Z (`chain.hpp`).** Cells with integer boundary coefficients and births. Validation checks that faces are earlier cells of one dimension less, born no later, and that d∘d = 0 exactly (128-bit accumulation, bigints when that overflows).

**5.2 Algebraic Morse reduction (`chain.hpp`, `algebraic_morse`).** A pair σ < τ with [τ : σ] = ±1 and equal births is cancelled by the change of basis ρ' = ρ − [ρ : σ][τ : σ] τ for the other cofaces ρ of σ, and σ' = dτ. Every term is born no later than the cell it replaces, so the complex splits, as a filtered complex over Z, into the residual complex and acyclic pieces τ → dτ. Every barcode over every field, the bad primes, H_*(K_s; Z) for every s and every map H_*(K_s; Z) → H_*(K_t; Z) are therefore those of the residual complex; a cancelled pair is a zero-length bar. The reduction runs in two phases. First the cancellations without fill-in: coreductions (τ has one live face) and collapses (σ has one live coface), from a work list. Then Gaussian elimination of the remaining unit pairs, least Markowitz cost first, from a lazily updated heap; each elimination may unlock the first phase again. Coefficients are int64 with overflow checks; a pair whose update would overflow is left alone.
Forman, *Morse theory for cell complexes*, Adv. Math. 1998; Sköldberg, *Morse theory from an algebraic viewpoint*, Trans. AMS 2006; Jöllenbeck, Welker, *Minimal resolutions via algebraic discrete Morse theory*, Mem. AMS 2009; Kozlov, *Combinatorial Algebraic Topology*, 2008. Coreductions: Mrozek, Batko, *Coreduction homology algorithm*, DCG 2009. For filtrations: Mischaikow, Nanda (1.11). Nothing here is new in principle.

**5.3 Integral outputs (`zlin.hpp`, `chain.hpp`).** Smith normal form by the classical algorithm, with a global least-modulus pivot and rounded quotients at every round (pivoting on the first nonzero entry blows up even on small matrices), in int64 with overflow checks and repeated with bigints when a check fails. Kernels by column echelon form with extended-gcd column operations and the transformation matrix. The image of H_k(K_s; Z) → H_k(K_t; Z) is the lattice quotient (Z_k(K_s) + B_k(K_t)) / B_k(K_t): an echelon basis of the numerator, the coordinates of the generators of B_k(K_t) in it, and the Smith form of that coordinate matrix. These are dense, so they run on the residual complex of 5.2. The F_p barcode of each bad prime is a rerun over F_p.
Smith, *On systems of linear indivisible equations and congruences*, Phil. Trans. 1861. Homology over Z by reduction followed by Smith forms is the approach of CHomP: Kaczynski, Mischaikow, Mrozek, *Computational Homology*, 2004. Polynomial-time Smith forms (Kannan, Bachem, SIAM J. Comput. 1979) are not used; entries stay small on the residual complexes we meet.

**5.4 Simplicial sets (`simplicial_set.hpp`).** Only nondegenerate simplices are stored. A face is a nondegenerate simplex y with a degeneracy word in its canonical form s_j1 … s_jm y, j1 > … > jm (every composite of degeneracies has exactly one such form). Faces of degenerate simplices follow from the simplicial identities, which are checked on every simplex. The chain complex is the normalised one: degenerate faces are dropped, repeated faces add up, so coefficients 0 and ±2 occur. It is quasi-isomorphic to the unnormalised complex (Eilenberg–Mac Lane), filtration by filtration. Constructions: the nerve of a finite category (a chain of composable non-identity morphisms is a nondegenerate simplex, inner faces compose, identities become degeneracies), the bar construction BG, quotients X/A, cones, suspensions, minimal spheres. In a weighted nerve the birth of a chain is the largest weight of a composite of consecutive morphisms in it. Every such composite of a face is one of the simplex, so faces are never born later.
May, *Simplicial Objects in Algebraic Topology*, 1967; Goerss, Jardine, *Simplicial Homotopy Theory*, 1999. Homology of finite simplicial sets is in Sage (`sage.topology.simplicial_set`, J. Palmieri) and in Kenzo (Sergeraert et al.). We know of no other program that computes persistence of filtered simplicial sets, or certifies its field dependence.

**5.5 Magnitude homology of graphs (`magnitude.hpp`).** MC_{k,l}(G) is free on tuples (x_0, …, x_k) with x_i ≠ x_{i+1} and total length l; d = Σ_{i=1}^{k−1} (−1)^i d_i, where d_i deletes x_i if it lies on a geodesic from x_{i−1} to x_{i+1}. The differential keeps the endpoints and the length, so the complex splits into pieces MC_{*,l}(s, t), reduced in parallel. Each piece is enumerated by depth-first search with ascending neighbours, which lists the tuples of each size in lexicographic order, so faces are found by binary search without sorting. The piece is stored once in compressed rows (int32 faces, int8 coefficients ±1: two different deletions never give the same tuple) with the compressed coboundary. The fill-in-free phase of 5.2 runs in place on these arrays, with live face and coface counts; it removes 90–99.9% of the cells on the torsion families below. The survivors go through the general reduction and the Smith forms. Peak memory is about 80 bytes per tuple.
Hepworth, Willerton, *Categorifying the magnitude of a graph*, Homology Homotopy Appl. 2017; Leinster, *The magnitude of a graph*, Math. Proc. Camb. Phil. Soc. 2019. Torsion: Kaneta, Yoshinaga, *Magnitude homology of metric spaces and order complexes*, Bull. LMS 2021 (for the Hasse diagram of P with a bottom 0̂ and a top 1̂ added, MH_{k, r+2}(0̂, 1̂) ≅ H̃_{k−2}(Δ(P)), so torsion in the order complex appears in magnitude homology); Sazdanovic, Summers, *Torsion in the magnitude homology of graphs*, J. Homotopy Relat. Struct. 2021. The splitting by endpoints is standard; gramag (T. Chaplin, Rust) uses it, with ranks over F_2 from the lophat reduction.

*Measured* (`bench/mag_families.cpp`, pair (0̂, 1̂) of the Kaneta–Yoshinaga graph, one thread): see the table in the README. The size limit is memory: each step of 2 in l multiplies the number of tuples by about 30–75.

**5.6 Laurent coefficients (`laurent.hpp`).** Chain complexes over Z[t, t^-1]. This ring is not a principal ideal domain, so its homology has no barcode and no Smith form. Perseus computes what the principal ideal domains around it see:

| ring | what is computed |
|---|---|
| Q(t) | the barcode (equal to the barcode over the Novikov field Q((t))) |
| Q[t]/(f), f irreducible, f ≠ t | which f give another barcode (bad irreducibles), and that barcode |
| F_p(t) | which p give another barcode (bad primes), and that barcode |
| Q[t, t^-1] | the homology of the whole complex: free rank and invariant factors |
| F_p[t, t^-1] | the primes p where it is not the reduction of the rational answer, and the homology there |

Arithmetic is fraction-free on integer polynomials: pseudo-division, gcds by the primitive remainder sequence, exact division. A nonzero integer is a unit in every ring above except F_p, and t is a unit in all of them, so rows and columns are scaled by integers and powers of t freely. The integers that F_p could see are recorded. Smith forms over Q[t] (a Euclidean ring: least-degree pivot, pseudo-division). The barcodes come from a column reduction over the field, with clearing, normalising each column by the gcd of its entries and tracking the column's scale, so the pivot values (ratios of minors, as in `CERTIFICATE.md` §1) are known exactly. The candidates are a gcd-free basis of the primitive parts of their numerators and denominators, with linear factors split off by the rational root test, and the prime factors of their contents. A candidate g is settled by one reduction over Q[t]/(g) as if it were a field. A pivot that is not invertible modulo g splits g by a gcd, and the parts are settled in turn; otherwise all irreducible factors of g give the same barcode. This is the D5 principle of dynamic evaluation. Before all this, pairs with unit incidence ±t^k are cancelled (5.2 over Z[t, t^-1]; only pairs with equal births for the barcodes, all of them for the homology of the whole complex).
Pseudo-remainder sequences: Collins, *Subresultants and reduced polynomial remainder sequences*, JACM 1967; Brown, Traub, JACM 1971; Knuth, TAOCP vol. 2, §4.6.1. Gcd-free (coprime) bases: Bach, Shallit, *Algorithmic Number Theory*, 1996; Bernstein, *Factoring into coprimes in essentially linear time*, J. Algorithms 2005. Dynamic evaluation: Della Dora, Dicrescenzo, Duval, *About a new method for computing in algebraic number fields*, EUROCAL 1985. Novikov homology: Novikov 1981; Farber, *Topology of Closed One-Forms*, AMS 2004. Alexander modules by Fox calculus (used in the tests): Fox, *Free differential calculus I*, Ann. Math. 1953.
*Not computed:* the Z[t, t^-1]-module itself (it is not determined by the fibres above), the residue fields F_p[t]/(h), torsion over the Novikov ring Z((t)) (Farber's torsion numbers), and images H(K_s) → H(K_t) over Q[t, t^-1].

**5.7 Finite solvable groups (`pc_group.hpp`).** Input: the multiplication table of G and a polycyclic generating sequence f_1, …, f_k with relative orders r_1, …, r_k (GAP's `Pcgs`, written by `tools/export_pcgroup.g`). Perseus checks that the table is a group (Latin square and associativity; above order 512 associativity is sampled on 10^6 triples), that the collected words f_1^a_1 ⋯ f_k^a_k (0 ≤ a_i < r_i) are distinct, and that f_i^r_i and f_i^-1 f_j f_i (j > i) lie in G_{i+1} = ⟨f_{i+1}, …, f_k⟩. The collected words are then the normal forms of a complete rewriting system (collection) whose tips, the minimal non-normal words, are f_i^r_i and f_j f_i for j > i.

*The resolution.* The rewriting system induces an algebraic Morse matching on the normalised bar resolution of Z over ZG: cut a bar entry where a tip ends inside it, merge two neighbouring entries whose product is still collected. Its critical cells are Anick's chains, and its Morse complex is Anick's free resolution (Brown 1989; Anick 1986; Sköldberg 2006; Jöllenbeck–Welker 2009). For these tips the n-chains are the non-increasing letter sequences, so there are C(n+k−1, k−1) of them in degree n. Following the gradient flow on the bar construction itself costs about 9× per degree on groups of order 64, whatever the coefficients. Perseus instead builds the resolution as right ZG-modules P_n = Z[chains_n] ⊗ ZG by Anick's recursion with the contracting homotopy
  d_{n+1}(c t) = c ⊗ t − h(d_n(c) · t),   h(c ⊗ g) = (c s) ⊗ u − h(d(c s) · u − c ⊗ g) when g = s u extends c to a chain (else 0),
memoised on the monomials c ⊗ g. That is at most |G| · C(n+k−1, k−1) states in degree n: polynomial in the degree. h is evaluated iteratively, with a guard that reports non-termination (it cannot happen for a pcgs).

*Integral homology.* For n ≥ 1, |G| annihilates H_n(G; Z) (the transfer), so only the primes q dividing |G| occur. At such a prime, with q^v exactly dividing |G|, every nonzero elementary divisor of the augmented differentials has q-valuation at most v. The resolution is therefore built over Z/M with M = ∏ q^(v+1), and the q-part is read off over the local ring Z/q^(v+1). This is exact: no rational arithmetic, no candidate primes, no certificate needed. Over the local ring, Gaussian elimination on a unit entry is a chain homotopy equivalence (the Gaussian elimination lemma, as in 5.2), so unit pivots are cancelled across the whole complex first, least Markowitz cost first. A cancelled pair deletes a row of the next differential and a column of the previous one. A Smith form by least-valuation pivots then runs on the residual, which has no unit entries left. The ranks found at the different primes must agree; that is checked.

*Mod-p cohomology rings.* The resolution again, over F_p. Cocycles are found as the orthogonal complement of the span of the augmented differentials, and a basis of H^n is a complement of the coboundaries in them (sparse echelon forms). Products are Yoneda composites: a cocycle α of degree a lifts through h to a chain map α~ : P_{a+*} → P_* (α~_0 = α · e, α~_j = h(α~_{j−1} ∘ d)), and (α · β)(c) = β(ε α~_b(c)) for c of degree a + b. This is the cup product up to the sign (−1)^(ab), which does not change which classes are decomposable. All products of basis classes with a ≤ b and a + b ≤ maxdim are computed, and the minimal generators in each degree are the complement of the span of the products landing there. Relations are not extracted.

What is new here is the combination, not the parts: Anick resolutions of pc presentations are classical (HAP builds resolutions of nilpotent groups by related perturbation methods; Ellis, *Computing group resolutions*, JSC 2004), and so are minimal resolutions over F_p G for mod-p rings (Carlson; Green and King, *The computation of the cohomology rings of all groups of order 128*, J. Algebra 2011). Perseus contributes the polynomial memoised recursion over Z/q^(v+1) with cross-degree unit cancellation, which gives exact integral homology directly, and an implementation fast enough to compete; see the benchmarks in the README. Anick, *On the homology of associative algebras*, Trans. AMS 1986. Brown, *The geometry of rewriting systems: a proof of the Anick–Groves–Squier theorem*, MSRI 1992. Sköldberg, *Morse theory from an algebraic viewpoint*, Trans. AMS 2006. Jöllenbeck, Welker, *Minimal resolutions via algebraic discrete Morse theory*, Mem. AMS 2009. Collection and pc presentations: Sims, *Computation with Finitely Presented Groups*, 1994. Persistent homology of groups along a series (not implemented): Ellis, King, *Persistent homology of groups*, arXiv:1006.2237.
*Not computed:* relations of the mod-p ring, integral cup products, non-solvable groups (use a Sylow subgroup and stable elements yourself), and maps induced by homomorphisms.

## 6. Not used yet, but planned

* Classifying candidate primes cheaply (2.5): restart each prime from the first column whose pivot it divides, from checkpoints of the Q reduction. Batching by a composite modulus helps only up to the first split, because every candidate divides some pivot. Classification dominates on torsion-heavy complexes such as Linial–Meshulam near c_2. `docs/ROADMAP.md`, step 2.
* Discrete Morse pre-reduction (1.11) for simplicial and Rips inputs; finite twisted coefficients Z[t]/(t^n − 1) for Milnor fibres.
* Parallel column reduction: Morozov, Nigmetov, *Towards lockfree persistent homology*, SPAA 2020; giotto-ph (Burella Pérez et al., 2021).
