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
Perseus 6 uses a filtered discrete gradient for cubical grids in dimensions 2–4 (4.6). It does not pre-reduce other complexes (see §5). Lower-star versions for images: Robins, Wood, Sheppard, *Theory and algorithms for constructing discrete Morse complexes from grayscale digital images*, IEEE TPAMI 2011 (ProcessLowerStars). Union-find on critical cells along gradient paths for the extreme degrees: Guillou, Vidal, Tierny, *Discrete Morse Sandwich: fast computation of persistence diagrams for scalar data*, IEEE TVCG 2023.

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

## 5. Not used yet, but planned

* Settling many candidate primes in one run modulo their product (2.5). Classification dominates on torsion-heavy complexes such as Linial–Meshulam near c_2.
* Discrete Morse pre-reduction (1.11) beyond cubical grids, and with unit incidences in Z[t^±1] for twisted coefficients (Milnor fibres, Novikov).
* Parallel column reduction: Morozov, Nigmetov, *Towards lockfree persistent homology*, SPAA 2020; giotto-ph (Burella Pérez et al., 2021).
