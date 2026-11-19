# The Q certificate

Perseus 6 reports the barcode over Q and the set of *bad primes*: the p for which the barcode over F_p
differs. This note says exactly what is guaranteed and why.

## 1. Candidates

Let D be the filtered boundary matrix and run any column reduction R = D V over Q (V upper unitriangular,
R reduced).

**Lemma 2.** If a prime p divides no numerator and no denominator of a pivot value, then the same reduction
taken mod p is a valid F_p reduction with the same pivots, hence the same barcode.

*Proof.* Every multiplier has a pivot value as denominator, so every intermediate entry is p-integral. An
addition whose multiplier vanishes mod p is the identity mod p. Each pivot is nonzero mod p and the entries
below it are zero, so every column keeps its lowest entry and the pairing is unchanged. ∎

So the bad primes lie among the prime factors of the pivot values, and each candidate is settled by a rerun
mod p. That rerun starts in the first degree where p divides a pivot, because lower degrees and their
clearing data are unchanged.

**Pivot values are invariants.** Let (i, j) be a pair. Let J be the paired columns c ≤ j with low(c) ≥ i,
and I their lows: the pairs nested in (i, j). Then R[I, J] is triangular, and Cauchy–Binet with the
unitriangular V gives det D[I, J] = det R[I, J] = ±∏_{c∈J} π(c). Removing (i, j) gives

    π(i, j) = ± det D[I, J] / det D[I \ i, J \ j].

The pivot values therefore depend only on D and the filtration order, not on the reduction. Lower-left
blocks map to lower-left blocks under anti-transposition, so cohomology gives the same values. The test
suite checks this exactly, in exact arithmetic, across all three explicit reduction directions
(`exact_oracle.py`: identical values on 150 random torsion-rich complexes). Two engines are annotation
algorithms rather than column reductions: forward cohomology, and the Rips H1 engine. For those the
identity is verified by these tests but not yet written up as a proof.

**The cubical Morse engine reduces a different matrix.** It reduces the Morse complex M of a filtered
discrete gradient (`cubical_morse.hpp`), not D. Its entries are integer sums over gradient paths. They need
not be ±1: in 4D, entries ±2 occur. M is filtered chain-homotopy equivalent to the cubical complex over Z,
because every matched pair lies inside one filtration step. So over every field F_p and over Q the two
barcodes are equal, and so the bad primes are the same. Lemma 2, applied to M, puts them among the prime
factors of M's pivot values. Those values are invariants of M, not of D. The *candidate* list can therefore
differ between the engines; for example, a pivot of 2 adds the candidate 2. Each candidate is settled by a
rerun, so the list of *bad* primes does not change. Fast mode bounds the minors of M by Hadamard, using
the norms of M's columns (the largest m of them, for pairs bounded by m).

## 2. Three ways to run the Q reduction

**proved** (`--mode proved`). Exact rational arithmetic. Each degree is tried with int64 numerators and
denominators, with every operation overflow-checked. A degree that overflows is rerun with arbitrary
precision. The barcode and the pivot values are exact, so Lemma 2 applies as stated.

**fast** (`--mode fast`). The reduction is carried modulo K primes drawn uniformly at random from
[2^61, 2^62). The primes are distinct. They come from ChaCha20 keyed by a fresh 256-bit seed from the
operating system's random source (`std::random_device`), and the seed is reported.

*Seeds.* The bound below is a probability over the seed, for a seed drawn after the input is fixed. It says
nothing about a seed chosen beforehand: for any given seed one can build an input that fools it. So
`--seed` reproduces a run (same primes, same answer) but does not certify it. Such runs report
`"barcode_certified": false`, `"seed_supplied": true`, `"complete": false`, and bad primes "conditional on
the supplied seed". The seed has 256 bits, so the seed space itself contributes at most 2^-256, and the
analysis treats ChaCha20's output as uniform, the standard assumption for a cryptographic generator.
Perseus 6.0's first build took 64-bit seeds into a non-cryptographic generator and certified runs with a
supplied seed; with a fixed seed an input can be built to defeat it.

*Faithfulness.* Each lane computes the Q reduction mod its prime, up to the first Q-pivot that the prime
divides: until then it drops exactly the entries that vanish mod p, and those additions are no-ops mod p.
If such a pivot vanishes in some lanes but not all, the run is flagged and repeated with fresh primes. It
goes wrong undetectably only if **every** lane prime divides a nonzero integer that the input fixes:

* the numerator of a Q-pivot, which makes the column look zero or its pivot move;
* numerator ∓ denominator of a non-unit Q-pivot, which makes it look like ±1, so it is never factored;
* a ≠ v·b for a wrong reconstruction a/b of a pivot v (see below).

*Bound.* Every minor of the matrix reduced in degree k has log2 |det| ≤ h_k, by Hadamard over rows or over
columns, whichever is smaller (engines record h_k). A positive integer below 2^(h+2) has at most
m = ⌊(h+2)/61⌋ prime factors in the window, and the window contains N ≥ 5.39·10^16 primes. This is a
rigorous lower bound from Dusart's estimates for π(x). Therefore

    Pr[wrong barcode or wrong candidate set] ≤ Σ_k 3 P_k (m_k / N)^K  +  Σ_reconstructed (m'/N)^2,

with P_k at most the number of pairs in degree k, the probability being over the fresh seed. The bound
holds for every input; there is no assumption about where the input came from. Perseus adds lanes until it
is below 2^-64 (`--target`) and reports it.

*Reconstruction.* A non-unit pivot is reconstructed by Wang's method from all lanes but two. It is accepted
only if it agrees with those two independent lanes. Degrees holding non-unit pivots are rerun with more
fresh lanes (4, 8, 16, 32) as needed, and each rerun must reproduce the pairing of the first run.

**auto** (default). First try exact int64 arithmetic. When it does not overflow, which is the common case,
the result is a proof and costs about as much as the multi-modular run. When it does overflow
(torsion-heavy inputs), the fast mode runs and the output says so.

## 3. Classifying candidates

Candidates below 2^63 are settled by reruns over F_p, in parallel. Every candidate is rerun from the first
degree k0 in which it divides a pivot, reusing the Q run's clearing data below k0. This is valid by Lemma 2:
the prime divides no pivot in degrees below k0, so the F_p run agrees with the Q run there. It is essential
that k0 be the *first* such degree. A number that turns up in several degrees (say 3N in degree 1 and N in
degree 2) keeps the earliest. A rerun from a later degree would miss every change below it. The first build
of Perseus 6.0 could keep the later degree for a composite cofactor and miss its primes entirely when the
later degree showed no change.

Pivot values are factored by trial division, then by dividing out primes already found from other pivots,
then by Pollard–Brent rho with batched gcds, then by gcds between the leftover cofactors. Large primes and
unsplit composites are **not** factored further. Each such number N is classified by a single run modulo N
itself, using multi-limb Montgomery arithmetic:

**Lemma 3.** Run the reduction over Z/N, for N odd and possibly composite. Suppose every pivot is either a
unit mod N or, if it shares a factor with N, the run stops and reports gcd(pivot, N), splitting N. If the
run finishes without a split, then for every prime p dividing N it is a valid F_p run with the same pivots.
Hence all prime factors of N have the same barcode.

*Proof.* As for Lemma 2. A pivot that is a unit mod N is nonzero mod every p | N. A dropped entry is zero
mod N, hence mod p. An addition with a multiplier that vanishes mod p is the identity mod p. ∎

So the certificate can say "every prime factor of N is bad, with these interval changes" without knowing
the factors (`"every_prime_factor_of": true` in the JSON). It is still complete. A split refines N, and
the parts are classified in turn. Lemma 3 applies whether or not N is prime, so primes above 2^64, which
Miller–Rabin only identifies as probable primes (`"probable_prime": true`), need no primality proof.

## 4. What the output says

* `mode`: `proved`, `fast`, or `unimodular` (the planar engine, §5: no arithmetic needed).
* `guarantee.barcode_certified`: false only if every set of trial primes was unlucky, or two independent
  runs disagreed. Then nothing is claimed.
* `guarantee.log2_error_bound` (fast mode): the bound of §2.
* `bad_primes_status`: `complete` (every candidate classified), `candidates only` (`--no-classify`), or
  `incomplete`.
* `complete`: true when the barcode is certified and the bad primes are completely described.

## 5. Planar inputs: a proof by total unimodularity

The planar engine (`cubical_2d.hpp`, `docs/TECHNIQUES.md` §4.7) handles top-cell filtrations of 2D grids
without periodic axes. Its certificate needs no arithmetic. The output says mode `unimodular`, with arithmetic
"none needed: totally unimodular boundary matrices". The argument has two steps.

**The engine reduces graph incidence matrices.**
* Degree 0 is the persistence of the filtered graph G0: its vertices are the present pixels and its edges join
  8-adjacent pixels. An edge enters with its later end. The closed squares of a set of pixels are connected
  exactly when the pixels are 8-connected, so H0 of the filtered complex X_s equals H0 of G0 at every s, over
  every coefficient ring, compatibly with inclusions.
* Degree 1 is, by Alexander duality in S², the degree-0 persistence of the complement, swept downwards. That is
  the filtered graph G1: its vertices are the pixels, the absent pixels and one node for the outside; its edges
  join 4-adjacent cells and border cells to the outside. Over any coefficient ring,
  H̃_1(X_s) ≅ H̃^0(S² − X_s), naturally in s.
* The union-finds are the standard reductions of the boundary matrices ∂_1 of G0 and G1. Each column has
  entries +1 and −1 in the rows of its two ends.

**Incidence matrices are totally unimodular.** The incidence matrix of a directed graph is totally unimodular:
every square minor is 0 or ±1 (Poincaré; Schrijver, *Theory of linear and integer programming*, §19.3). Every
pivot value of a reduction is a ratio of minors of the blocks of nested pairs (§1, "Pivot values are
invariants"). Here those are
±1, so no prime divides a pivot value. The rank of every block is the same over Q and over every F_p. The
barcode is therefore the same over Q and every F_p, and the set of bad primes is empty. The Q barcode itself is
exact: the union-find never forms a non-unit coefficient.

The same conclusion has a homological proof. For a filtration of subcomplexes of the plane, every relative
group H_q(X, A; Z) is torsion-free. By Alexander duality, H_1(A) → H_1(X) is dual to restriction
H̃^0(S² − A) → H̃^0(S² − X). That restriction is the pullback of locally constant functions along the map of
component sets, and its cokernel is free. H_2(X, A) ≅ ker(H_1(A) → H_1(X)) is free, and so is H_0(X, A).

Periodic 2D grids (a torus, cylinder or annulus) are not planar. They go to the classic engine, whose
certificate does the arithmetic of §2. A torus has no torsion, but its essential classes are not given by
duality in S². Three-dimensional inputs can have torsion: a (2,1) curve on the core of a solid torus gives Z/2.

## 6. Chain complexes over Z, magnitude homology and Laurent coefficients

**Integral outputs** (`--integral`, `--image`, `integral=True`, `image=(s, t)`). The integral homology of the whole
complex and the image of H_*(K_s; Z) → H_*(K_t; Z) are exact: Smith forms and lattice quotients in int64 with
overflow checks, repeated with bigints when a check fails. Nothing is random. They are computed on the residual
complex of the algebraic Morse reduction (`TECHNIQUES.md` 5.2), which is filtered chain homotopy equivalent to
the input over Z; the tests check on random complexes that the reduction changes neither the groups nor any
image. The F_p barcodes of the bad primes are reruns over F_p of the primes the Q certificate found, so they
carry the Q certificate's guarantee (and its seed, in fast mode). Primes that the certificate only knows as
"every prime factor of N" get no separate barcode; their intervals are in `bad_primes` as before.

**Magnitude homology.** MH_{k,l}(G; Z) is computed exactly: the complexes are enumerated in full, reduced over Z
by unit pairs, and the residual is put in Smith form. The torsion is exact, not inferred from ranks over fields.

**Laurent coefficients.** Lemma 2 holds over any discrete valuation ring O with fraction field K and residue
field k, with the same proof: if every pivot value of a reduction over K is a unit in O, the reduction is defined
over O and its reduction to k is a valid reduction with the same pivots. Perseus uses it twice.

* O = Q[t]_(f), f ∈ Q[t] irreducible, f ≠ t: residue field Q[t]/(f). A pivot value is a unit there unless f
  divides its numerator or its denominator. So the bad irreducibles are among the irreducible factors of
  the numerators and denominators of the pivot values over Q(t). Like the rational pivot values, these are
  ratios of minors of the filtered matrix, independent of the reduction.
* O = Z[t]_(p), the polynomials over Q whose denominators have content prime to p: residue field F_p(t).
  A pivot value c·a/b, with a and b primitive integer polynomials and c ∈ Q, is a unit unless p divides the
  numerator or the denominator of c (Gauss's lemma). So the bad primes of the barcode over F_p(t) are among the
  prime factors of those contents.

The polynomial candidates are not factored into irreducibles. Their gcd-free basis is computed, with linear factors
split off by the rational root test. Each basis element g is settled by a reduction over Q[t]/(g), run as if it
were a field. This is the analogue of Lemma 3. If a pivot is not invertible modulo g, the gcd splits g and the
parts are settled separately. If the run finishes, then for every irreducible factor f of g it is a valid run
over Q[t]/(f) with the same pivots, so all of them have the barcode it found. An entry with `"irreducible": true`
(degree 1, or degree 2 with a non-square discriminant) names one irreducible. `"every_irreducible_factor_of": true`
means the answer holds for every irreducible factor of the polynomial given.

For the homology of the whole complex over Q[t, t^-1], the Smith form is computed by row and column operations
that multiply by nonzero integers and powers of t, then divide by integer contents. If p divides none of these
integers and none of the contents of the final pivots, the same operations are invertible over Z_(p)[t, t^-1].
Then H_*( ; F_p[t, t^-1]) is the reduction mod p of the rational answer: the same free ranks, and the invariant
factors (primitive integer polynomials, so none vanishes mod p, and the divisibility chain survives) reduced mod p,
with units dropped. Every prime dividing a recorded integer is a candidate. Each candidate is settled by a Smith
form over F_p[t]. `homology_bad_primes` lists those where the answer is not the reduction. `complete` is false
only if some recorded integer could not be factored (`unresolved`).

What this is not: a computation of the Z[t, t^-1]-module. That ring has dimension 2, and its modules are not
determined by their fibres over Q and the F_p. Nor is it a computation over the residue fields F_p[t]/(h) or
over the Novikov ring Z((t)).
