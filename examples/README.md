# Getting started with Perseus: examples

Nine short examples, each in two forms: a shell script for the command-line program (`cli/`) and a Python
script for the `perseus` module (`python/`). They print what they compute and say what to look for. Every
input they read is in `data/`, and each runs in a second or two.

| | topic | shows |
|---|---|---|
| 01 | [point clouds](python/01_point_cloud.py) ([cli](cli/01_point_cloud.sh)) | Vietoris–Rips of points or a distance matrix; long and short intervals; output files; `--coeff`, `--threshold` |
| 02 | [bad primes](python/02_bad_primes.py) ([cli](cli/02_bad_primes.sh)) | the barcode over Q, every prime p for which F_p differs, and which intervals change (a Moore space, RP², a point cloud on RP²) |
| 03 | [images and volumes](python/03_images.py) ([cli](cli/03_images.sh)) | cubical complexes from numpy or text; sub- and superlevel sets; periodic axes; absent pixels; engines |
| 04 | [your own filtration](python/04_your_own_filtration.py) ([cli](cli/04_your_own_filtration.sh)) | simplicial complexes (a lower-star filtration on a torus) and clique filtrations of weighted graphs |
| 05 | [integral homology](python/05_integral_homology.py) ([cli](cli/05_integral_homology.sh)) | H_*(K; Z), the image of H_*(K_s; Z) → H_*(K_t; Z), and the torsion of H_1(K_t; Z) along a random 2-complex |
| 06 | [simplicial sets](python/06_simplicial_sets.py) ([cli](cli/06_simplicial_sets.sh)) | RP² with three simplices; nerves of categories; classifying spaces BG |
| 07 | [magnitude homology](python/07_magnitude.py) ([cli](cli/07_magnitude.sh)) | MH_{k,l} of graphs over Z, with torsion |
| 08 | [Laurent coefficients](python/08_laurent.py) ([cli](cli/08_laurent.sh)) | chain complexes over Z[t, t⁻¹]: the trefoil's Alexander polynomial, bad irreducibles, integral content |
| 09 | [finite groups](python/09_finite_groups.py) ([cli](cli/09_finite_groups.sh)) | H_*(G; Z) and mod-p cohomology rings (a pc group, A₅ through its Sylow subgroups), persistence along a series |

## Running them

With the executable from a [release](https://github.com/viditnanda/perseus/releases) on your `PATH` (or
built with `make`):

```
sh examples/cli/01_point_cloud.sh
PERSEUS=./perseus sh examples/cli/02_bad_primes.sh       # if perseus is not on PATH
```

With the Python module (`pip install perseus-tda`, which needs numpy):

```
python3 examples/python/01_point_cloud.py
```

`sh examples/run_all.sh` (or `make examples`) runs all of them and stops at the first failure.

## Reading the output

* On the command line, intervals go to stdout in Ripser's format, `[birth,death)`, with an empty death for
  intervals that never end, and a summary of the certificate goes to stderr:

  ```
  Q barcode: proved (exact arithmetic: int64)
  bad primes: 2  (2 non-unit pivots)
  ```

  With an output prefix, `perseus points cloud.txt out` writes `out_0.txt`, `out_1.txt`, … (one `birth death`
  per line, death `-1` if infinite), `out_betti.txt` and, over Q, `out_cert.json`.
* In Python every function returns a dictionary. `r["barcode"][k]` is an (m, 2) numpy array of degree-k
  intervals, with death `inf` if infinite; `r["certificate"]` holds the certificate (`"proof"`, the
  `"bad_primes"` with the intervals present only over F_p and only over Q, the non-unit pivots, timings);
  `perseus.bad_primes(r)` lists the bad primes.

## The inputs in `data/`

| file | format | what it is |
|---|---|---|
| `circle.txt` | `points` | 100 noisy points near the unit circle |
| `rp2_points.txt` | `points` | 250 points on the Veronese surface (RP²) in R⁶ |
| `rp2_simplices.txt` | `simplices` | the 6-vertex triangulation of RP² |
| `moore12_cells.txt` | `cells` | the Moore space M(Z/12, 1), filtered |
| `rp2_two_stage.cells` | `cells` | RP¹ at time 0, RP² at time 1 |
| `image_cubical.txt` | `cubical` | a 3 × 4 image as text, with an absent pixel |
| `blobs.npy`, `volume.npy` | `npy` | a 64 × 64 image and a 32³ volume |
| `torus_height.txt` | `simplices` | a triangulated torus, lower-star filtration of a height function |
| `square_flag.txt` | `flag` | a weighted 4-cycle with a late chord, labelled vertices |
| `lm40.txt` | `simplices` | a Linial–Meshulam process on 40 vertices, through its torsion |
| `rp2.sset` | `sset` | RP² as a simplicial set (one simplex per dimension) |
| `cycle7.txt`, `ky_rp2.txt` | `magnitude` | the 7-cycle; the Kaneta–Yoshinaga graph of RP² |
| `trefoil.laurent` | `laurent` | the trefoil's infinite cyclic cover over Z[t, t⁻¹] |
| `q16_pcgroup.txt`, `a5_sylow.txt`, `d64_lower_central.txt` | `pcgroup`, `sylow`, `pcseries` | Q₁₆, A₅ and a series of dihedral groups, written from GAP by `tools/export_*.g` |

`data/make_data.py` regenerates the generated ones (numpy only, fixed seeds). The full list of input formats
and options is in the [main README](../README.md#inputs) and `perseus --help`.
