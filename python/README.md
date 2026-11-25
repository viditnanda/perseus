# perseus-tda

Python interface to [Perseus 6](https://people.maths.ox.ac.uk/nanda/perseus/): persistent homology over Q,
together with every prime p for which coefficients in F_p give a different barcode, and a proof that the
list is complete. Or simply over the prime field of your choice. Beyond barcodes it computes integral
homology and torsion along filtrations, magnitude homology of graphs, homology over Z[t, t⁻¹], and the
homology and cohomology rings of finite groups.

```python
import numpy as np, perseus
r = perseus.rips(points=np.random.rand(500, 3), maxdim=2)
r["barcode"][1]                  # H1 intervals as an (m, 2) array, death = inf if essential
r["certificate"]["proof"]        # True: the barcode and the list of bad primes are proved
perseus.bad_primes(r)            # the primes for which F_p gives a different barcode
```

| function | input |
|---|---|
| `rips`, `flag` | point clouds, distance matrices, weighted graphs |
| `cubical` | images and volumes as numpy arrays (2D: planar engine; 3D, 4D: discrete gradient) |
| `simplices`, `cells`, `simplicial_set` | filtered simplicial complexes, chain complexes over Z, simplicial sets; `integral=True`, `image=(s, t)`, `torsion=True` |
| `nerve`, `bar_construction` | nerves of finite categories, classifying spaces of finite monoids and groups |
| `magnitude` | magnitude homology of graphs over Z |
| `laurent` | chain complexes over Z[t, t⁻¹] |
| `pc_group`, `finite_group`, `pc_series` | finite groups: H_*(G; Z), mod-p cohomology rings, persistence along series |

Getting-started examples for every function, with their data, are in
[examples/](https://github.com/viditnanda/perseus/tree/main/examples); the full reference is the
[README](https://github.com/viditnanda/perseus#readme) and the [web page](https://people.maths.ox.ac.uk/nanda/perseus/).
