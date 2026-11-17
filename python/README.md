# perseus-tda

Python interface to [Perseus 6](https://people.maths.ox.ac.uk/nanda/perseus/): persistent homology over Q,
together with every prime p for which coefficients in F_p give a different barcode, and a certificate
that the list is complete. Or simply over the prime field of your choice.

```python
import numpy as np, perseus
r = perseus.rips(points=np.random.rand(500, 3), maxdim=2)
r["barcode"][1]                  # H1 intervals, death = inf if essential
r["certificate"]["mode"]         # "proved" (exact arithmetic) or "fast" (random primes, error bound reported)
perseus.bad_primes(r)            # the primes for which F_p gives a different barcode
```

Also `perseus.cubical`, `perseus.flag`, `perseus.simplices` and `perseus.cells`. See the web page for details.
