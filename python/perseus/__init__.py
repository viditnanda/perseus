"""Perseus 6: persistent homology over Q with certified bad primes, or over any prime field.

    import perseus
    r = perseus.rips(points, maxdim=2)              # Q barcode + every bad prime
    r["barcode"][1]                                 # (m, 2) array of H1 intervals, death = inf if essential
    r["certificate"]["mode"]                        # "proved" (exact arithmetic) or "fast" (random primes)
    r["certificate"]["bad_primes"]                  # {p: {"only_over_Fp": ..., "only_over_Q": ..., ...}}
    perseus.rips(points, maxdim=2, coeff=2)         # just the F_2 barcode
    perseus.rips(points, maxdim=2, mode="proved")   # insist on exact arithmetic throughout
    perseus.cubical(image), perseus.flag([(0, 1, 0.5), ...]), perseus.simplices([((0, 1, 2), 3.0), ...]), perseus.cells([...])

Certificates: mode="auto" (default) runs the Q reduction in exact int64 arithmetic, which is a proof;
if that overflows it uses random 62-bit primes and reports a bound on the error probability
("log2_error_bound"). An entry of "bad_primes" with "every_prime_factor_of" set applies to every
prime factor of its key. See docs/CERTIFICATE.md.
"""
from ._perseus import rips, flag, cells, simplices, cubical

__all__ = ["rips", "flag", "cells", "simplices", "cubical", "bad_primes"]
__version__ = "6.0.0"


def bad_primes(result):
    """Sorted list of the bad primes (and numbers all of whose prime factors are bad) of a Q result."""
    return sorted(result.get("certificate", {}).get("bad_primes", {}).keys())
