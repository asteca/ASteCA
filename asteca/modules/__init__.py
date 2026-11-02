# numba is optional: modules check HAS_NUMBA and fall back to their NumPy
# implementations when it is not installed
try:
    from numba import njit

    HAS_NUMBA = True
except ImportError:
    njit = None
    HAS_NUMBA = False
