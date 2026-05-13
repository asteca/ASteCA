# ASteCA – Copilot Instructions

## Build, test, and lint

```bash
# Install in development mode
pip install -e .

# Quick import/version smoke check
python -c "import asteca; print(asteca.__version__)"

# Build docs (if docs are being changed)
pip install -r docs/requirements.txt
make -C docs html
```

There is currently no maintained top-level automated test suite in this checkout (no `tests/` directory), and no repository-configured lint command in `pyproject.toml`.

## Architecture

ASteCA is a Python package (Python ≥ 3.12) that automates stellar cluster analysis. The public API is five classes, all importable from `asteca`:

| Class | File | Role |
|---|---|---|
| `Cluster` | `asteca/cluster.py` | Container for observed star arrays (ra, dec, mag, color, pmra, pmde, plx, and their `e_*` uncertainties). Can be initialised from raw NumPy arrays or from a UCC members DataFrame. |
| `Isochrones` | `asteca/isochrones.py` | Loads and preprocesses theoretical isochrones (PARSEC, MIST, BASTI). |
| `Synthetic` | `asteca/synthetic.py` | Generates synthetic clusters from `Isochrones`; fits fundamental parameters. |
| `Membership` | `asteca/membership.py` | Estimates per-star membership probabilities (`bayesian`, `fastmp` algorithms). |
| `Likelihood` | `asteca/likelihood.py` | Scores how closely a synthetic cluster matches an observed one (Hess-diagram-based). |

All heavy implementation lives in `asteca/modules/`. Module files whose names end in `_priv.py` are the private backends for the corresponding public class (e.g. `synth_cluster_priv.py` backs `Synthetic`). Other modules are shared helpers (`imfs.py`, `stellar_stats_funcs.py`, `cluster_mass_funcs.py`).

Typical execution flow across classes is:

1. Load observed data in `Cluster`.
2. Estimate center and member count in `Cluster` (`get_center`, `get_nmembers`).
3. Run `Membership` methods (these depend on `Cluster` attributes created in step 2).
4. Load theoretical grids with `Isochrones`.
5. Create `Synthetic`, then `calibrate(cluster)` before higher-level inference (`get_models`, `stellar_parameters`, `cluster_masses`).
6. Compare observed vs synthetic data with `Likelihood`.

The `helpers/` directory contains standalone scripts for manual checks/profiling and a vendored SPISEA tree; it is not part of the installed package.

## Key conventions

### Public-class pattern
Each public class accepts a `verbose: int` parameter (default `1`). Output is controlled via the `_vp(msg, level)` helper defined on the class — level `0` messages always print at `verbose > 0`, level `1` at `verbose > 1`, etc. Passing `verbose=0` silences all output.

### Uncertainty columns are always paired
Every observable column has a matching uncertainty column prefixed with `e_`: `mag`/`e_mag`, `color`/`e_color`, `pmra`/`e_pmra`, `pmde`/`e_pmde`, `plx`/`e_plx`. Input validation raises `ValueError` whenever one is given without the other.

### Color handling stays split across `color` and `color2`
`Cluster` stores first and second color dimensions separately (`color`, `color2`) and `Synthetic` later combines these into list-based structures for downstream calculations. Keep this split/list transition consistent when changing data flow.

### pandas is an optional dependency
`pandas` is **not** in `dependencies`; use `from __future__ import annotations` and the `TYPE_CHECKING` guard when typing with `pd.DataFrame` to avoid a hard import at runtime.

### Private module naming
All functions that implement the internals of a public class live in `asteca/modules/` with the `_priv` suffix. The public class imports from these modules; do not put algorithm logic directly in the public class file.

### Isochrone model names
The three supported models are identified by the strings `"parsec"`, `"mist"`, and `"basti"` (case-insensitive at the API level; stored uppercase internally in `isochrones_priv.py`). The per-model metadata dict `phot_systs_data` in `isochrones_priv.py` is the single source of truth for column names, comment characters, and parsing flags.

### UCC cluster name normalisation
When loading from a UCC members file, cluster names are normalised by stripping spaces, underscores, hyphens, dots, and replacing `+` with `p` before lookup. Keep this logic consistent if adding new load paths.

### Membership prerequisites are attribute-based
`Membership` expects derived attributes already present in the `Cluster` object (especially `N_cluster`, and for methods also `radec_c`, `radius`, `pms_c`, `plx_c` as needed). Preserve this contract when changing cluster workflows.

### Synthetic lifecycle is stateful
`Synthetic.generate()` can run without calibration, but `get_models()`, `stellar_parameters()`, and `cluster_masses()` require calibration-driven attributes to exist. Preserve these guardrails and error messages when changing call flow.

### Recent API changes noted in repo docs
`pandas` was removed as a required dependency in v0.6.6, and `Synthetic.stellar_masses()` was renamed to `Synthetic.stellar_parameters()` in v0.6.5.
