# Scripts for the documentation images

Each image in `docs/_static/` used by the documentation is made by one of these
scripts. Run them from this folder (e.g. `python bss_region.py`): they read the data
and isochrones in `docs/_static/`, use the same style as the tutorials
(`../asteca.mplstyle`), and write the image directly to `docs/_static/`.

| Image(s)                                                       | Script           |
|----------------------------------------------------------------|------------------|
| `asteca_icon.webp`, `favicon.ico`                              | `asteca_icon.py` (uses `asteca_icon.dat`) |
| `binar_distr.webp`                                             | `binar_distr.py` |
| `bss_region.webp`                                              | `bss_region.py`  |
| `dr_models_geometry.webp`, `dr_models_distr.webp`, `dr_models_cmd.webp` | `dr_models.py` |
| `ext_law.webp`                                                 | `ext_law.py`     |
| `IMFs.webp`                                                    | `IMFs.py`        |
| `isoch_models.webp`                                            | `isoch_models.py` |
| `mag_correct.webp`                                             | `mag_correct.py` |
| `qdist_mass.webp`, `qdist_unif.webp`                           | `qdist.py`       |
| `synthplot.webp`                                               | `synthplot.py`   |

`gen_bss_csv.py` makes the data file `cluster_bss.csv` (used by the post-processing
tutorial), adding synthetic BSSs to `cluster.csv`. Running it again overwrites that
file.
