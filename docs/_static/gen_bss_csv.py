import pandas as pd
import numpy as np

df = pd.read_csv("cluster.csv")


def add_synthetic_members(df, n_new=30, seed=0):
    rng = np.random.default_rng(seed)

    # cluster statistics
    ra_med, ra_std = df["RA_ICRS"].median(), df["RA_ICRS"].std()
    de_med, de_std = df["DE_ICRS"].median(), df["DE_ICRS"].std()

    plx_med, plx_std = df["Plx"].median(), df["Plx"].std()
    pmra_med, pmra_std = df["pmRA"].median(), df["pmRA"].std()
    pmde_med, pmde_std = df["pmDE"].median(), df["pmDE"].std()

    rv_med = df["RV"].median(skipna=True)
    rv_std = df["RV"].std(skipna=True)

    g_min = df["Gmag"].min()
    bp_rp_min = df["BP-RP"].min()

    # generate synthetic stars
    new_rows = pd.DataFrame({
        "Source": [f"syn_{i}" for i in range(n_new)],

        "RA_ICRS": rng.normal(ra_med, 0.5 * ra_std, n_new),
        "DE_ICRS": rng.normal(de_med, 0.5 * de_std, n_new),

        "Plx": rng.normal(plx_med, 0.5 * plx_std, n_new),
        "e_Plx": rng.uniform(0.01, 0.05, n_new),

        "pmRA": rng.normal(pmra_med, 0.5 * pmra_std, n_new),
        "e_pmRA": rng.uniform(0.01, 0.05, n_new),

        "pmDE": rng.normal(pmde_med, 0.5 * pmde_std, n_new),
        "e_pmDE": rng.uniform(0.01, 0.05, n_new),

        "RV": rng.normal(rv_med, 0.5 * rv_std, n_new),
        "e_RV": rng.uniform(0.5, 2.0, n_new),

        # brighter and bluer sequence
        "Gmag": rng.normal(g_min + 1, 2, n_new),
        "BP-RP": rng.normal(bp_rp_min - 0.5, 0.1, n_new),

        "e_Gmag": rng.uniform(0.001, 0.01, n_new),
        "e_BP-RP": rng.uniform(0.002, 0.02, n_new),
    })

    return pd.concat([df, new_rows], ignore_index=True)

N = 20

# Remove N random entries
rng = np.random.default_rng(42)
drop_indices = rng.choice(df.index, size=N, replace=False)
df = df.drop(index=drop_indices).reset_index(drop=True)


df = add_synthetic_members(df, n_new=N)

# # Drop Source column
# df = df.drop(columns=["Source"])



df.to_csv("cluster_bss.csv", index=False)
