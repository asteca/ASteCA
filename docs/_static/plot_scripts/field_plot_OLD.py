import matplotlib.pyplot as plt
import pandas as pd

import asteca

df = pd.read_csv("field.csv")

# Generate a `cluster` object
my_cluster = asteca.cluster(
    ra=df['RA_ICRS'],
    dec=df['DE_ICRS'],
    magnitude=df["Gmag"],
    e_mag=df["e_Gmag"],
    color=df["BP-RP"],
    e_color=df['e_BP-RP'],
    pmra=df["pmRA"],
    e_pmra=df["e_pmRA"],
    pmde=df["pmDE"],
    e_pmde=df["e_pmDE"],
    plx=df["Plx"],
    e_plx=df["e_Plx"],
    verbose=2
)

my_cluster.get_center()

breakpoint()
# Access estimated center values
ra_c, dec_c = my_cluster.radec_c
pmra_c, pmde_c = my_cluster.pms_c
plx_c = my_cluster.plx_c

# Create figure and a custom GridSpec layout
fig = plt.figure(figsize=(10, 10))

# Create the first plot in the top-left
ax1 = fig.add_axes([0.1, 0.5, 0.35, 0.35])
ax1.scatter(my_cluster.ra, my_cluster.dec, c="k", alpha=0.25, s=5)
ax1.scatter(ra_c, dec_c, marker="x", s=50, c="r")
ax1.set_xlabel("dec")
ax1.set_ylabel("ra")

# Create the second plot in the top-right
ax2 = fig.add_axes([0.55, 0.5, 0.35, 0.35])
ax2.scatter(my_cluster.pmra, my_cluster.pmde, c="k", alpha=0.25, s=5)
ax2.scatter(pmra_c, pmde_c, marker="x", s=50, c="r")
ax2.set_xlabel("pmra")
ax2.set_ylabel("pmde")
ax2.set_xlim(-5, 0)
ax2.set_ylim(-2, 2)

# Create the third plot in the middle of the second row
ax3 = fig.add_axes([0.325, 0.1, 0.35, 0.35])  # Centered below
ax3.hist(my_cluster.plx, 150)
ax3.axvline(plx_c, c="r", ls=":", lw=3)
ax3.set_xlabel("plx")
ax3.set_xlim(-0.5, 2)

plt.tight_layout()
# plt.savefig("field.webp", bbox_inches="tight")
plt.savefig("field2.webp", bbox_inches="tight")

