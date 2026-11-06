import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gammaln
from math import factorial

n = np.arange(0, 21)   # observed stars per bin
m = np.arange(0, 21)   # synthetic stars per bin
N, M = np.meshgrid(n, m)

log_nfact = np.array([np.log(factorial(k)) for k in n])  # shape (21,)

# Per-bin log-likelihood (B=1 term)
log_p = (
    gammaln(N + M + 0.5)
    - gammaln(M + 0.5)
    - (N + M + 0.5) * np.log(2)
    - log_nfact[N]
)

fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(
    log_p,
    origin="lower",
    extent=[-0.5, 20.5, -0.5, 20.5],
    aspect="equal",
    cmap="viridis",
)
cb = fig.colorbar(im, ax=ax, label=r"$\log\,p_i$")
ax.set_xlabel(r"$n_i$  (observed)", fontsize=12)
ax.set_ylabel(r"$m_i$  (synthetic)", fontsize=12)
ax.set_xticks(np.arange(0, 21, 5))
ax.set_yticks(np.arange(0, 21, 5))
plt.tight_layout()
plt.savefig("../lkl_heatmap.webp", dpi=300)

