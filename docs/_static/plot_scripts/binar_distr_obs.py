import matplotlib.pyplot as plt
import numpy as np

# Since the original data is not available, we will generate synthetic data
# that visually resembles the distribution shown in the histogram.
# We estimate the properties from the image:
# - The distribution looks roughly normal or slightly skewed right.
# - The peak (mode) is around 0.36.
# - The data ranges roughly from 0.20 to 0.62.
# - The standard deviation seems to be around 0.06 - 0.07.
# - The maximum frequency is around 30. Let's estimate the total number of samples.
#   Summing approximate bin heights: 5+2+8+7+20+21+30+28+24+14+11+8+4+2+5+2+3+1+1 ~ 230? Let's try N=200.

# Set a seed for reproducibility of the random data
np.random.seed(42)

# Generate data resembling the histogram
mean_val = 0.36
std_dev = 0.065
num_points = 200
data = np.random.normal(loc=mean_val, scale=std_dev, size=num_points)

# Optional: Clip data to match the approximate range observed, though hist range handles this
# data = np.clip(data, 0.20, 0.62)

# Create the histogram plot
fig, ax = plt.subplots(figsize=(8, 4), dpi=150) # Adjust figsize and dpi as needed

# Define bins based on visual inspection (approx 21 bins from 0.20 to 0.62, width=0.02)
bin_edges = np.linspace(0.20, 0.62, 22) # 22 edges make 21 bins

# Plot the histogram
ax.hist(data, bins=bin_edges, color='steelblue', alpha=0.55) # Adjust alpha for transparency

ax.axvline(np.median(data), c='r', ls='--', lw=3)
ax.axvline(np.median(data)-np.std(data), c='orange', ls=':', lw=3)
ax.axvline(np.median(data)+np.std(data), c='orange', ls=':', lw=3)
print(np.median(data), np.std(data))

# Set labels and title (matching the original image)
ax.set_xlabel(r'$b_{fr}$', fontsize=12) # Use LaTeX for subscript
# The original image does not have a y-axis label or a title.
# ax.set_ylabel('Frequency', fontsize=12)
# ax.set_title('Histogram of b_fr', fontsize=14)

# Set axis limits and ticks to closely match the original image
ax.set_xlim(0.19, 0.63)
ax.set_ylim(0, 32) # A bit higher than the max frequency
ax.set_xticks(np.arange(0.20, 0.61, 0.05))
ax.set_yticks(np.arange(0, 31, 5))

# Adjust layout
plt.tight_layout()

# Show the plot
plt.savefig('binar_distr_obs.webp', bbox_inches='tight')