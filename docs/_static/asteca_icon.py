import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

cmd_xy = []
with open("asteca_icon.dat", "r") as f:
    for line in f:
        line_split = line.split()
        cmd_xy.append([float(line_split[1]), float(line_split[0])])
cmd_xy = np.array(cmd_xy).T

fig = plt.figure(figsize=(25, 25))
ax = fig.add_subplot(111, aspect="auto")
ax.set_axis_off()
plt.scatter(cmd_xy[0], cmd_xy[1], c="grey", s=300, edgecolor="none")
e = Ellipse(
    xy=(0.77, 19.3), width=4.5, height=9, angle=0.0, lw=45, fc="None", edgecolor="k"
)
ax.add_artist(e)
plt.xlim(-3.0, 6)
plt.ylim(28, 10)
plt.savefig("asteca_icon.png", dpi=150, bbox_inches="tight", transparent=True)


# # Favicon
# fig = plt.figure(figsize=(10, 10))
# ax = fig.add_subplot(111, aspect="auto")
# ax.set_axis_off()
# plt.scatter(cmd_xy[0], cmd_xy[1], c="k", s=300, edgecolor="none")
# xc = (cmd_xy[0].min() + cmd_xy[0].max()) * 0.5
# yc = (cmd_xy[1].min() + cmd_xy[1].max()) * 0.5
# # 0.68, 20.37
# plt.xlim(xc - 1.5, xc + 1.5)
# plt.ylim(yc + 4, yc - 4)
# # plt.gca().invert_yaxis()
# # plt.show()
# plt.savefig("asteca_icon.png", dpi=300, bbox_inches="tight", transparent=True)
