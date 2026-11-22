"""ASteCA logo (asteca_icon.webp) and favicon (favicon.ico), from the CMD stored in
asteca_icon.dat (run from this folder)."""

import io

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse
from PIL import Image


def main():
    # Magnitude and color of the CMD
    mag, col = np.loadtxt("asteca_icon.dat", unpack=True)

    logo(col, mag)
    favicon(col, mag)


def logo(col, mag):
    """CMD in grey inside a black ellipse, with a transparent background"""
    fig = plt.figure(figsize=(25, 25))
    ax = fig.add_subplot(111)
    ax.set_axis_off()
    ax.scatter(col, mag, c="grey", s=300, edgecolor="none")
    e = Ellipse(
        xy=(0.77, 19.3), width=4.5, height=9, angle=0.0, lw=45, fc="None", edgecolor="k"
    )
    ax.add_artist(e)
    ax.set_xlim(-3.0, 6)
    ax.set_ylim(28, 10)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, transparent=True)
    plt.close(fig)

    # Crop to the visible part (the ring), and save at ~1550 px
    buf.seek(0)
    img = Image.open(buf).convert("RGBA")
    img = img.crop(img.getchannel("A").getbbox())
    img.thumbnail((1550, 1550))
    img.save("../asteca_icon.webp")


def favicon(col, mag):
    """CMD only, in black and zoomed in, with a transparent background. Rendered at a
    large size and downscaled to the favicon sizes"""
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    ax.set_axis_off()
    ax.scatter(col, mag, c="k", s=300, edgecolor="none")
    xc = (col.min() + col.max()) * 0.5
    yc = (mag.min() + mag.max()) * 0.5
    ax.set_xlim(xc - 1.5, xc + 1.5)
    ax.set_ylim(yc + 4, yc - 4)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=30, bbox_inches="tight", transparent=True)
    plt.close(fig)

    buf.seek(0)
    img = Image.open(buf).convert("RGBA")
    # Square canvas, so that the icon is not distorted
    side = max(img.size)
    square = Image.new("RGBA", (side, side), (0, 0, 0, 0))
    square.paste(img, ((side - img.width) // 2, (side - img.height) // 2))
    square.save("../favicon.ico", sizes=[(16, 16), (32, 32), (48, 48)])


if __name__ == "__main__":
    main()
