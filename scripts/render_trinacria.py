"""Render the Trinacria attractor gallery.

A new palette for the atlas's first dihedral member: coral-ember (deep wine
to coral to cream), chosen for a fossil-sea diatom and distinct from Aurelia
gold, Naiad cyan, Cassiopea emerald and Mobula violet-rose.

Trinacria carries the chiral dihedral group D3: the 3-fold rotation about z
plus the half-turn kappa:(x,y,z)->(x,-y,-z), which swaps the upper and lower
halves and reverses their winding. Color runs with z so the two
kappa-exchanged halves read as warm-over-cool in every view.

Outputs in gallery/:
  trinacria_hero.png        three-quarter view
  trinacria_triskelion.png  top-down view of the 3-fold (six-petal) pattern
  trinacria_profile.png     side profile (the two counter-rotating halves)
"""

from __future__ import annotations

import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from trinacria.dynamics import trajectory

GALLERY = pathlib.Path(__file__).resolve().parents[1] / "gallery"

# Coral-ember: deep wine -> crimson -> coral -> apricot -> cream.
TRINACRIA_CMAP = LinearSegmentedColormap.from_list(
    "trinacria",
    ["#15060a", "#6e1230", "#d2304c", "#ff8c5a", "#ffe6c0"],
)


def _frame(ax, pts, pad=0.04):
    for axis, col in zip("xyz", pts.T):
        lo, hi = col.min(), col.max()
        m = pad * (hi - lo)
        getattr(ax, f"set_{axis}lim")(lo - m, hi + m)


# Point size and alpha are larger than Mobula's (0.06, 0.42): Trinacria's orbit
# lies on thin, widely spaced sheets, so at the Mobula settings the 2400 px
# frame came out nearly black (viewed and rejected).
def render(pts, elev, azim, path, size=12, dpi=200, ps=0.35, alpha=0.7):
    t = (pts[:, 2] - pts[:, 2].min()) / (np.ptp(pts[:, 2]) or 1)
    c = 0.30 + 0.70 * t
    fig = plt.figure(figsize=(size, size), facecolor="black")
    ax = fig.add_subplot(111, projection="3d", facecolor="black")
    ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], s=ps, c=c, cmap=TRINACRIA_CMAP,
               alpha=alpha, linewidths=0)
    ax.view_init(elev=elev, azim=azim)
    ax.set_axis_off()
    ax.set_box_aspect((1, 1, 1.0), zoom=1.25)
    _frame(ax, pts)
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    fig.savefig(path, dpi=dpi, facecolor="black")
    plt.close(fig)
    print(f"wrote {path}")


def main():
    GALLERY.mkdir(exist_ok=True)
    print("integrating 2.5M points...")
    pts = trajectory(n_steps=2_500_000, dt=0.004, transient=40_000)
    span = np.ptp(pts, axis=0)
    print(f"spans {span.round(3)}, iso {span.min() / span.max():.3f}")
    render(pts, elev=22, azim=-52, path=GALLERY / "trinacria_hero.png")
    render(pts, elev=87, azim=-90, path=GALLERY / "trinacria_triskelion.png")
    render(pts, elev=2, azim=-90, path=GALLERY / "trinacria_profile.png")


if __name__ == "__main__":
    main()
