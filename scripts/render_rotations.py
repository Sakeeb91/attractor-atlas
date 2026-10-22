"""Render compact, looping rotation videos of four of the five certified attractors.

One full orbit per attractor (Aurelia, Naiad, Cassiopea, Mobula; the names in
SPECS), encoded with ffmpeg. Three outputs per attractor:

    gallery/<name>_orbit.mp4    h264, committed
    gallery/<name>_orbit.webm   vp9, committed
    gallery/<name>_orbit.gif    palettegen/paletteuse fallback, gitignored
                                (regenerable; the PillowWriter route produced
                                ~6 MB files, this one is a fraction of that)

Trinacria (the fifth specimen, D3) has no orbit video yet: it is not in SPECS.
Adding one is separate content work for that specimen.

Usage: python scripts/render_rotations.py [names...]   (default: every name in SPECS)
"""

from __future__ import annotations

import pathlib
import subprocess
import sys
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
GALLERY = ROOT / "gallery"

from aurelia.dynamics import trajectory as aurelia_traj
from cassiopea.dynamics import trajectory as cassiopea_traj
from mobula.dynamics import trajectory as mobula_traj
from naiad.dynamics import trajectory as naiad_traj
from scripts.render_cassiopea import CASSIOPEA_CMAP
from scripts.render_gallery import AURELIA_CMAP
from scripts.render_mobula import MOBULA_CMAP
from scripts.render_naiad import NAIAD_CMAP

# name -> (trajectory, cmap, box_aspect, elev, base_azim, color_bias)
SPECS = {
    "aurelia":   (aurelia_traj,   AURELIA_CMAP,   (1, 1, 0.9),  26, -52, 0.14),
    "naiad":     (naiad_traj,     NAIAD_CMAP,     (1, 1, 1.15), 24, -60, 0.16),
    "cassiopea": (cassiopea_traj, CASSIOPEA_CMAP, (1, 1, 1.15), 28, -44, 0.22),
    "mobula":    (mobula_traj,    MOBULA_CMAP,    (1, 1, 1.1),  24, -52, 0.32),
}

N_FRAMES = 60
N_POINTS = 400_000
PX = 600
FPS = 18
ALPHA = 0.46
PT_SIZE = 0.085


def _limits(ax, pts, pad=0.04):
    for axis, col in zip("xyz", pts.T):
        lo, hi = col.min(), col.max()
        m = pad * (hi - lo)
        getattr(ax, f"set_{axis}lim")(lo - m, hi + m)


def render(name):
    traj, cmap, box, elev, base, bias = SPECS[name]
    print(f"[{name}] integrating {N_POINTS} points...", flush=True)
    pts = traj(n_steps=N_POINTS, dt=0.004, transient=30_000)
    t = (pts[:, 2] - pts[:, 2].min()) / (np.ptp(pts[:, 2]) or 1)
    c = bias + (1 - bias) * t

    fig = plt.figure(figsize=(PX / 100, PX / 100), facecolor="black")
    ax = fig.add_subplot(111, projection="3d", facecolor="black")
    ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], s=PT_SIZE, c=c, cmap=cmap,
               alpha=ALPHA, linewidths=0)
    ax.set_axis_off()
    ax.set_box_aspect(box, zoom=1.26)
    _limits(ax, pts)
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)

    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        for i in range(N_FRAMES):
            ax.view_init(elev=elev, azim=base + 360.0 * i / N_FRAMES)
            fig.savefig(td / f"f{i:03d}.png", dpi=100, facecolor="black")
        plt.close(fig)
        seq = str(td / "f%03d.png")
        # looping video (primary web format): mp4 (h264) + webm (vp9). Small, smooth,
        # and pausable under prefers-reduced-motion.
        mp4 = GALLERY / f"{name}_orbit.mp4"
        webm = GALLERY / f"{name}_orbit.webm"
        gif = GALLERY / f"{name}_orbit.gif"
        print(f"[{name}] encoding mp4...", flush=True)
        subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(FPS),
                        "-i", seq, "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "23",
                        "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-movflags", "+faststart",
                        str(mp4)], check=True)
        print(f"[{name}] encoding webm...", flush=True)
        subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(FPS),
                        "-i", seq, "-c:v", "libvpx-vp9", "-crf", "36", "-b:v", "0",
                        "-row-mt", "1", str(webm)], check=True)
        print(f"[{name}] encoding compact gif (fallback)...", flush=True)
        pal = td / "pal.png"
        subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(FPS),
                        "-i", seq, "-vf", "scale=460:-1:flags=lanczos,palettegen=stats_mode=full",
                        str(pal)], check=True)
        subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(FPS),
                        "-i", seq, "-i", str(pal),
                        "-lavfi", "scale=460:-1:flags=lanczos[x];[x][1:v]paletteuse=dither=bayer:bayer_scale=3",
                        "-loop", "0", str(gif)], check=True)
        for f in (mp4, webm, gif):
            print(f"[{name}] wrote {f.name}  ({f.stat().st_size/1024:.0f} KB)", flush=True)


def main():
    names = sys.argv[1:] or list(SPECS)
    for n in names:
        render(n)


if __name__ == "__main__":
    main()
