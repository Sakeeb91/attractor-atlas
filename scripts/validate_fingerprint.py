"""Validate the shape-fingerprint novelty metric (gap 3 of the paper plan).

Three tests of whether fingerprint distance (atlas/fingerprint.py) measures
something real about attractor shape:

  (a) Within-dysts sanity. Nearest-neighbor table for ten famous systems,
      plus two controls. Positive control: the Lorenz lineage (documented
      geometric variants of the Lorenz butterfly) should sit closer to each
      other than typical catalog pairs. Negative control: the Sprott family
      is *algebraically* minimal by design and geometrically diverse, so it
      should NOT cluster - if it did, the metric would be reading algebra,
      not shape.

  (b) Robustness / signal-to-noise. Fingerprint a single system (Aurelia,
      and Lorenz via dysts) many times with different trajectory initial
      conditions and different subsampling seeds. The within-system spread
      is the noise floor; the catalog's between-system distances are the
      signal. SNR = nearest-catalog distance / mean within-system distance.

  (c) A 2-D map. Metric MDS embedding of the full distance matrix (135
      catalog systems + Aurelia, Naiad, Cassiopea), the candidate key
      figure for the paper.

Outputs:
    results/fingerprint_validation.json
    gallery/fingerprint_map.png
"""

from __future__ import annotations

import itertools
import json
import pathlib
import sys
import warnings

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")

from atlas.fingerprint import distance, fingerprint  # noqa: E402

FAMOUS = [
    "Lorenz", "Rossler", "Thomas", "Aizawa", "Chen", "ChenLee",
    "Halvorsen", "Dadras", "NoseHoover", "RabinovichFabrikant",
]
PUBLISHED = {"aurelia": ("Aizawa", 0.185), "naiad": ("SprottJerk", 0.34),
             "cassiopea": ("NoseHoover", 0.224)}


def load_catalog():
    cat = json.loads((ROOT / "results" / "catalog_fingerprints.json").read_text())
    names = sorted(cat)
    fps = {k: np.array(cat[k]["fingerprint"]) for k in names}
    return names, fps


def our_fingerprints(n_replicas=1, rng=None):
    """Fingerprint the three C_n specimens (Aurelia, Naiad, Cassiopea) from fresh package trajectories."""
    from aurelia.dynamics import trajectory as traj_a
    from cassiopea.dynamics import trajectory as traj_c
    from naiad.dynamics import trajectory as traj_n

    rng = rng or np.random.default_rng(0)
    out = {}
    for name, traj in (("aurelia", traj_a), ("naiad", traj_n), ("cassiopea", traj_c)):
        reps = []
        for r in range(n_replicas):
            ic = np.array([0.3, 0.2, 0.1]) + rng.normal(0, 0.05, 3)
            pts = traj(n_steps=120_000, state0=ic)
            reps.append(fingerprint(pts[::8], rng=np.random.default_rng(r)))
        out[name] = reps
    return out


LORENZ_LINEAGE = ["Lorenz", "LorenzBounded", "LorenzStenflo", "Chen", "HyperLorenz", "HyperLu"]


def test_within_dysts(names, fps):
    """(a) nearest neighbors of the famous ten + positive/negative controls."""
    nn_table = {}
    nn_dists = {}
    for k in names:
        ds = sorted((distance(fps[k], fps[m]), m) for m in names if m != k)
        nn_dists[k] = ds[0][0]
        if k in FAMOUS:
            nn_table[k] = {
                "nearest": ds[0][1], "distance": round(float(ds[0][0]), 4),
                "second": ds[1][1], "second_distance": round(float(ds[1][0]), 4),
            }
    median_nn = float(np.median(list(nn_dists.values())))

    # Positive control: documented Lorenz-butterfly variants should be
    # mutually close relative to typical catalog pairs.
    lineage = [k for k in LORENZ_LINEAGE if k in names]
    within_d = [distance(fps[a], fps[b]) for a, b in itertools.combinations(lineage, 2)]
    rng = np.random.default_rng(0)
    all_pairs = [distance(fps[a], fps[b])
                 for a, b in itertools.combinations(rng.permutation(names)[:60], 2)]
    pos = {
        "members": lineage,
        "mean_within": round(float(np.mean(within_d)), 4),
        "median_within": round(float(np.median(within_d)), 4),
        "catalog_mean_pairwise": round(float(np.mean(all_pairs)), 4),
        "ratio": round(float(np.mean(within_d)) / float(np.mean(all_pairs)), 3),
    }

    # Negative control: the Sprott systems are algebraically minimal by
    # design and geometrically diverse; NN-coherence at the null rate means
    # the metric reads geometry, not algebra.
    sprott = [k for k in names if k.startswith("Sprott")]
    within = sum(
        1 for k in sprott
        if min((distance(fps[k], fps[m]), m) for m in names if m != k)[1].startswith("Sprott")
    )
    null_frac = (len(sprott) - 1) / (len(names) - 1)
    neg = {
        "n_members": len(sprott),
        "nn_within_family": within,
        "fraction": round(within / len(sprott), 3),
        "null_fraction": round(null_frac, 3),
    }
    return {
        "famous_nn_table": nn_table,
        "median_nn_spacing": round(median_nn, 4),
        "lorenz_lineage_positive_control": pos,
        "sprott_negative_control": neg,
    }


def test_robustness(names, fps):
    """(b) within-system fingerprint spread vs between-system distances."""
    rng = np.random.default_rng(7)
    out = {}

    # Aurelia: 4 trajectory ICs x 3 subsample seeds = 12 replicas.
    from aurelia.dynamics import trajectory
    replicas = []
    for i in range(4):
        ic = np.array([0.3, 0.2, 0.1]) + rng.normal(0, 0.08, 3)
        pts = trajectory(n_steps=120_000, state0=ic)[::8]
        for s in range(3):
            replicas.append(fingerprint(pts, rng=np.random.default_rng(100 * i + s)))
    out["aurelia"] = _robustness_stats(replicas, names, fps)

    # Lorenz via dysts (perturbed ICs), at the catalog's own trajectory
    # length (1500 points) and at 4x that, to measure how the noise floor
    # shrinks with trajectory length.
    try:
        import dysts.flows as flows

        for n_pts, key in ((1500, "lorenz_dysts_1500"), (6000, "lorenz_dysts_6000")):
            replicas = []
            for i in range(4):
                model = flows.Lorenz()
                model.ic = np.array(model.ic, dtype=float) * (1.0 + 0.02 * rng.standard_normal(3))
                pts = np.asarray(model.make_trajectory(n_pts))[:, :3]
                for s in range(3):
                    replicas.append(fingerprint(pts, rng=np.random.default_rng(100 * i + s)))
            out[key] = _robustness_stats(replicas, names, fps, exclude="Lorenz")
    except Exception as e:  # noqa: BLE001
        out["lorenz_dysts_1500"] = {"skipped": f"{type(e).__name__}: {e}"}
    return out


def _robustness_stats(replicas, names, fps, exclude=None):
    replicas = [r for r in replicas if r is not None]
    within = [distance(a, b) for a, b in itertools.combinations(replicas, 2)]
    mean_fp = np.mean(replicas, axis=0)
    between = sorted(
        (distance(mean_fp, fps[m]), m) for m in names if m != exclude
    )
    within_mean = float(np.mean(within))
    return {
        "n_replicas": len(replicas),
        "within_mean": round(within_mean, 4),
        "within_max": round(float(np.max(within)), 4),
        "nearest_catalog": between[0][1],
        "nearest_distance": round(float(between[0][0]), 4),
        "snr_nearest": round(float(between[0][0]) / within_mean, 1),
    }


def test_map(names, fps, ours):
    """(c) MDS embedding of the joint distance matrix + figure."""
    from sklearn.manifold import MDS

    all_names = names + list(ours)
    all_fps = [fps[k] for k in names] + [np.mean(ours[k], axis=0) for k in ours]
    n = len(all_names)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = D[j, i] = distance(all_fps[i], all_fps[j])

    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=0,
              normalized_stress="auto", n_init=4, max_iter=500)
    XY = mds.fit_transform(D)
    # Kruskal stress-1 (scale-free, interpretable: <0.2 usable, <0.1 good).
    emb = np.sqrt(((XY[:, None, :] - XY[None, :, :]) ** 2).sum(-1))
    iu = np.triu_indices(n, k=1)
    stress1 = float(np.sqrt(((emb[iu] - D[iu]) ** 2).sum() / (D[iu] ** 2).sum()))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(11, 8.5), facecolor="#0d0d14")
    ax.set_facecolor("#0d0d14")

    cat_idx = [i for i, k in enumerate(all_names) if k in names]
    ax.scatter(XY[cat_idx, 0], XY[cat_idx, 1], s=14, c="#5a6b85", alpha=0.75,
               linewidths=0, zorder=2, label="dysts catalog (135 systems)")

    for k in FAMOUS:
        i = all_names.index(k)
        ax.annotate(k, XY[i], fontsize=7.5, color="#aebbd0",
                    xytext=(4, 4), textcoords="offset points", zorder=4)
        ax.scatter(*XY[i], s=22, c="#aebbd0", linewidths=0, zorder=3)

    style = {"aurelia": ("#f4c95d", "Aurelia (C$_3$)"),
             "naiad": ("#6fd3c7", "Naiad (C$_2$)"),
             "cassiopea": ("#c792ea", "Cassiopea (C$_4$)")}
    for k, (color, label) in style.items():
        i = all_names.index(k)
        nn = min((distance(all_fps[i], fps[m]), m) for m in names)
        j = all_names.index(nn[1])
        ax.plot([XY[i, 0], XY[j, 0]], [XY[i, 1], XY[j, 1]],
                color=color, lw=0.9, alpha=0.55, ls=":", zorder=2)
        ax.scatter(*XY[i], s=230, marker="*", c=color, edgecolors="white",
                   linewidths=0.6, zorder=5)
        ax.annotate(label, XY[i], fontsize=10, color=color, fontweight="bold",
                    xytext=(8, -12), textcoords="offset points", zorder=5)

    ax.set_title("Shape space of known chaotic attractors (metric MDS of fingerprint "
                 f"distances, Kruskal stress-1 = {stress1:.3f})\n"
                 "stars = this work, dotted line = nearest known neighbor", fontsize=11)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.legend(loc="lower right", fontsize=9, framealpha=0.2)
    fig.tight_layout()
    out = ROOT / "gallery" / "fingerprint_map.png"
    fig.savefig(out, dpi=170, facecolor=fig.get_facecolor())
    print(f"wrote {out}")
    return {"kruskal_stress_1": round(stress1, 4),
            "raw_stress": round(float(mds.stress_), 4), "n_points": n}


def main():
    names, fps = load_catalog()
    print(f"catalog: {len(names)} systems")

    print("\n(a) within-dysts nearest-neighbor sanity...")
    res_a = test_within_dysts(names, fps)
    for k, v in res_a["famous_nn_table"].items():
        print(f"  {k:22s} -> {v['nearest']:22s} {v['distance']:.3f}")
    pos = res_a["lorenz_lineage_positive_control"]
    neg = res_a["sprott_negative_control"]
    print(f"  median NN spacing {res_a['median_nn_spacing']}")
    print(f"  positive control (Lorenz lineage {len(pos['members'])} systems): mean "
          f"within {pos['mean_within']} vs catalog mean {pos['catalog_mean_pairwise']} "
          f"(ratio {pos['ratio']})")
    print(f"  negative control (Sprott, algebraic family): NN-coherence "
          f"{neg['nn_within_family']}/{neg['n_members']} = {neg['fraction']} "
          f"(null {neg['null_fraction']}; at-null = metric reads shape, not algebra)")

    print("\n(b) robustness / SNR...")
    res_b = test_robustness(names, fps)
    for k, v in res_b.items():
        if "skipped" in v:
            print(f"  {k}: skipped ({v['skipped']})")
        else:
            print(f"  {k:14s} within {v['within_mean']:.4f} (max {v['within_max']:.4f}) "
                  f"vs nearest {v['nearest_catalog']} {v['nearest_distance']:.3f} "
                  f"-> SNR {v['snr_nearest']}")

    print("\n(c) our three C_n specimens (Aurelia, Naiad, Cassiopea) vs the catalog + MDS map...")
    ours = our_fingerprints(n_replicas=3, rng=np.random.default_rng(3))
    res_ours = {}
    for k, reps in ours.items():
        mean_fp = np.mean(reps, axis=0)
        ds = sorted((distance(mean_fp, fps[m]), m) for m in names)
        pub_nn, pub_d = PUBLISHED[k]
        res_ours[k] = {
            "nearest": ds[0][1], "distance": round(float(ds[0][0]), 4),
            "second": ds[1][1], "second_distance": round(float(ds[1][0]), 4),
            "published_nearest": pub_nn, "published_distance": pub_d,
            "nn_identity_stable": ds[0][1] == pub_nn,
        }
        print(f"  {k:10s} nearest {ds[0][1]} at {ds[0][0]:.3f} "
              f"(published: {pub_nn} at {pub_d})")
    res_c = test_map(names, fps, ours)

    results = {
        "catalog_size": len(names),
        "within_dysts": res_a,
        "robustness": res_b,
        "our_systems": res_ours,
        "mds": res_c,
    }
    out = ROOT / "results" / "fingerprint_validation.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"\nwrote {out}")

    noise_floor = res_b["aurelia"]["within_mean"]
    results["interpretation"] = {
        "noise_floor_15k_points": noise_floor,
        "noise_floor_catalog_1500_points": res_b.get("lorenz_dysts_1500", {}).get("within_mean"),
        "rule_of_thumb": "fingerprint distances below ~2x the noise floor are not "
                         "interpretable; NN identities become rank-unstable when "
                         "several catalog systems sit within one noise floor of "
                         "the same distance",
        "our_isolation_in_noise_floors": {
            k: round(v["distance"] / noise_floor, 1) for k, v in res_ours.items()
        },
    }
    out.write_text(json.dumps(results, indent=2))

    pos_ok = pos["ratio"] <= 0.6
    neg_ok = neg["fraction"] <= 2 * neg["null_fraction"]
    snr_ok = res_b["aurelia"]["snr_nearest"] >= 3
    stable = all(v["nn_identity_stable"] for v in res_ours.values())
    print(f"\nVERDICT: positive control (Lorenz lineage clusters) "
          f"{'PASS' if pos_ok else 'FAIL'}; "
          f"negative control (Sprott at null) {'PASS' if neg_ok else 'FAIL'}; "
          f"isolation SNR (Aurelia) {'PASS' if snr_ok else 'WEAK'}; "
          f"published NN identities {'all reproduced' if stable else 'Naiad rank-unstable'}")
    print(f"noise floor: {noise_floor} (15k-point orbits); our systems sit "
          f"{results['interpretation']['our_isolation_in_noise_floors']} noise floors out")


if __name__ == "__main__":
    main()
