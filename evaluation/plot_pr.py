"""Precision-recall curves for the three detectors on set 1 and set 2.

    python plot_pr.py                 all three detectors
    python plot_pr.py --no-network    the two cross-correlation detectors only

Needs pr_curve.csv (written by sweep_pr.m: the original detector and the
prototype over a range of settings) and the trained networks in cnn/models.
The network curve comes from sweeping the detection threshold on the average
of the copies trained on train1 + train2 (train12_seed*); the average of the
copies trained on train1 alone is drawn as a dashed line for comparison.

    recall    = shots found / shots marked
    precision = detections that hold a shot / all detections

Writes pr_curve_cnn.csv and figures/PR_curves.png.
"""
import csv
import os
import sys

import numpy as np
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "cnn"))
import data          # noqa: E402
import score         # noqa: E402
import ensemble      # noqa: E402

SETS = [("set1", "Set 1 (22 windows, 313 shots)"), ("set2", "Set 2 (28 windows, 546 shots)")]
COL = {"original": "#D55E00", "prototype": "#0072B2", "cnn": "#009E73"}
THRS = np.round(np.arange(0.05, 0.96, 0.05), 2)


def cnn_curve(prefix):
    models, thrs = ensemble.load(prefix)
    if not models:
        return None, None
    used = round(float(np.mean(thrs)) / 0.05) * 0.05
    out = {}
    for s, _ in SETS:
        ws = data.cache_set(s)
        probs = {}
        with torch.no_grad():
            for w in ws:
                x = torch.from_numpy(w["S"]).unsqueeze(0)
                probs[w["window"]] = np.mean([torch.sigmoid(m(x))[0].numpy() for m in models], axis=0)
        pts = []
        for thr in THRS:
            dets = {k: score.pick(p, thr) for k, p in probs.items()}
            _, rows = score.score_set(ws, dets)
            pts.append(dict(thr=float(thr), shots=sum(r["shots"] for r in rows), detections=sum(r["detections"] for r in rows),
                            found=sum(r["found"] for r in rows), spare=sum(r["spare"] for r in rows)))
        out[s] = pts
    return out, used


def pr(p):
    rec = p["found"] / max(p["shots"], 1)
    prec = (p["detections"] - p["spare"]) / p["detections"] if p["detections"] else np.nan
    return rec, prec


def main():
    mat = {}
    for r in csv.DictReader(open(os.path.join(HERE, "pr_curve.csv"))):
        p = {k: int(float(r[k])) for k in ("shots", "detections", "found", "spare", "step", "current")}
        mat.setdefault((r["detector"], r["set"]), []).append(p)
    net = "--no-network" not in sys.argv
    new, thr_new = cnn_curve("train12_seed") if net else (None, None)
    old, thr_old = cnn_curve("train1_seed") if net else (None, None)
    main_cnn, main_thr, main_label = (new, thr_new, "Network (train1 + train2)") if new else (old, thr_old, "Network (train1)")

    with open(os.path.join(HERE, "pr_curve_cnn.csv" if net else os.devnull), "w", newline="") as f:
        c = csv.writer(f)
        c.writerow(["network", "set", "threshold", "current", "shots", "detections", "found", "spare", "recall", "precision"])
        for name, cur, used in (("train1 + train2", new, thr_new), ("train1", old, thr_old)):
            if cur is None:
                continue
            for s, _ in SETS:
                for p in cur[s]:
                    c.writerow([name, s, "%.2f" % p["thr"], int(abs(p["thr"] - used) < 1e-6), p["shots"], p["detections"],
                                p["found"], p["spare"], "%.3f" % pr(p)[0], "%.3f" % pr(p)[1]])

    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2), sharey=True)
    for ax, (s, title) in zip(axes, SETS):
        for det, label, mk in (("original", "Original", "o"), ("prototype", "Shape-based cross-correlation", "s")):
            pts = sorted(mat[(det, s)], key=lambda p: p["step"])
            xy = np.array([pr(p) for p in pts])
            ax.plot(xy[:, 0], xy[:, 1], "-", color=COL[det], lw=2, marker=mk, ms=4, label=label)
            k = [i for i, p in enumerate(pts) if p["current"]][0]
            ax.plot(xy[k, 0], xy[k, 1], mk, ms=12, mfc="none", mec=COL[det], mew=2)
        if new and old:
            xy = np.array([pr(p) for p in old[s]])
            ax.plot(xy[:, 0], xy[:, 1], "--", color=COL["cnn"], lw=1.5, label="Network (train1 only)")
        if net:
            xy = np.array([pr(p) for p in main_cnn[s]])
            ax.plot(xy[:, 0], xy[:, 1], "-", color=COL["cnn"], lw=2, marker="^", ms=4, label=main_label)
            k = int(np.argmin([abs(p["thr"] - main_thr) for p in main_cnn[s]]))
            ax.plot(xy[k, 0], xy[k, 1], "^", ms=13, mfc="none", mec=COL["cnn"], mew=2)
        ax.set_title(title, fontsize=11, loc="left")
        ax.set_xlabel("Recall (shots found / shots marked)")
        ax.set_xlim(0, 1.02)
        ax.set_ylim(0, 1.02)
        ax.grid(color="0.9", lw=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].set_ylabel("Precision (detections on a shot / all detections)")
    axes[0].plot([], [], "o", ms=11, mfc="none", mec="0.3", mew=2, label="Setting in use")
    axes[0].legend(loc="lower left", frameon=False, fontsize=9)
    fig.tight_layout()
    os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
    out = os.path.join(HERE, "figures", "PR_curves.png" if net else "PR_curves_xcorr.png")
    fig.savefig(out, dpi=170)
    print("saved", out)
    for s, _ in SETS:
        for det in ("original", "prototype"):
            p = [q for q in mat[(det, s)] if q["current"]][0]
            print("%s %-9s in use: recall %.2f precision %.2f" % (s, det, *pr(p)))
        if not net:
            continue
        p = main_cnn[s][int(np.argmin([abs(q["thr"] - main_thr) for q in main_cnn[s]]))]
        print("%s %-9s in use (threshold %.2f): recall %.2f precision %.2f" % (s, "network", main_thr, *pr(p)))


if __name__ == "__main__":
    main()
