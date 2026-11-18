"""Figure for the report: what goes into the network, what it is made of, and what comes out.

    python plot_network.py [set] [window] [start_s] [length_s]

Top:    the input, a normalised spectrogram of a real window, with the reviewed shot marks.
Middle: the layers of the network, with the size of the data after each one.
Bottom: the output, the probability of a shot start every 0.1 s (average of the
        train12_seed* copies), the threshold, and the peaks taken as detections.
Writes figures/M_network.png.
"""
import os
import sys

import numpy as np
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "cnn"))
import data          # noqa: E402
import score         # noqa: E402
import ensemble      # noqa: E402

SET = sys.argv[1] if len(sys.argv) > 1 else "set1"
WIN = sys.argv[2] if len(sys.argv) > 2 else "W05"
T0 = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
LEN = float(sys.argv[4]) if len(sys.argv) > 4 else 120.0
MARK, DET, INK, MUTED = "#D55E00", "#0072B2", "#222222", "#666666"
FILL2D, FILL1D, FILLIO = "#DCEBF5", "#DDF0E6", "#F1F1F1"


def box(ax, x, y, w, h, title, lines, fill):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.004,rounding_size=0.012", fc=fill, ec="#8A8A8A", lw=1))
    ax.text(x + w / 2, y + h - 0.06, title, ha="center", va="top", fontsize=9.5, fontweight="bold", color=INK)
    ax.text(x + w / 2, y + h - 0.25, "\n".join(lines), ha="center", va="top", fontsize=8.3, color=INK, linespacing=1.45)


def arrow(ax, x0, x1, y):
    ax.annotate("", xy=(x1, y), xytext=(x0, y), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.3))


def main():
    w = [q for q in data.cache_set(SET) if q["window"] == WIN][0]
    models, thrs = ensemble.load("train12_seed")
    thr = round(float(np.mean(thrs)) / 0.05) * 0.05
    with torch.no_grad():
        x = torch.from_numpy(w["S"]).unsqueeze(0)
        p = np.mean([torch.sigmoid(m(x))[0].numpy() for m in models], axis=0)
    det = score.pick(p, thr)[:, 0]
    t = np.arange(w["S"].shape[1]) * data.HOP_S
    marks = w["marks"]
    inside = lambda v: v[(v >= T0) & (v <= T0 + LEN)]

    fig = plt.figure(figsize=(11, 9.2))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.25, 1.35, 0.9], hspace=0.42, left=0.07, right=0.93, top=0.955, bottom=0.065)

    # ---- input ----
    ax = fig.add_subplot(gs[0])
    im = ax.imshow(w["S"] * 20.0, origin="lower", aspect="auto", cmap="viridis", vmin=-10, vmax=40,
                   extent=[t[0], t[-1] + data.HOP_S, data.FMIN, data.FMAX], interpolation="nearest")
    ax.plot(inside(marks), np.full(len(inside(marks)), 240.0), "v", ms=7, mfc=MARK, mec="w", mew=0.8, clip_on=False)
    ax.set_xlim(T0, T0 + LEN)
    ax.set_ylabel("Frequency (Hz)")
    ax.set_title("1. Input: spectrogram, 99 frequency rows (5 to 250 Hz) by one column every 0.1 s",
                 loc="left", fontsize=11, fontweight="bold")
    cb = fig.colorbar(im, ax=ax, pad=0.012, fraction=0.03)
    cb.set_label("dB above the median\nlevel of each frequency", fontsize=8.5)
    ax.text(0.005, 0.04, "triangles: shot marks from the reviewer (not given to the network when it runs)",
            transform=ax.transAxes, fontsize=8.3, color="w")

    # ---- architecture ----
    ax = fig.add_subplot(gs[1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("2. Network: about 208,000 learned values", loc="left", fontsize=11, fontweight="bold")
    y, h = 0.17, 0.76
    xs = [0.004, 0.152, 0.436, 0.586, 0.862]
    ws = [0.113, 0.250, 0.115, 0.242, 0.134]
    box(ax, xs[0], y, ws[0], h, "Spectrogram", ["1 layer", "99 rows", "T columns", "", "(T = 3,001 for", "5 minutes)"], FILLIO)
    box(ax, xs[1], y, ws[1], h, "Shape-finding part",
        ["4 blocks, each: 3 x 3 filters over", "frequency and time, then", "frequency rows are pooled",
         "", "features:  16 > 32 > 48 > 64", "rows:  99 > 49 > 24 > 12 > 4", "time columns kept: T"], FILL2D)
    box(ax, xs[2], y, ws[2], h, "Combine", ["64 features x", "4 rows = 256", "per time step", "", "reduced to 64", "per time step"], FILLIO)
    box(ax, xs[3], y, ws[3], h, "Context part",
        ["7 blocks along time only:", "5-wide filters whose taps are", "spaced 0.1, 0.2, 0.4, 0.8,", "1.6, 3.2 and 6.4 s apart",
         "", "each output can use about", "25 s of sound on either side"], FILL1D)
    box(ax, xs[4], y, ws[4], h, "Output layer", ["64 features > 1", "value per time", "step, squashed", "to 0 - 1", "", "T probabilities"], FILLIO)
    for k in range(4):
        arrow(ax, xs[k] + ws[k] + 0.004, xs[k + 1] - 0.004, y + h / 2)
    ax.text(xs[1] + ws[1] / 2, 0.09, "learns what one shot looks like", ha="center", fontsize=8.5, color=MUTED, style="italic")
    ax.text(xs[3] + ws[3] / 2, 0.09, "learns rhythm and what surrounds a shot", ha="center", fontsize=8.5, color=MUTED, style="italic")
    ax.text(0.5, -0.02, "Training: each reviewed mark becomes a 0.4 s wide bump on a target curve; the learned values are adjusted until "
            "the output matches it.", ha="center", fontsize=8.5, color=INK)

    # ---- output ----
    ax = fig.add_subplot(gs[2])
    ax.plot(t, p, color=INK, lw=1.3)
    ax.axhline(thr, color=MUTED, lw=1, ls="--")
    ax.text(T0 + 0.24 * LEN, thr + 0.03, "threshold %.2f" % thr, ha="center", va="bottom", fontsize=8.5, color=MUTED)
    d = inside(det)
    ax.plot(d, p[np.round(d / data.HOP_S).astype(int)], "o", ms=7, mfc=DET, mec="w", mew=0.8, label="detection (peak above threshold)")
    ax.plot(inside(marks), np.full(len(inside(marks)), 1.08), "v", ms=7, mfc=MARK, mec="w", mew=0.8, clip_on=False,
            label="reviewer's shot mark")
    ax.set_xlim(T0, T0 + LEN)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Time in window (s)")
    ax.set_ylabel("Probability")
    ax.set_title("3. Output: probability that a shot starts, every 0.1 s",
                 loc="left", fontsize=11, fontweight="bold", pad=16)
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.10), ncol=2, frameon=False, fontsize=8.5, borderaxespad=0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    # line the time axis up with the spectrogram above (which has a colourbar)
    top = fig.axes[0].get_position()
    pos = ax.get_position()
    ax.set_position([top.x0, pos.y0, top.width, pos.height])

    os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
    out = os.path.join(HERE, "figures", "M_network.png")
    fig.savefig(out, dpi=170)
    print("saved %s  (%s %s, %d marks and %d detections shown)" % (out, SET, WIN, len(inside(marks)), len(d)))


if __name__ == "__main__":
    main()
