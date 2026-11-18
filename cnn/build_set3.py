"""Draw the third test set: random windows at the three sites results are needed for.

    python build_set3.py

  - GC_15, MC_17 and DC_15, 20 five-minute windows each.
  - Random within each site. No detector takes part in choosing them, so the
    score reflects what a run over the whole deployment would give.
  - Only days not used by any labelled set (set 1, set 2, train1, train2),
    and at most one window per day.
  - GC_15: only blocks with normal recording level (the level collapses after
    mid-October 2022; those windows were dropped as bad data in train1).

Writes set3/windows.csv and the network's detections on those windows (average
of the train12_seed* copies) for the starting marks. The original and the
shape-based detector are then run by score_detector.m and all three are pooled
by make_proposals.py (see build_set3.sh). The set is for testing only: nothing
is to be trained or tuned on it.
"""
import csv
import os

import numpy as np
import scipy.io as sio
import torch

import data
import ensemble
import score
from select_windows import datenum_to_dt

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = data.TESTS
OUT = os.path.join(TESTS, "set3")
PER_SITE = 20
PROPOSAL_THR = 0.30                      # lenient, so faint shots get a starting mark
rng = np.random.default_rng(20261118)
#        site     scan file              minimum background level (dB counts)
SITES = [("GC_15", "scan_GOM_GC_15.mat", 20.0),
         ("MC_17", "scan_GOM_MC_17.mat", None),
         ("DC_15", "scan_GOM_DC_15.mat", None)]


def used_days():
    used = set()
    for folder in ("set1", "set2", "train1", "train2"):
        for r in csv.DictReader(open(os.path.join(TESTS, folder, "windows.csv"))):
            used.add((r["site"], r["start_time"][:10]))
    return used


def main():
    models, _ = ensemble.load("train12_seed")
    used = used_days()
    chosen = []
    for site, scan, min_bg in SITES:
        m = sio.loadmat(os.path.join(TESTS, "scans", scan), squeeze_me=True, struct_as_record=False)
        R, fl = m["R"], m["fl"]
        ok = R[:, -1] == data.SEG
        if min_bg is not None:
            ok &= R[:, 6] >= min_bg
        days = np.array([datenum_to_dt(d).strftime("%Y-%m-%d") for d in R[:, 2]])
        ok &= np.array([(site, d) not in used for d in days])
        paths = data.file_index({site})
        taken, n = set(), 0
        for i in rng.permutation(np.where(ok)[0]):
            if days[i] in taken:
                continue
            name = fl[int(R[i, 0]) - 1].name
            w = dict(site=site, file=name, path=paths[name], seg=int(R[i, 1]),
                     start=datenum_to_dt(R[i, 2]).strftime("%Y-%m-%d %H:%M:%S"))
            try:
                S = data.spectrogram(data.read_window(w))
            except Exception as err:                      # unreadable stretch: skip it
                print("  skipped %s seg %d: %s" % (name, w["seg"], err))
                continue
            with torch.no_grad():
                x = torch.from_numpy(S).unsqueeze(0)
                p = np.mean([torch.sigmoid(mm(x))[0].numpy() for mm in models], axis=0)
            w["det"] = score.pick(p, PROPOSAL_THR)
            chosen.append(w)
            taken.add(days[i])
            n += 1
            if n == PER_SITE:
                break
        print("%s: %d windows from %d eligible blocks on %d unused days" %
              (site, n, int(ok.sum()), len(set(days[ok]))))

    os.makedirs(OUT, exist_ok=True)
    chosen.sort(key=lambda w: (w["site"], w["start"]))
    with open(os.path.join(OUT, "windows.csv"), "w", newline="") as f:
        c = csv.writer(f)
        c.writerow(["window", "condition_guess", "site", "file", "start_segment", "start_time", "proposed_shots"])
        for k, w in enumerate(chosen):
            w["window"] = "Y%02d" % (k + 1)
            c.writerow([w["window"], "random", w["site"], w["file"], w["seg"], w["start"], 0])
    with open(os.path.join(OUT, "detections_proposals_cnn.csv"), "w", newline="") as f:
        c = csv.writer(f)
        c.writerow(["window", "start_s", "end_s"])
        for w in chosen:
            c.writerows([[w["window"], "%.2f" % a, "%.2f" % b] for a, b in w["det"]])
    with open(os.path.join(OUT, "truth.csv"), "w", newline="") as f:
        f.write("window,time_s,source\n")
    print("%d windows written to %s" % (len(chosen), OUT))


if __name__ == "__main__":
    main()
