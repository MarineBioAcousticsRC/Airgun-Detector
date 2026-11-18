"""Average several trained networks and score the average.

    python ensemble.py train12_seed set1 set2
    python ensemble.py train1_seed set1 set2 --threshold 0.45

The first argument is the start of the model names (all cnn/models/<prefix>*.pt
are averaged). The probabilities of the copies are averaged before peaks are
picked. The threshold is the mean of the thresholds each copy chose on its own
held-back training windows, rounded to 0.05, unless --threshold is given.
Writes cnn/models/<prefix>avg_detections_<set>.csv and prints the scores.
"""
import argparse
import csv
import glob
import os

import numpy as np
import torch

import data
import score
from model import ShotNet

HERE = os.path.dirname(os.path.abspath(__file__))


def load(prefix):
    models, thr = [], []
    for p in sorted(glob.glob(os.path.join(HERE, "models", prefix + "*.pt"))):
        ck = torch.load(p)
        m = ShotNet(ck["n_freq"])
        m.load_state_dict(ck["state"])
        m.eval()
        models.append(m)
        thr.append(ck["threshold"])
    return models, thr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prefix")
    ap.add_argument("sets", nargs="+")
    ap.add_argument("--threshold", type=float, default=None)
    a = ap.parse_args()

    models, thrs = load(a.prefix)
    thr = a.threshold if a.threshold is not None else round(float(np.mean(thrs)) / 0.05) * 0.05
    print("%d networks averaged (own thresholds %s), threshold used %.2f" %
          (len(models), ", ".join("%.2f" % t for t in thrs), thr))
    for s in a.sets:
        ws = data.cache_set(s)
        dets = {}
        with torch.no_grad():
            for w in ws:
                x = torch.from_numpy(w["S"]).unsqueeze(0)
                p = np.mean([torch.sigmoid(m(x))[0].numpy() for m in models], axis=0)
                dets[w["window"]] = score.pick(p, thr)
        tot, rows = score.score_set(ws, dets)
        print("RESULT %-5s %s" % (s, score.describe(tot)))
        print("DETAIL " + " | ".join("%s %d/%d s%d" % (r["window"], r["found"], r["shots"], r["spare"]) for r in rows))
        with open(os.path.join(HERE, "models", "%savg_detections_%s.csv" % (a.prefix, s)), "w", newline="") as f:
            c = csv.writer(f)
            c.writerow(["window", "start_s", "end_s"])
            for w in ws:
                c.writerows([[w["window"], "%.2f" % d0, "%.2f" % d1] for d0, d1 in dets[w["window"]]])


if __name__ == "__main__":
    main()
