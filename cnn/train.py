"""Train the spectrogram network and score it.

    python train.py --train set1 --test set2 --epochs 40 --name smoke

--train / --test take one or more set names (set1, set2, train1).
A fifth of the training windows are held back to choose the detection
threshold, so the test sets play no part in any choice.

--mix P (optional): for a fraction P of the training crops, a stretch of a
training window without airguns is added on top (as sound power, at a random
level), to make more examples of shots in ship noise and other loud background.
Writes cnn/models/<name>.pt and cnn/models/<name>_detections_<set>.csv.
"""
import argparse
import csv
import os
import time

import numpy as np
import torch
import torch.nn as nn

import data
import score
from model import ShotNet, count_parameters

HERE = os.path.dirname(os.path.abspath(__file__))


def mix_in(x, n, gain_dB):
    """Add background n (at gain_dB relative to x's own background) to x; both are normalised spectrograms."""
    g = 10.0 ** (gain_dB / 10.0)
    p = 10.0 ** (2.0 * x) + g * 10.0 ** (2.0 * n)          # spectrogram values are dB / 20
    out = 10.0 * np.log10(p / (1.0 + g))                   # back to dB above the new background
    return (np.clip(out, -10.0, 40.0) / 20.0).astype(np.float32)


def batches(windows, batch_size, crop, rng, mix=0.0):
    """Random time crops of the training windows."""
    quiet = [w for w in windows if len(w["marks"]) == 0]
    order = rng.permutation(len(windows))
    for i in range(0, len(order), batch_size):
        X, Y = [], []
        for j in order[i:i + batch_size]:
            S, y = windows[j]["S"], windows[j]["y"]
            a = rng.integers(0, S.shape[1] - crop + 1)
            x = S[:, a:a + crop].copy()
            if quiet and rng.random() < mix:            # add background from a window without airguns
                N = quiet[rng.integers(len(quiet))]["S"]
                b = rng.integers(0, N.shape[1] - crop + 1)
                x = mix_in(x, N[:, b:b + crop], rng.uniform(-6.0, 6.0))
            if rng.random() < 0.5:                      # hide a random band of frequencies
                f0 = rng.integers(0, S.shape[0] - 8)
                x[f0:f0 + rng.integers(2, 8)] = 0.0
            X.append(x)
            Y.append(y[a:a + crop])
        yield torch.from_numpy(np.stack(X)), torch.from_numpy(np.stack(Y))


def predict(model, windows):
    model.eval()
    out = {}
    with torch.no_grad():
        for w in windows:
            out[w["window"]] = torch.sigmoid(model(torch.from_numpy(w["S"]).unsqueeze(0)))[0].numpy()
    return out


def best_threshold(windows, probs):
    """Threshold that gives the most shots found minus spare detections on held-back training windows."""
    best = (None, -1e9)
    for thr in np.arange(0.2, 0.91, 0.05):
        tot, _ = score.score_set(windows, {w["window"]: score.pick(probs[w["window"]], thr) for w in windows})
        value = tot["found"] - tot["extra"] - tot["false"] - tot["merged"]
        if value > best[1]:
            best = (float(thr), value)
    return best[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", nargs="+", required=True)
    ap.add_argument("--test", nargs="*", default=[])
    ap.add_argument("--epochs", type=int, default=40)
    ap.add_argument("--name", default="model")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--mix", type=float, default=0.0)
    a = ap.parse_args()

    torch.manual_seed(a.seed)
    rng = np.random.default_rng(a.seed)
    t0 = time.time()
    pool = [w for s in a.train for w in data.cache_set(s) if w["status"] != "unknown"]
    order = rng.permutation(len(pool))
    n_val = max(2, len(pool) // 5)
    val = [pool[i] for i in order[:n_val]]
    trn = [pool[i] for i in order[n_val:]]
    print("training windows %d (%d shots), held back %d (%d shots); loaded in %.0f s" %
          (len(trn), sum(len(w["marks"]) for w in trn), len(val), sum(len(w["marks"]) for w in val), time.time() - t0))

    model = ShotNet(pool[0]["S"].shape[0])
    print("network parameters: %d" % count_parameters(model))
    opt = torch.optim.AdamW(model.parameters(), lr=2e-3, weight_decay=1e-3)
    sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=2e-3, total_steps=a.epochs * int(np.ceil(len(trn) / 8)))
    loss_fn = nn.BCEWithLogitsLoss(pos_weight=torch.tensor(3.0))
    t0 = time.time()
    for ep in range(a.epochs):
        model.train()
        tot = 0.0
        for X, Y in batches(trn, 8, 1500, rng, a.mix):
            opt.zero_grad()
            loss = loss_fn(model(X), Y)
            loss.backward()
            opt.step()
            sched.step()
            tot += loss.item()
        if (ep + 1) % 10 == 0 or ep == 0:
            print(flush=True, end="")
            print("epoch %3d  loss %.4f  (%.0f s)" % (ep + 1, tot / max(1, int(np.ceil(len(trn) / 8))), time.time() - t0))

    thr = best_threshold(val, predict(model, val))
    print("threshold chosen on held-back training windows: %.2f" % thr)
    os.makedirs(os.path.join(HERE, "models"), exist_ok=True)
    torch.save(dict(state=model.state_dict(), n_freq=pool[0]["S"].shape[0], threshold=thr),
               os.path.join(HERE, "models", a.name + ".pt"))

    for s in a.test:
        ws = data.cache_set(s)
        probs = predict(model, ws)
        dets = {w["window"]: score.pick(probs[w["window"]], thr) for w in ws}
        tot, rows = score.score_set(ws, dets)
        print("RESULT %-5s %s" % (s, score.describe(tot)))
        print("DETAIL " + " | ".join("%s %d/%d s%d" % (r["window"], r["found"], r["shots"], r["spare"]) for r in rows))
        with open(os.path.join(HERE, "models", "%s_detections_%s.csv" % (a.name, s)), "w", newline="") as f:
            c = csv.writer(f)
            c.writerow(["window", "start_s", "end_s"])
            for w in ws:
                c.writerows([[w["window"], "%.2f" % d0, "%.2f" % d1] for d0, d1 in dets[w["window"]]])


if __name__ == "__main__":
    main()
