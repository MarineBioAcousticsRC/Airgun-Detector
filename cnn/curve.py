"""How a trained network's results change with its detection threshold.

    python curve.py train1_seed1 set1 set2

For description only: the threshold the model uses is the one chosen on
held-back training windows and stored with the model.
"""
import os
import sys

import numpy as np
import torch

import data
import score
from model import ShotNet

HERE = os.path.dirname(os.path.abspath(__file__))
name, sets = sys.argv[1], sys.argv[2:]
ck = torch.load(os.path.join(HERE, "models", name + ".pt"))
model = ShotNet(ck["n_freq"])
model.load_state_dict(ck["state"])
model.eval()
print("model %s, stored threshold %.2f" % (name, ck["threshold"]))
for s in sets:
    ws = data.cache_set(s)
    with torch.no_grad():
        probs = {w["window"]: torch.sigmoid(model(torch.from_numpy(w["S"]).unsqueeze(0)))[0].numpy() for w in ws}
    print(s)
    for thr in (0.2, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9):
        tot, _ = score.score_set(ws, {w["window"]: score.pick(probs[w["window"]], thr) for w in ws})
        print("  threshold %.2f: %s" % (thr, score.describe(tot)))
