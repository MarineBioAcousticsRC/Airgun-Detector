"""Choose the next windows to label, using the network's own output.

    python select_windows.py train2

Steps:
  1. Take a random sample of five-minute blocks from the site scans, on days
     not used by any labelled set (set 1, set 2, train1).
  2. Run the trained networks (all models/train1_seed*.pt, averaged) on each.
  3. Pick about 100 where a label should teach the network the most:
       - it fires confidently, but the scan saw no regular pulsing
       - it is unsure (many mid-probability peaks)
       - the separately trained copies of the network disagree
       - the scan saw regular pulsing, but the network finds little
  4. Write <set>/windows.csv and the network's detections on those windows.
The original detector and the prototype are then run on the chosen windows
(score_detector.m) and the three are combined into proposals (make_proposals.py).
"""
import csv
import datetime
import glob
import os
import sys
import time

import numpy as np
import scipy.io as sio
import torch

import data
import score
from model import ShotNet

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = data.TESTS
OUT = sys.argv[1] if len(sys.argv) > 1 else "train2"
PROPOSAL_THR = 0.35 if OUT == "train2" else 0.30

if OUT == "train2":      # second round: chosen with the networks trained on train1
    SEED, MODELS, PREFIX = 20261119, "train1_seed*.pt", "U"
    USED = ("set1", "set2", "train1")
    UNSURE = (0.2, 0.6)
    #        site     scan file              blocks to try   minimum background level (dB counts)
    SITES = [("GC_15", "scan_GOM_GC_15.mat", 400, 20.0),   # GC_15 level collapses after mid-Oct 2022; those windows were dropped in train1
             ("MC_17", "scan_GOM_MC_17.mat", 400, None),
             ("DC_15", "scan_GOM_DC_15.mat", 400, None),
             ("CE_03", "scan_ce03.mat", 150, None),
             ("Y2B_01", "scan_GOM_Y2B_01.mat", 150, None),
             ("SL_03", "scan_GOM_SL_03.mat", 150, None),
             ("DT_15", "scan_GOM_DT_15.mat", 150, None),
             ("EB_01", "scan_GOM_EB_01.mat", 100, None)]
    #           reason                                                  how many
    QUOTA = {"network confident, scan saw no regular pulsing": 30,
             "network unsure": 30,
             "network copies disagree": 20,
             "scan saw regular pulsing, network finds little": 20}
    # site groups and the share of each reason's windows they get
    GROUPS = [({"GC_15", "MC_17", "DC_15"}, 2 / 3), ({"CE_03", "Y2B_01", "SL_03", "DT_15", "EB_01"}, 1 / 3)]
else:                    # third round: chosen with the networks trained on train1 + train2 (with and without
    #                      background mixing). Aimed at what is still missed: faint shots and fast shooting in
    #                      loud background, above all at DC_15. Days of the frozen test set 3 are not touched.
    SEED, MODELS, PREFIX = 20261120, "train12*seed*.pt", "V"
    USED = ("set1", "set2", "set3", "train1", "train2")
    UNSURE = (0.1, 0.5)
    SITES = [("DC_15", "scan_GOM_DC_15.mat", 700, None),
             ("MC_17", "scan_GOM_MC_17.mat", 350, None),
             ("GC_15", "scan_GOM_GC_15.mat", 250, 20.0),
             ("CE_03", "scan_ce03.mat", 60, None),
             ("Y2B_01", "scan_GOM_Y2B_01.mat", 60, None),
             ("SL_03", "scan_GOM_SL_03.mat", 60, None),
             ("DT_15", "scan_GOM_DT_15.mat", 60, None),
             ("EB_01", "scan_GOM_EB_01.mat", 60, None)]
    QUOTA = {"network sees pulses but is not confident": 40,
             "network unsure": 20,
             "network copies disagree": 20,
             "scan saw regular pulsing, network finds little": 20}
    GROUPS = [({"DC_15"}, 0.5), ({"GC_15", "MC_17"}, 0.3), ({"CE_03", "Y2B_01", "SL_03", "DT_15", "EB_01"}, 0.2)]
rng = np.random.default_rng(SEED)


def datenum_to_dt(dn):
    return datetime.datetime(1, 1, 1) + datetime.timedelta(days=float(dn) - 367.0)


def used_days():
    used = set()
    for folder in USED:
        for r in csv.DictReader(open(os.path.join(TESTS, folder, "windows.csv"))):
            used.add((r["site"], r["start_time"][:10]))
    return used


def load_models():
    models = []
    for p in sorted(glob.glob(os.path.join(HERE, "models", MODELS))):
        ck = torch.load(p)
        m = ShotNet(ck["n_freq"])
        m.load_state_dict(ck["state"])
        m.eval()
        models.append(m)
    return models


def main():
    models = load_models()
    print("networks averaged: %d" % len(models))
    used = used_days()
    cands = []
    t0 = time.time()
    for site, scan, n_try, min_bg in SITES:
        m = sio.loadmat(os.path.join(TESTS, "scans", scan), squeeze_me=True, struct_as_record=False)
        R, fl = m["R"], m["fl"]
        scan_score = np.minimum(R[:, 3], R[:, 5]) if R.shape[1] == 7 else R[:, 3]
        ok = R[:, -1] == data.SEG
        if min_bg is not None:
            ok &= R[:, 6] >= min_bg
        days = [datenum_to_dt(d).strftime("%Y-%m-%d") for d in R[:, 2]]
        ok &= np.array([(site, d) not in used for d in days])
        idx = rng.permutation(np.where(ok)[0])[:n_try]
        paths = data.file_index({site})
        for i in idx:
            name = fl[int(R[i, 0]) - 1].name
            w = dict(site=site, file=name, path=paths[name], seg=int(R[i, 1]),
                     start=datenum_to_dt(R[i, 2]).strftime("%Y-%m-%d %H:%M:%S"), scan=float(scan_score[i]))
            try:
                S = data.spectrogram(data.read_window(w))
            except Exception as err:                      # unreadable stretch: skip it
                print("  skipped %s seg %d: %s" % (name, w["seg"], err))
                continue
            with torch.no_grad():
                P = np.stack([torch.sigmoid(mm(torch.from_numpy(S).unsqueeze(0)))[0].numpy() for mm in models])
            p = P.mean(0)
            pk = (score.pick(p, 0.2)[:, 0] / data.HOP_S).round().astype(int)
            w["n_conf"] = int((p[pk] >= 0.6).sum())
            pk = (score.pick(p, UNSURE[0])[:, 0] / data.HOP_S).round().astype(int)
            w["n_unsure"] = int(((p[pk] >= UNSURE[0]) & (p[pk] < UNSURE[1])).sum())
            w["disagree"] = float(P[:, pk].std(0).mean()) if len(pk) else 0.0
            w["det"] = score.pick(p, PROPOSAL_THR)
            cands.append(w)
        print("%s: %d blocks tried (%.0f s so far)" % (site, len(idx), time.time() - t0))

    # choose windows: each reason takes its best-ranked blocks, at most 2 per site-day,
    # shared between the site groups as set in GROUPS
    rank = {"network sees pulses but is not confident":       # many low-probability peaks, few confident ones
                lambda w: (w["n_unsure"] - w["n_conf"]) if (w["n_unsure"] >= 8 and w["n_conf"] <= 0.5 * w["n_unsure"]) else -1,
            "network confident, scan saw no regular pulsing":
                lambda w: w["n_conf"] if (w["scan"] < 0.3 and w["n_conf"] >= 4) else -1,
            "network unsure": lambda w: w["n_unsure"],
            "network copies disagree": lambda w: w["disagree"] if (w["n_conf"] + w["n_unsure"]) >= 3 else -1,
            "scan saw regular pulsing, network finds little":
                lambda w: (w["scan"] - 0.02 * w["n_conf"]) if (w["scan"] >= 0.5 and w["n_conf"] <= 5) else -1}
    chosen, per_day = [], {}
    for reason, quota in QUOTA.items():
        pool = sorted([w for w in cands if rank[reason](w) > -1 and "reason" not in w], key=rank[reason], reverse=True)
        want = [int(round(quota * share)) for _, share in GROUPS]
        want[0] += quota - sum(want)
        got = [0] * len(GROUPS)
        per_site = {}
        for w in pool:
            grp = [k for k, (sites, _) in enumerate(GROUPS) if w["site"] in sites][0]
            cap = int(np.ceil(want[grp] / len(GROUPS[grp][0]))) + 1     # spread over the sites of the group
            key = (w["site"], w["start"][:10])
            if got[grp] >= want[grp] or per_site.get(w["site"], 0) >= cap or per_day.get(key, 0) >= 2:
                continue
            w["reason"] = reason
            chosen.append(w)
            got[grp] += 1
            per_site[w["site"]] = per_site.get(w["site"], 0) + 1
            per_day[key] = per_day.get(key, 0) + 1
        print("%-48s chosen %d (by site group: %s)" % (reason, sum(got), got))

    out = os.path.join(TESTS, OUT)
    os.makedirs(out, exist_ok=True)
    chosen.sort(key=lambda w: (w["site"], w["start"]))
    with open(os.path.join(out, "windows.csv"), "w", newline="") as f:
        c = csv.writer(f)
        c.writerow(["window", "condition_guess", "site", "file", "start_segment", "start_time", "proposed_shots"])
        for k, w in enumerate(chosen):
            w["window"] = "%s%03d" % (PREFIX, k + 1)
            c.writerow([w["window"], w["reason"].replace(",", ";"), w["site"], w["file"], w["seg"], w["start"], 0])
    with open(os.path.join(out, "detections_proposals_cnn.csv"), "w", newline="") as f:
        c = csv.writer(f)
        c.writerow(["window", "start_s", "end_s"])
        for w in chosen:
            c.writerows([[w["window"], "%.2f" % a, "%.2f" % b] for a, b in w["det"]])
    with open(os.path.join(out, "truth.csv"), "w", newline="") as f:
        f.write("window,time_s,source\n")
    with open(os.path.join(out, "candidates.csv"), "w", newline="") as f:
        c = csv.writer(f)
        c.writerow(["site", "file", "start_segment", "start_time", "scan_score", "confident_peaks", "unsure_peaks",
                    "disagreement", "chosen_as"])
        for w in cands:
            c.writerow([w["site"], w["file"], w["seg"], w["start"], "%.2f" % w["scan"], w["n_conf"], w["n_unsure"],
                        "%.3f" % w["disagree"], w.get("window", "")])
    by_site = {}
    for w in chosen:
        by_site[w["site"]] = by_site.get(w["site"], 0) + 1
    print("%d windows chosen from %d blocks: %s" % (len(chosen), len(cands), by_site))


if __name__ == "__main__":
    main()
