"""Run the airgun network over whole deployments.

    python run_network.py GOM_MC_17 --out E:/AirgunDetections
    python run_network.py GOM_GC_15 GOM_MC_17 GOM_DC_15 --out E:/AirgunDetections
    python run_network.py D:/somewhere/GOM_MC_17 --out E:/AirgunDetections

Each deployment is a folder name under --root (default the lab server,
//frosty.ucsd.edu/GofMX_Decimated_8) or a full path. Every *df100.x.wav in
its disk folders (or directly in the folder) is processed, as is every
*.d100.x.wav (the naming of the older GofMX deployments, whose df100 files sit
in folders named *_df20).

Output has the same layout as the original detector (xcorr_airgun_df100.m):
a folder <out>/<deployment>_AirgunDetections with one .mat per xwav, named
like the xwav, holding the same variables:
  allSmpPts   [start end] sample of each detection (1-based, in the xwav)
  allExp      [start end] time as MATLAB datenum
  bt          [start-5000, end+5000, 0, start datenum, end datenum]
  rawStart, rawDur, parm
  allCorrVal, allDur, allRmsNBefore, allRmsNAfter, allRmsDet, allPpNBefore,
  allPpNAfter, allPpDet      NaN: the network does not measure these
The network gives a start time only, so each detection's end is set equal to
its start. One extra variable, allProb, holds the network's probability.
Give --out a folder of its own, so results of the original detector are not
overwritten.

Three extra files are written next to the .mat files:
  <deployment>_detections.csv      all detections of the deployment
  <deployment>_minutes_per_day.csv per day: minutes recorded and minutes with
                                   at least one detection
  run_info.txt                     model, threshold and settings used

Files that already have a result are skipped, so an interrupted run can simply
be started again. Use --overwrite to redo them.

How it works: each xwav is split into stretches of continuous recording. A
stretch is read in 300 s pieces that overlap by 60 s (the network was trained
on 300 s pieces, and the level of each frequency is measured over the piece).
The probabilities of the pieces are joined, leaving out 30 s at the inner
edges of each piece, and peaks above the threshold at least 2.5 s apart are
the detections. The network gives a start time for each shot, not a duration.

Data with instrument noise should be removed from the results afterwards (or
not run): the network was not trained to reject every kind of it.
"""
import argparse
import csv
import datetime
import glob
import os
import struct
import sys
import time

import numpy as np
import scipy.io as sio
import soundfile as sf
import torch

import model

ROOT = "//frosty.ucsd.edu/GofMX_Decimated_8"
PIECE_S, MARGIN_S = 300.0, 30.0          # length of a piece and what is dropped at its inner edges
HOP = int(round(model.HOP_S * model.FS))   # samples per output step (200)
PIECE = int(PIECE_S * model.FS)
MARGIN = int(round(MARGIN_S / model.HOP_S))        # in output steps
STEP = PIECE - 2 * MARGIN * HOP                   # samples between piece starts
MIN_RUN_S = 10.0                         # stretches shorter than this are skipped
EPOCH = datetime.datetime(1, 1, 1)


def read_xwav_header(path):
    """Start time, first sample and length of every raw file in an xwav (same layout as readxwavhd.m)."""
    with open(path, "rb") as f:
        head = f.read(12)
        assert head[:4] == b"RIFF" and head[8:12] == b"WAVE", "not a wav file: " + path
        cid, size = struct.unpack("<4sI", f.read(8))
        assert cid == b"fmt ", "unexpected header in " + path
        fmt = f.read(size)
        n_ch, fs = struct.unpack("<HI", fmt[2:8])
        bytes_per = struct.unpack("<H", fmt[14:16])[0] // 8
        cid, size = struct.unpack("<4sI", f.read(8))
        assert cid == b"harp", "no HARP header (not an xwav): " + path
        h = f.read(size)
    n_raw = struct.unpack("<H", h[36:38])[0]
    raws = []
    for i in range(n_raw):
        r = h[56 + 32 * i: 56 + 32 * (i + 1)]
        yr, mo, dy, hr, mi, sec, ticks, byte_loc, byte_len = struct.unpack("<6BHII", r[:16])
        t = datetime.datetime(2000 + yr, mo, dy, hr, mi, sec) + datetime.timedelta(milliseconds=ticks)
        raws.append((t, byte_loc, byte_len))
    loc0 = raws[0][1]
    frame = bytes_per * n_ch
    return fs, n_ch, [dict(time=t, first=(loc - loc0) // frame, n=ln // frame) for t, loc, ln in raws]


def continuous_runs(raws, fs):
    """Group raw files into stretches of continuous recording: (start time, first sample, samples)."""
    runs = []
    for r in raws:
        if runs:
            t0, first, n = runs[-1]
            gap = (r["time"] - (t0 + datetime.timedelta(seconds=n / fs))).total_seconds()
            if abs(gap) < 0.05 and r["first"] == first + n:
                runs[-1] = (t0, first, n + r["n"])
                continue
        runs.append((r["time"], r["first"], r["n"]))
    return runs


def datenum(t):
    """MATLAB datenum of a datetime."""
    return (t - EPOCH).total_seconds() / 86400.0 + 367.0


def run_probability(snd, first, n, models, device, batch):
    """Shot-start probability for one continuous stretch, one value every HOP samples."""
    n = (n // HOP) * HOP
    n_out = n // HOP + 1
    prob = np.zeros(n_out, dtype=np.float32)
    if n <= PIECE:
        starts = [0]
    else:
        starts = list(range(0, n - PIECE, STEP)) + [n - PIECE]
    pieces = []          # (start sample, first output step kept, last + 1)
    for k, s0 in enumerate(starts):
        o0 = s0 // HOP
        o1 = o0 + min(PIECE, n) // HOP + 1
        a = 0 if k == 0 else pieces[-1][2]             # start exactly where the piece before stopped
        b = o1 if k == len(starts) - 1 else o1 - MARGIN
        assert k == 0 or a - o0 >= MARGIN              # so every kept value is at least 30 s from a piece edge
        pieces.append((s0, a, b))
    for i in range(0, len(pieces), batch):
        group = pieces[i:i + batch]
        X = []
        for s0, _, _ in group:
            snd.seek(first + s0)
            y = snd.read(min(PIECE, n), dtype="int16", always_2d=True)[:, 0].astype(np.float32)
            X.append(model.spectrogram(y))
        with torch.no_grad():
            x = torch.from_numpy(np.stack(X)).to(device)
            p = torch.stack([torch.sigmoid(m(x)) for m in models]).mean(0).cpu().numpy()
        for (s0, a, b), pk in zip(group, p):
            o0 = s0 // HOP
            prob[a:b] = pk[a - o0:b - o0]
    return prob


def process_file(path, out_mat, models, thr, device, batch, model_name, segments=None):
    """Detections of one xwav, saved to out_mat in the original detector's layout."""
    fs, n_ch, raws = read_xwav_header(path)
    assert fs == model.FS, "sample rate %d, expected %d: %s" % (fs, model.FS, path)
    all_raws = raws
    if segments:
        raws = raws[segments[0] - 1: segments[0] - 1 + segments[1]]
    runs = continuous_runs(raws, fs)
    smp, dn, pr = [], [], []
    with sf.SoundFile(path) as snd:
        for t0, first, n in runs:
            if n < MIN_RUN_S * fs:
                continue
            prob = run_probability(snd, first, n, models, device, batch)
            for ts in model.pick(prob, thr):
                k = int(round(ts / model.HOP_S))
                smp.append(first + k * HOP + 1)
                dn.append(datenum(t0 + datetime.timedelta(seconds=float(ts))))
                pr.append(float(prob[k]))
    n = len(smp)
    col = lambda v: np.asarray(v, dtype=float).reshape(-1, 1)
    blank = np.full((n, 1), np.nan) if n else np.zeros((0, 0))
    two = lambda v: np.hstack([col(v), col(v)]) if n else np.zeros((0, 0))      # end = start
    out = dict(allSmpPts=two(smp), allExp=two(dn), allProb=col(pr) if n else np.zeros((0, 0)),
               allCorrVal=blank, allDur=blank, allRmsNBefore=blank, allRmsNAfter=blank, allRmsDet=blank,
               allPpNBefore=blank, allPpNAfter=blank, allPpDet=blank,
               rawStart=np.array([[r["time"].year, r["time"].month, r["time"].day, r["time"].hour, r["time"].minute,
                                   r["time"].second + r["time"].microsecond / 1e6] for r in all_raws], dtype=float),
               rawDur=np.array([[r["n"] / fs for r in all_raws]], dtype=float),
               parm=dict(detector="spectrogram network (cnn/run_network.py)", model=model_name, threshold=float(thr),
                         minSep_s=2.5, piece_s=PIECE_S, edge_s=MARGIN_S, fmin_Hz=model.FMIN, fmax_Hz=model.FMAX))
    out["bt"] = (np.hstack([col(smp) - 5000, col(smp) + 5000, np.zeros((n, 1)), col(dn), col(dn)])
                 if n else np.zeros((0, 0)))
    tmp = out_mat + ".part"
    with open(tmp, "wb") as f:
        sio.savemat(f, out, oned_as="column")
    os.replace(tmp, out_mat)                         # a result file only exists once it is complete
    return n, [(t0, k / fs) for t0, _, k in runs]


def minutes_per_day(folder, name, effort):
    """Combine the per-file results into one detections file and a per-day summary."""
    det_min, rows = {}, []
    for p in sorted(glob.glob(os.path.join(folder, "*.mat"))):
        m = sio.loadmat(p)
        if m["allExp"].size == 0:
            continue
        for dn, smp, pr in zip(m["allExp"][:, 0], m["allSmpPts"][:, 0], m["allProb"][:, 0]):
            t = (EPOCH + datetime.timedelta(days=float(dn) - 367.0)).strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]
            rows.append(dict(start_time_utc=t, start_datenum="%.8f" % dn, start_sample=int(smp),
                             probability="%.3f" % pr, file=os.path.basename(p)[:-4] + ".x.wav"))
            det_min.setdefault(t[:10], set()).add(t[:16])
    rows.sort(key=lambda r: r["start_time_utc"])
    with open(os.path.join(folder, name + "_detections.csv"), "w", newline="") as f:
        c = csv.DictWriter(f, fieldnames=["start_time_utc", "start_datenum", "start_sample", "probability", "file"])
        c.writeheader()
        c.writerows(rows)
    rec_min = {}
    for t0, dur in effort:                           # minutes that have any recording in them
        t = t0.replace(second=0, microsecond=0)
        end = t0 + datetime.timedelta(seconds=dur)
        while t < end:
            rec_min.setdefault(t.strftime("%Y-%m-%d"), set()).add(t)
            t += datetime.timedelta(minutes=1)
    with open(os.path.join(folder, name + "_minutes_per_day.csv"), "w", newline="") as f:
        c = csv.writer(f)
        c.writerow(["date", "minutes_recorded", "minutes_with_airguns", "detections"])
        counts = {}
        for r in rows:
            counts[r["start_time_utc"][:10]] = counts.get(r["start_time_utc"][:10], 0) + 1
        for day in sorted(set(rec_min) | set(det_min)):
            c.writerow([day, len(rec_min.get(day, ())), len(det_min.get(day, ())), counts.get(day, 0)])
    return len(rows)


def main():
    ap = argparse.ArgumentParser(description="Run the airgun network over whole deployments.")
    ap.add_argument("deployments", nargs="+", help="deployment folder names under --root, or full paths")
    ap.add_argument("--out", required=True, help="where results go")
    ap.add_argument("--root", default=ROOT)
    ap.add_argument("--model", default="network_airgun")
    ap.add_argument("--threshold", type=float, default=None, help="default: the one stored with the model")
    ap.add_argument("--mask", nargs="+", default=["*df100.x.wav", "*.d100.x.wav"], help="which files to process")
    ap.add_argument("--batch", type=int, default=8)
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--segments", type=int, nargs=2, metavar=("FIRST", "COUNT"),
                    help="for testing: only these raw files (75 s segments) of each xwav")
    a = ap.parse_args()

    models, stored_thr = model.load(a.model)
    assert models, "model not found: " + a.model
    thr = a.threshold if a.threshold is not None else stored_thr
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    models = [m.to(device) for m in models]
    print("model %s (%d copies), threshold %.2f, running on %s" % (a.model, len(models), thr, device))

    for dep in a.deployments:
        base = dep if os.path.isdir(dep) else os.path.join(a.root, dep)
        name = os.path.basename(os.path.normpath(base))
        files = sorted({f for m in a.mask for f in glob.glob(os.path.join(base, "*", m)) + glob.glob(os.path.join(base, m))})
        if not files:
            print("%s: no files matching %s found under %s" % (name, " or ".join(a.mask), base))
            continue
        folder = os.path.join(a.out, name + "_AirgunDetections")
        os.makedirs(folder, exist_ok=True)
        with open(os.path.join(folder, "run_info.txt"), "w") as f:
            f.write("model: %s (%d copies)\nthreshold: %.2f\npiece length: %.0f s, edge left out: %.0f s\n"
                    "minimum spacing of detections: 2.5 s\ncommand: %s\nstarted: %s\n" %
                    (a.model, len(models), thr, PIECE_S, MARGIN_S, " ".join(sys.argv),
                     datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))
        print("%s: %d files -> %s" % (name, len(files), folder))
        effort, t0, done_s = [], time.time(), 0.0
        for i, path in enumerate(files):
            out_mat = os.path.join(folder, os.path.basename(path).split(".df100")[0].split(".x.wav")[0] + ".mat")
            if os.path.exists(out_mat) and not a.overwrite:
                fs, _, raws = read_xwav_header(path)
                effort += [(t, n / fs) for t, _, n in continuous_runs(raws, fs)]
                print("  %3d/%d %s: already done" % (i + 1, len(files), os.path.basename(path)))
                continue
            try:
                n_det, runs = process_file(path, out_mat, models, thr, device, a.batch, a.model, a.segments)
            except Exception as err:                 # one unreadable file should not stop the deployment
                print("  %3d/%d %s: FAILED (%s)" % (i + 1, len(files), os.path.basename(path), err))
                continue
            effort += runs
            done_s += sum(d for _, d in runs)
            el = time.time() - t0
            print("  %3d/%d %s: %d detections in %.1f h of recording (%.0f s; %.0f x real time)" %
                  (i + 1, len(files), os.path.basename(path), n_det, sum(d for _, d in runs) / 3600.0, el,
                   done_s / max(el, 1e-9)), flush=True)
        n = minutes_per_day(folder, name, effort)
        print("%s: %d detections in all; summary written to %s" % (name, n, folder))


if __name__ == "__main__":
    main()
