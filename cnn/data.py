"""Data handling for the spectrogram (CNN) airgun detector.

A "set" is a folder with windows.csv (one five-minute window per row) and,
once reviewed, review.csv (what each window is), truth.csv (the shot marks)
and optionally humps.csv. The sets are folders in ../evaluation (set1, set2, set3,
train1, train2, ...).
"""
import csv
import glob
import os

import numpy as np
import soundfile as sf
from scipy import signal

FS = 2000                 # sample rate of the df100 data (Hz)
SEG = 150000              # samples in one 75 s segment
WIN_SEGS = 4              # a window is four segments = 300 s
HOP_S = 0.1               # time step of the spectrogram and of the network output (s)
NPERSEG = 800             # 0.4 s analysis window -> 2.5 Hz resolution
FMIN, FMAX = 5.0, 250.0   # frequency range given to the network (Hz)
ROOT = "//frosty.ucsd.edu/GofMX_Decimated_8"

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.join(os.path.dirname(HERE), "evaluation")     # where the labelled sets live
INDEX = os.path.join(HERE, "file_index.csv")


def set_dir(name):
    """Folder of a set: '' or 'set1' is the first test set."""
    return os.path.join(TESTS, name or "set1")


def file_index(sites):
    """Map recording file name -> full path, cached in file_index.csv."""
    idx = {}
    if os.path.exists(INDEX):
        idx = {r["file"]: r["path"] for r in csv.DictReader(open(INDEX))}
    have = {os.path.basename(p).split("_")[1] + "_" + os.path.basename(p).split("_")[2] for p in idx.values()}
    for site in sorted(set(sites) - have):
        for p in glob.glob("%s/GOM_%s/*_df100/*df100.x.wav" % (ROOT, site)):
            idx[os.path.basename(p)] = p.replace(os.sep, "/")
    with open(INDEX, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["file", "path"])
        w.writerows(sorted(idx.items()))
    return idx


def load_set(name):
    """Windows of a set with their review status and shot marks.

    Returns a list of dicts: window, site, path, seg, status, marks (seconds),
    hump_join (seconds, or None). Windows that were dropped or never reviewed
    are left out. "unknown" windows are returned with status "unknown".
    """
    d = set_dir(name)
    rows = list(csv.DictReader(open(os.path.join(d, "windows.csv"))))
    idx = file_index({r["site"] for r in rows})
    status = {}
    rv = os.path.join(d, "review.csv")
    if os.path.exists(rv):
        status = {r["window"]: r["status"] for r in csv.DictReader(open(rv))}
    marks = {}
    tr = os.path.join(d, "truth.csv")
    if os.path.exists(tr):
        for r in csv.DictReader(open(tr)):
            marks.setdefault(r["window"], []).append(float(r["time_s"]))
    humps = {}
    hp = os.path.join(d, "humps.csv")
    if os.path.exists(hp):
        humps = {r["window"]: float(r["join_s"]) for r in csv.DictReader(open(hp))}
    out = []
    for r in rows:
        st = status.get(r["window"])
        if st is None:      # a set reviewed before review.csv existed
            st = "no airguns" if "no airguns" in r["condition_guess"] else "airguns"
        if st in ("dropped", "not reviewed", "airguns unclear"):   # unclear: airguns present, shots not markable
            continue
        m = sorted(marks.get(r["window"], [])) if st == "airguns" else []
        if st == "airguns" and not m:
            st = "no airguns"
        out.append(dict(window=r["window"], site=r["site"], path=idx[r["file"]], seg=int(r["start_segment"]),
                        status=st, marks=np.array(m), hump_join=humps.get(r["window"])))
    return out


def read_audio(path, first_sample, n_samples):
    """Read raw counts from an xwav (first_sample is 0-based)."""
    y, fs = sf.read(path, start=first_sample, frames=n_samples, dtype="int16")
    assert fs == FS, "unexpected sample rate %d" % fs
    return y.astype(np.float32)


def read_window(w):
    return read_audio(w["path"], (w["seg"] - 1) * SEG, WIN_SEGS * SEG)


def spectrogram(y):
    """Normalised log spectrogram, shape (n_freq, n_time), one column every HOP_S.

    Each frequency row has its own median level (over the stretch given)
    removed, which takes out the site's background spectrum and steady tones,
    so the network sees how far the sound stands above the background.
    """
    y = y - y.mean()
    f, t, Z = signal.stft(y, fs=FS, window="hann", nperseg=NPERSEG, noverlap=NPERSEG - int(HOP_S * FS),
                          boundary="even", padded=True)
    keep = (f >= FMIN) & (f <= FMAX)
    S = 10.0 * np.log10(np.abs(Z[keep]) ** 2 + 1e-6)
    S = S - np.median(S, axis=1, keepdims=True)
    S = np.clip(S, -10.0, 40.0) / 20.0
    return S.astype(np.float32)


def target(marks, n_time, sigma_s=0.4):
    """Per-time-step target: a bump of height 1 centred on each shot mark."""
    t = np.arange(n_time) * HOP_S
    y = np.zeros(n_time, dtype=np.float32)
    for m in marks:
        y = np.maximum(y, np.exp(-0.5 * ((t - m) / sigma_s) ** 2))
    return y


def cache_set(name):
    """Spectrograms and targets for every usable window of a set (cached on disk)."""
    cdir = os.path.join(HERE, "cache")
    os.makedirs(cdir, exist_ok=True)
    out = []
    for w in load_set(name):
        p = os.path.join(cdir, "%s_%s.npy" % (name or "set1", w["window"]))
        if os.path.exists(p):
            S = np.load(p)
        else:
            S = spectrogram(read_window(w))
            np.save(p, S)
        w = dict(w)
        w["S"] = S
        w["y"] = target(w["marks"], S.shape[1])
        out.append(w)
    return out
