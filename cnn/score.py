"""Scoring with the same rules as score_detector.m.

found   an event (shot) with a detection on it (within TOL seconds)
merged  an event whose only detections already hold another event
extra   a detection that holds no event, in a window with airguns
false   a detection in a window without airguns
In a hump window, marks closer than hump_join form one event, a detection
anywhere in the hump counts, and every detection in the hump belongs to it.
"""
import numpy as np
from scipy.signal import find_peaks

from data import HOP_S

TOL = 2.0


def pick(prob, threshold, min_sep_s=2.5):
    """Detections from a probability curve: peaks at least min_sep_s apart."""
    pk, _ = find_peaks(prob, height=threshold, distance=max(1, int(round(min_sep_s / HOP_S))))
    t = pk * HOP_S
    return np.stack([t, t + 0.5], axis=1) if len(t) else np.zeros((0, 2))


def score_window(marks, det, hump_join=None):
    marks = np.sort(np.asarray(marks, dtype=float))
    events = [[m] for m in marks]
    if hump_join and len(marks):
        events = [[marks[0]]]
        for m in marks[1:]:
            if m - events[-1][-1] > hump_join:
                events.append([m])
            else:
                events[-1].append(m)
    used = np.zeros(len(det), dtype=bool)
    found = merged = 0
    for ev in events:
        if len(ev) == 1:
            k = np.where((ev[0] >= det[:, 0] - TOL) & (ev[0] <= det[:, 1] + TOL))[0]
        else:
            k = np.where((det[:, 1] >= ev[0] - TOL) & (det[:, 0] <= ev[-1] + TOL))[0]
        if len(k) == 0:
            continue
        free = k[~used[k]]
        if len(free) == 0:
            merged += 1
            continue
        found += 1
        if len(ev) == 1:
            used[free[0]] = True
        else:
            used[k] = True
    return dict(shots=len(events), detections=len(det), found=found, merged=merged, spare=int((~used).sum()))


def score_set(windows, dets):
    """windows: from data.load_set; dets: {window: array of [start, end]}. Returns totals and per-window rows."""
    rows = []
    for w in windows:
        r = score_window(w["marks"], dets[w["window"]], w["hump_join"])
        r.update(window=w["window"], site=w["site"], status=w["status"])
        rows.append(r)
    air = [r for r in rows if r["shots"] > 0]
    quiet = [r for r in rows if r["shots"] == 0]
    tot = dict(shots=sum(r["shots"] for r in air), found=sum(r["found"] for r in air),
               merged=sum(r["merged"] for r in air), extra=sum(r["spare"] for r in air),
               false=sum(r["spare"] for r in quiet), n_air=len(air), n_quiet=len(quiet))
    return tot, rows


def describe(tot):
    return ("found %d of %d (%.0f%%), merged %d, extra %d in %d airgun windows; false %d in %d other windows" %
            (tot["found"], tot["shots"], 100.0 * tot["found"] / max(tot["shots"], 1), tot["merged"], tot["extra"],
             tot["n_air"], tot["false"], tot["n_quiet"]))
