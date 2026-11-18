"""Build Airgun_detector_report.pdf from the results in this folder.

    python make_report.py

Needs reportlab. All scores in the tables are read from the result files, so
the report cannot drift from them:
  setN/scores_reviewed_original.csv, setN/scores_reviewed_prototype.csv
      written by score_both_sets.m (original and shape-based detector)
  setN/detections_network.csv
      the averaged network's detections (cnn/ensemble.py), scored here with
      the same rules (cnn/score.py)
Figures come from figures/ (plot_methods.m, plot_network.py, plot_pr.py) and
setN/plots_compare/ (plot_comparison.m).
"""
import csv
import os
import sys
from collections import defaultdict

import numpy as np
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (Image, KeepTogether, Paragraph, SimpleDocTemplate, Spacer, Table,
                                TableStyle)

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "cnn"))
import data     # noqa: E402
import score    # noqa: E402

OUT = os.path.join(HERE, "Airgun_detector_report.pdf")
FIGS = os.path.join(HERE, "figures")
TITLE = "Airgun detection: three detectors compared"
DETS = ("original", "shape", "network")

INK = colors.HexColor("#1f2933")
MUTED = colors.HexColor("#52606d")
RULE = colors.HexColor("#cbd2d9")
HEAD_BG = colors.HexColor("#e4e7eb")

ss = getSampleStyleSheet()
body = ParagraphStyle("body", parent=ss["Normal"], fontName="Helvetica", fontSize=10,
                      leading=14, textColor=INK, spaceAfter=6)
small = ParagraphStyle("small", parent=body, fontSize=8.5, leading=11.5, textColor=MUTED)
cell = ParagraphStyle("cell", parent=body, fontSize=8.5, leading=11, spaceAfter=0)
cellb = ParagraphStyle("cellb", parent=cell, fontName="Helvetica-Bold")
h1 = ParagraphStyle("h1", parent=body, fontName="Helvetica-Bold", fontSize=14, leading=18,
                    spaceBefore=12, spaceAfter=6)
h2 = ParagraphStyle("h2", parent=body, fontName="Helvetica-Bold", fontSize=11, leading=14,
                    spaceBefore=8, spaceAfter=4)
title = ParagraphStyle("title", parent=body, fontName="Helvetica-Bold", fontSize=19, leading=24,
                       spaceAfter=4)
bullet = ParagraphStyle("bullet", parent=body, leftIndent=14, bulletIndent=3, spaceAfter=3)
caption = ParagraphStyle("caption", parent=small, spaceBefore=3, spaceAfter=10)


def P(text, style=body):
    return Paragraph(text, style)


def bullets(items):
    return [Paragraph(t, bullet, bulletText="•") for t in items]


def table(rows, widths):
    rows_p = [[Paragraph(str(c), cellb if r == 0 else cell) for c in row] for r, row in enumerate(rows)]
    t = Table(rows_p, colWidths=[w * inch for w in widths], repeatRows=1)
    t.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LINEBELOW", (0, 0), (-1, -1), 0.4, RULE),
        ("TOPPADDING", (0, 0), (-1, -1), 3.5), ("BOTTOMPADDING", (0, 0), (-1, -1), 3.5),
        ("LEFTPADDING", (0, 0), (-1, -1), 5), ("RIGHTPADDING", (0, 0), (-1, -1), 5),
        ("BACKGROUND", (0, 0), (-1, 0), HEAD_BG), ("LINEBELOW", (0, 0), (-1, 0), 0.8, MUTED),
    ]))
    return t


def figure(path, text, width=7.0):
    img = Image(path)
    w = width * inch
    img.drawHeight = w * img.imageHeight / img.imageWidth
    img.drawWidth = w
    return KeepTogether([img, P(text, caption)])


def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(MUTED)
    canvas.drawString(0.75 * inch, 0.5 * inch, TITLE)
    canvas.drawRightString(letter[0] - 0.75 * inch, 0.5 * inch, "Page %d" % doc.page)
    canvas.restoreState()


# ---------------------------------------------------------------- results, read from the files
def load_results(name):
    """Per-window results of the three detectors on one set: {window: {...}} in window order."""
    folder = os.path.join(HERE, name)
    res = {}
    for det, f in (("original", "scores_reviewed_original.csv"), ("shape", "scores_reviewed_prototype.csv")):
        for r in csv.DictReader(open(os.path.join(folder, f))):
            w = res.setdefault(r["window"], dict(window=r["window"], site=r["site"], status=r["status"],
                                                 shots=int(r["shots"])))
            w[det] = dict(found=int(r["found"]), merged=int(r["merged"]), spare=int(r["false"]))
    dets = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(folder, "detections_network.csv"))):
        dets[r["window"]].append((float(r["start_s"]), float(r["end_s"])))
    for w in data.load_set(name):
        r = score.score_window(w["marks"], np.array(dets[w["window"]]).reshape(-1, 2), w["hump_join"])
        assert r["shots"] == res[w["window"]]["shots"], (name, w["window"])     # both scorers agree on the events
        res[w["window"]]["network"] = dict(found=r["found"], merged=r["merged"], spare=r["spare"])
    labels = {r["window"]: r["condition_guess"] for r in csv.DictReader(open(os.path.join(folder, "windows.csv")))}
    for w in res.values():
        w["label"] = labels[w["window"]]
    return res


def totals(res, keep=lambda w: True):
    t = {}
    ws = [w for w in res.values() if keep(w)]
    for d in DETS:
        air = [w for w in ws if w["shots"] > 0]
        t[d] = dict(found=sum(w[d]["found"] for w in air), merged=sum(w[d]["merged"] for w in air),
                    extra=sum(w[d]["spare"] for w in air), false=sum(w[d]["spare"] for w in ws if w["shots"] == 0))
    t["shots"] = sum(w["shots"] for w in ws)
    t["n_air"] = sum(1 for w in ws if w["shots"] > 0)
    t["n_quiet"] = sum(1 for w in ws if w["shots"] == 0)
    return t


def pct(a, b):
    return "%d%%" % round(100.0 * a / b) if b else ""


def found(t, d):
    return "%d (%s)" % (t[d]["found"], pct(t[d]["found"], t["shots"]))


def prec(t, d):
    hit = t[d]["found"]
    return "%.2f" % (hit / (hit + t[d]["extra"] + t[d]["false"]))


def one(w, d):
    if w["shots"] == 0:
        return "F %d" % w[d]["spare"]
    return "%d" % w[d]["found"] + (", F %d" % w[d]["spare"] if w[d]["spare"] else "")


def window_table(res, notes=None):
    """One row per airgun window; windows without airguns are grouped."""
    notes = notes or {}
    rows = [["Window", "Site", "Shots", "Original", "Shape-based", "Network", "Note"]]
    for w in sorted((w for w in res.values() if w["shots"] > 0), key=lambda w: w["window"]):
        rows.append([w["window"], w["site"].replace("_", " "), w["shots"]] + [one(w, d) for d in DETS] +
                    [notes.get(w["window"], "")])
    quiet = [w for w in res.values() if w["shots"] == 0]
    for w in sorted((w for w in quiet if any(w[d]["spare"] for d in DETS)), key=lambda w: w["window"]):
        note = "not sure whether it contains airguns" if w["status"] == "unknown" else notes.get(w["window"], "no airguns")
        rows.append([w["window"], w["site"].replace("_", " "), "0"] + [one(w, d) for d in DETS] + [note])
    clean = sorted(w["window"] for w in quiet if not any(w[d]["spare"] for d in DETS))
    if clean:
        rows.append([", ".join(clean), "", "0", "F 0", "F 0", "F 0", "no airguns, no detections"])
    return table(rows, [1.55, 0.75, 0.5, 0.8, 0.9, 0.8, 1.7])


def cap(res, wid):
    w = res[wid]
    if w["shots"] == 0:
        return "detections: original %d, shape-based %d, network %d" % tuple(w[d]["spare"] for d in DETS)
    return "of %d shots: original %d, shape-based %d, network %d" % ((w["shots"],) + tuple(w[d]["found"] for d in DETS))


R1, R2, R3 = load_results("set1"), load_results("set2"), load_results("set3")
T1, T2, T3 = totals(R1), totals(R2), totals(R3)


def compare(name, wid):
    return os.path.join(HERE, name, "plots_compare", wid + ".png")


s = []

# ---------------------------------------------------------------- title and summary
s += [P(TITLE, title),
      P("Vanessa ZoBell &nbsp;|&nbsp; 18 November 2026 &nbsp;|&nbsp; "
        "Repository: MarineBioAcousticsRC/Airgun-Detector", small),
      P("The code review, testing and the two new detectors described here were done with Claude Code.", small),
      Spacer(1, 6)]

s += [P("Summary", h1)]
s += bullets([
    "<b>Three detectors are compared:</b> the original cross-correlation detector, a shape-based "
    "cross-correlation detector, and a small neural network that reads the spectrogram. All aim at one "
    "detection per airgun shot.",
    "<b>They were tested on 107 five-minute windows</b> in three sets, with every shot marked by eye "
    "(%d shots in all). Set 3 is the cleanest test: 57 random windows at GOM_GC_15, MC_17 and DC_15, "
    "labelled once and never used to tune or train anything." % (T1["shots"] + T2["shots"] + T3["shots"]),
    "<b>On set 3 the original detector finds %s of shots; the shape-based detector and the network "
    "each find %s.</b> Precision is similar for all three (%s, %s and %s)." %
    (pct(T3["original"]["found"], T3["shots"]), pct(T3["shape"]["found"], T3["shots"]),
     prec(T3, "original"), prec(T3, "shape"), prec(T3, "network")),
    "<b>The original detector misses faint shots and merges close ones.</b> Its threshold is a fixed "
    "number in raw recording units, and it finds almost nothing at the quieter site (13 of 265 shots at "
    "DC_15 in set 3). Changing its thresholds does not fix this.",
    "<b>The two new detectors are level on set 3.</b> The network is ahead on set 2 (%s against %s) and "
    "is less affected by instrument noise. The shape-based detector needs no training and runs in MATLAB "
    "like the original." % (pct(T2["network"]["found"], T2["shots"]), pct(T2["shape"]["found"], T2["shots"])),
    "<b>Both new detectors still miss shots in the hardest windows:</b> fast shooting in a loud "
    "background, and faint shots at DC_15.",
    "<b>Three small bugs in the original were fixed</b> and are in review as pull request #2. Nothing "
    "else described here has been proposed for the main repository yet.",
])

# ---------------------------------------------------------------- 1 original
s += [P("1. The original detector", h1),
      P("The detector (<font face='Courier'>xcorr_original/xcorr_airgun_df100.m</font>) works on decimated "
        "recordings (2 kHz sample rate). Each recording file is made of 75-second segments, and each segment "
        "is processed separately:")]
s += [table([
    ["Step", "What it does"],
    ["1. Filter", "Removes the DC offset and low-passes the audio at 150 Hz."],
    ["2. Envelope", "Takes the envelope (the outline of the waveform's amplitude)."],
    ["3. Match curve", "Slides the envelope of a 6-second airgun template along the segment and records "
                       "how strongly it matches at each position (a cross-correlation, then squared). "
                       "The result rises with anything loud, whatever its shape."],
    ["4. Threshold", "A detection starts wherever the match curve is above a threshold, which is the "
                     "higher of (a) 15% above the segment's median and (b) a fixed floor of 7e15 in raw "
                     "recording units. Stretches less than 1 s apart are joined."],
    ["5. Start, end, duration", "The start and end of each stretch are shifted 3 s later (half the "
                                "template) to line up with the audio. Duration is the time the match "
                                "curve stayed above the threshold."],
    ["6. Noise tests", "Five tests compare the detection with the noise before and after it (rms and "
                       "peak-to-peak level, each way) and check its duration. A detection is rejected "
                       "when two or more tests fail."],
], [1.5, 5.5])]
s += [P("A separate script, <font face='Courier'>merge_airguns_one_file.m</font>, groups detections into "
        "bouts for upload.", body)]
s += [figure(os.path.join(FIGS, "M_template.png"), "<b>The template.</b> Both cross-correlation detectors "
             "compare the audio with the envelope (black) of this 6-second recording of an airgun shot.")]

# ---------------------------------------------------------------- 2 bug fixes
s += [P("2. Bug fixes to the original (pull request #2)", h1),
      P("These are in review at github.com/MarineBioAcousticsRC/Airgun-Detector/pull/2, one commit per fix. "
        "Each commit message carries the full explanation and test results.")]
s += [table([
    ["Fix", "What was wrong", "What changed", "Effect and test"],
    ["Milliseconds in segment start times (readxwavhd.m)",
     "Segment start times were built from the header's seconds only, dropping the milliseconds. "
     "Detection times could be early by up to 999 ms.",
     "Start times now include milliseconds, using the same formula as Triton's rdxwavhd.",
     "Identical to Triton on all 13 files of GOM_AC_02 disk 1 (65,064 segments). The milliseconds "
     "are 0 in those files, so nothing changes for that deployment."],
    ["One-sample offset in the match curve (xcorr_airgun_df100.m)",
     "The segment length used to trim the match curve was one sample short, so every detection was "
     "one sample (0.5 ms) late.",
     "The trim now uses the lags returned by xcorr.",
     "On one full file: 7,091 detections before and after; 6,496 moved one sample earlier; "
     "11 borderline detections dropped out and 11 appeared. A correctness cleanup."],
    ["Bout end times (merge_airguns_one_file.m)",
     "Each bout's end was read from the start-time column, so bouts ended at the start of their last "
     "detection and single-detection bouts had zero length.",
     "Bout ends now use the end-time column.",
     "Synthetic test with three known bouts: 1 of 3 ends correct before, 3 of 3 after."],
], [1.45, 1.95, 1.5, 2.1])]
s += [P("The copies of these files inside Triton (Remoras/Airgun-Detector, and the header reader in "
        "Explosion-Detector and EchoDet) have the same problems and were not changed. None of the fixes "
        "changes which shots are detected in any meaningful way.", small)]

# ---------------------------------------------------------------- 3 labelled data
s += [P("3. Labelled data and scoring", h1),
      P("Everything is measured on five-minute windows in which the shots were marked by eye on the "
        "spectrogram. Marks start from automatic proposals and are corrected in a point-and-click editor "
        "(<font face='Courier'>evaluation/review_marks.m</font>).")]
s += [table([
    ["Set", "Windows", "Marked shots", "Sites", "How the windows were chosen", "Used for"],
    ["Set 1", "%d (%d with airguns)" % (len(R1), T1["n_air"]), T1["shots"], "CE_03, MC_17, DC_15, Y2B_01, AC_02",
     "By hand, to cover known conditions: loud, paired, fast, overlapping and faint shots; ship noise; "
     "whale calls; quiet.", "Testing. Also used to tune the shape-based detector."],
    ["Set 2", "%d (%d with airguns)" % (len(R2), T2["n_air"]), T2["shots"],
     "CE_03, MC_17, DC_15, Y2B_01, SL_03, DT_15, EB_01",
     "At random, four per site.", "Testing. Later used to diagnose two problems of the shape-based detector."],
    ["Set 3", "%d (%d with airguns)" % (len(R3), T3["n_air"]), T3["shots"], "GC_15, MC_17, DC_15",
     "At random, 20 per site, on days used by no other set. No detector took part in choosing them.",
     "Testing only. Labelled once and frozen."],
    ["Train 1", "150 (134 usable)", "2,539", "The seven set 2 sites plus GC_15",
     "From a coarse scan of each deployment, to cover pulsing and non-pulsing stretches.", "Training the network."],
    ["Train 2", "100", "2,638", "The same eight sites",
     "Chosen by the network: windows where it was unsure, or where its copies disagreed.", "Training the network."],
], [0.6, 1.0, 0.65, 1.5, 2.0, 1.25])]
s += bullets([
    "<b>One event per shot cycle.</b> Where a shot reaches the recorder as several overlapping arrivals, "
    "it counts as one shot. In one window (X21) it was not clear whether the pulses in each hump were "
    "separate shots or arrivals of one, so each hump counts as one event.",
    "<b>Scoring.</b> A shot is <i>found</i> if a detection falls within 2 s of it. Two shots under one "
    "detection count as one found and one <i>merged</i>. A detection that holds no marked shot is "
    "<i>extra</i> if the window has airguns and <i>false</i> if it does not. Recall is shots found over "
    "shots marked; precision is detections on a shot over all detections.",
    "<b>Windows the reviewer could not decide.</b> X25 and X28 (set 2) are scored as having no airguns; "
    "the only detection in either is one by the network in X28. In Y01 (set 3) airguns are clearly "
    "present but the shots could not be marked one by one, so it is left out of the shot scores; the "
    "three detectors make 17, 18 and 10 detections in it.",
    "<b>Instrument noise is not tested.</b> Such data is removed before any detector is run. Two set 3 "
    "windows (Y58, Y59) turned out to be instrument noise and were dropped after the first scoring; "
    "they held 15, 55 and 11 detections by the original, shape-based detector and network.",
    "<b>Starting marks for set 3 pooled all three detectors</b> without showing which detector made "
    "which mark, so the labels do not lean towards one of them.",
])

# ---------------------------------------------------------------- 4 original results
by = lambda ids: totals(R1, lambda w: w["window"] in ids)
loud, close, faint = by({"W01", "W02", "W03", "W11"}), by({"W05", "W07", "W08"}), by({"W04", "W18", "W22"})
s += [P("4. What testing showed about the original detector", h1)]
s += [table([
    ["Set 1, by condition", "Shots", "Original", "Shape-based", "Network"],
    ["Loud, regular (W01, W02, W03, W11)", loud["shots"]] + [found(loud, d) for d in DETS],
    ["Close-spaced: fast or overlapping (W05, W07, W08)", close["shots"]] + [found(close, d) for d in DETS],
    ["Faint (W04, W18, W22)", faint["shots"]] + [found(faint, d) for d in DETS],
    ["False detections in the 12 windows without airguns", ""] + [T1[d]["false"] for d in DETS],
], [3.2, 0.6, 1.05, 1.1, 1.05])]
s += [P("Three weaknesses account for the original detector's losses:", body)]
s += bullets([
    "<b>Faint shots never reach the fixed floor.</b> The floor is in raw recording units, and the "
    "background level differs by about 15 dB between sites, so one number cannot suit them all.",
    "<b>Close shots merge.</b> The 6-second template spreads each shot into a wide hump in the match "
    "curve. Shots a few seconds apart run together into one long detection.",
    "<b>Loud backgrounds bypass the floor.</b> When ship noise is loud enough, the \"15% above median\" "
    "rule takes over, and ordinary fluctuations in the noise cross it. All of its false detections in "
    "set 1 are in ship noise.",
])
s += [P("<b>Changing the thresholds is not enough.</b> Lowering the floor recovers faint shots but adds "
        "false detections; raising the offset above the median removes false detections in ship noise but "
        "loses most close-spaced shots. With the floor removed entirely, recall reaches only 0.88 in set 1 "
        "and about 0.70 in set 2, at a precision near 0.8 (the orange curves in section 7). The threshold is "
        "doing two jobs, deciding whether a shot is present and keeping neighbouring shots apart, and no "
        "single value does both.", body)]

# ---------------------------------------------------------------- 5 shape-based
s += [P("5. The shape-based cross-correlation detector", h1),
      P("<font face='Courier'>xcorr_shape/xcorr_airgun_shape.m</font> keeps the filter, the envelope, the "
        "template and the level tests. It changes how candidate shots are found.")]
s += [table([
    ["", "Original", "Shape-based"],
    ["What is measured", "How strongly the signal matches the template. Rises with loudness.",
     "How closely the signal's <i>shape</i> matches the template (a correlation coefficient from "
     "-1 to 1). Does not depend on loudness."],
    ["What counts as a shot", "Any stretch where the match is above the threshold.",
     "Each <i>peak</i> in the shape score. Peaks must be at least 2.5 s apart, so neighbouring shots "
     "stay separate."],
    ["Thresholds", "15% above the segment median, or a fixed floor of 7e15, whichever is higher.",
     "Accepted outright if the shape score is 0.65 or more, or the pulse is loud, or it stands at "
     "least 8 dB above the noise before it. Otherwise, down to a score of 0.40, accepted only if at "
     "least three other detections line up with it at a steady spacing (4 to 40 s). No fixed floor."],
    ["Several arrivals of one shot", "Smeared into one detection by the 6 s template.",
     "Candidates no more than 3.5 s apart are joined into one event, from the first arrival to the "
     "end of the last."],
    ["Segment boundaries", "Each 75 s segment is processed alone, so a shot at a boundary can be cut.",
     "Each segment is read with a few seconds of its neighbours when the recording is continuous."],
    ["Start, end, duration", "From where the match curve crosses the threshold. Depends on the "
     "threshold and the 6 s template.",
     "From the sound itself: where the envelope is within 10 dB of the pulse's peak. Durations are "
     "much shorter as a result, and not comparable with the original's."],
    ["Noise tests", "Four level tests plus duration (too short or too long); rejected when two or "
     "more fail.",
     "The same four level tests, with the noise before measured ahead of the first arrival. Short "
     "duration is no longer a test; only events of 15 s or more count against."],
], [1.35, 2.5, 3.15])]
s += [P("How the two cross-correlation detectors see the same audio", h2),
      P("Each figure below shows one 75-second segment three ways: (a) the low-passed audio and its "
        "envelope, (b) the original detector's match curve with its threshold, and (c) the shape score "
        "with its two thresholds. Shaded bands are the detections each method made. In these panels the "
        "shape-based detector is labelled \"prototype\", its working name.")]
s += [figure(os.path.join(FIGS, "M_fast.png"), "<b>Fast shooting (W07).</b> In (b) the match curve only just "
             "rises above the threshold at each shot and does not fall back between some of them, so shots "
             "merge into long detections. In (c) every shot is a separate peak of the shape score."),
      figure(os.path.join(FIGS, "M_ship.png"), "<b>Ship noise, no airguns (W10).</b> In (b) the threshold sits "
             "15% above the median of the noise, and ordinary fluctuations cross it several times in this "
             "segment. In (c) the shape score never reaches 0.65, and the few peaks above 0.40 are not "
             "part of a train."),
      figure(os.path.join(FIGS, "M_faint.png"), "<b>Faint shots (W22).</b> In (b) the threshold is the fixed "
             "floor, and the match curve stays below it, so the shots are not detected. In (c) the shape "
             "score does not depend on loudness, and the shots appear as regular peaks.")]

# ---------------------------------------------------------------- 6 network
s += [P("6. The spectrogram network", h1),
      P("The third detector (<font face='Courier'>cnn/</font>, in Python with PyTorch) uses no template. A "
        "small convolutional neural network is shown spectrograms together with the reviewer's shot marks, "
        "and learns to output, for every 0.1 s, the probability that a shot starts there.")]
s += [figure(os.path.join(FIGS, "M_network.png"), "<b>Input, layers and output of the network,</b> on two "
             "minutes of window W05. The shot near 44 s scores just under the threshold and is missed; "
             "the other 16 are detected.")]
s += [table([
    ["Part", "What it is"],
    ["Input", "A spectrogram from 5 to 250 Hz in 2.5 Hz steps (99 rows), one column every 0.1 s. The "
              "median level of each frequency over the stretch is subtracted, so the values are dB above "
              "background (clipped to -10 to 40 dB). Site-to-site level differences therefore do not "
              "matter. Nothing else is given: no template, no site, no time of day."],
    ["Shape-finding layers", "Four blocks of 3 x 3 filters over frequency and time. After each block the "
                             "frequency rows are pooled (99 to 49, 24, 12, 4) while every time step is "
                             "kept. They learn what a single shot looks like."],
    ["Context layers", "Seven blocks of filters along time only, with taps spaced 0.1 s to 6.4 s apart, "
                       "so each output can use about 25 s of sound on either side. They learn rhythm and "
                       "what surrounds a shot, and let a long multi-pulse shot count once."],
    ["Output", "One probability per 0.1 s. Peaks above a threshold and at least 2.5 s apart are the "
               "detections. It gives a start time for each shot but no end time or duration."],
    ["Size", "About 208,000 learned values. One copy trains in under half an hour on a desktop processor."],
], [1.5, 5.5])]
s += [P("How it was trained", h2)]
s += bullets([
    "<b>Labels.</b> Each reviewed mark becomes a bump 0.4 s wide on a target curve, and the network is "
    "adjusted until its output matches. 234 labelled windows with 5,177 shots were used (train 1 and "
    "train 2). The test sets are never trained on.",
    "<b>Threshold.</b> A fifth of the training windows are held back, and the threshold that works best "
    "on them is used. The test sets play no part in it.",
    "<b>Averaging.</b> Four copies trained from different random starts are averaged (threshold 0.40). "
    "Single copies vary: on set 2 they find between 422 and 465 shots.",
    "<b>Labelling in rounds.</b> For train 2 the network itself chose the windows it was least sure "
    "about, and the reviewer corrected the pooled output of all three detectors. That second round "
    "raised its score on set 2 from 359 to %d of %d shots." % (T2["network"]["found"], T2["shots"]),
    "<b>Mixing in background noise</b> from windows without airguns during training gave a small gain "
    "(set 2: 451 shots found, with 7 more extra detections) and made the copies more alike. It is not "
    "used in the results reported here.",
])

# ---------------------------------------------------------------- 7 results
s += [P("7. Results", h1)]
rows = [["", "Shots", "Original", "Shape-based", "Network"]]
for name, t in (("Set 1", T1), ("Set 2", T2), ("Set 3", T3)):
    rows += [[name + ": shots found", t["shots"]] + [found(t, d) for d in DETS],
             [name + ": merged / extra / false detections", ""] +
             ["%d / %d / %d" % (t[d]["merged"], t[d]["extra"], t[d]["false"]) for d in DETS],
             [name + ": precision", ""] + [prec(t, d) for d in DETS]]
s += [table(rows, [3.0, 0.7, 1.1, 1.1, 1.1])]
s += [P("How to read the three sets", h2)]
s += bullets([
    "<b>Set 3 is the fair comparison.</b> Random windows at the three sites where results are needed, "
    "labelled once, and not used to tune or train anything.",
    "<b>Sets 1 and 2 favour the shape-based detector.</b> Its settings were adjusted over several rounds "
    "on set 1, and set 2 was used to diagnose two of its problems. The network was never trained or "
    "tuned on any of the three.",
    "<b>The labels of sets 1 and 2 were re-checked</b> with prompts at places where the shape-based "
    "detector or one copy of the network had fired without a mark, so changes made there lean towards "
    "those two detectors.",
])
site_rows = [["Set 3 by site", "Shots", "Original", "Shape-based", "Network"]]
for site in ("GC_15", "MC_17", "DC_15"):
    t = totals(R3, lambda w, site=site: w["site"] == site)
    site_rows.append([site.replace("_", " "), t["shots"]] +
                     ["%s; %d extra, %d false" % (found(t, d), t[d]["extra"], t[d]["false"]) for d in DETS])
s += [table(site_rows, [1.3, 0.6, 1.7, 1.7, 1.7])]
s += [P("At DC_15 the original detector finds almost nothing, and both new detectors are at their "
        "weakest. At GC_15 the network finds nearly every shot.", body)]
s += [figure(os.path.join(FIGS, "PR_curves.png"), "<b>Precision and recall as each detector's settings go "
             "from strict to lenient,</b> on sets 1 and 2. The large open marker is the setting in use. "
             "Original: the fixed floor is lowered from 1e17 to 0, then the offset above the median. "
             "Shape-based: its four acceptance settings are moved together, so this is a path through the "
             "settings and not a single-threshold curve. Network: the probability threshold; the dashed "
             "line is the network before the second labelling round.")]
s += [P("Set 3 by window (found; F = extra or false detections):", body), window_table(R3)]
s += [P("Set 2 by window:", body),
      window_table(R2, {"X04": "fast shooting in a loud background", "X05": "long shots, several arrivals",
                        "X06": "long shots, several arrivals", "X07": "faint, fast, low-frequency",
                        "X21": "scored one event per hump", "X13": "ship noise", "X15": "ship noise"})]
s += [P("Set 1 by window:", body),
      window_table(R1, {w["window"]: w["label"] for w in R1.values()})]

# ---------------------------------------------------------------- 8 examples
s += [P("8. Examples", h1),
      P("Each figure shows one five-minute window. The strip above the spectrogram has four rows: the "
        "marked shots (black triangles), the original detector's detections "
        "(<font color='#d55e00'><b>orange</b></font> bars), the shape-based detector's "
        "(<font color='#0072b2'><b>blue</b></font> bars) and the network's "
        "(<font color='#009e73'><b>green</b></font> ticks; it gives start times only). A bar runs from a "
        "detection's start to its end. The number in each row label is the count. Dotted vertical lines "
        "are segment boundaries.")]
examples = [
    ("set3", R3, "Y12", "Faint shots at DC_15 (Y12, set 3)", "The original's fixed floor is never reached."),
    ("set2", R2, "X07", "Faint, fast shots (X07, MC_17, set 2)", "Again below the original's floor."),
    ("set1", R1, "W07", "Fast shooting (W07, CE_03, set 1)", "Shots are 6.5 s apart, and the original joins "
     "several pairs into one long detection."),
    ("set2", R2, "X15", "Ship noise, no airguns (X15, Y2B_01, set 2)", "Only the original fires."),
    ("set2", R2, "X06", "Long shots with several arrivals (X06, MC_17, set 2)", "Each shot lasts several "
     "seconds. The shape-based detector looks for one sharp pulse and does worst here."),
    ("set2", R2, "X04", "Fast shooting in a loud background (X04, CE_03, set 2)", "The hardest window for "
     "both new detectors; the shots barely stand out."),
    ("set3", R3, "Y02", "Hard for all three (Y02, DC_15, set 3)", "Faint shots that the network scores "
     "below its threshold."),
]
for name, res, wid, head, text in examples:
    s.append(figure(compare(name, wid), "<b>%s.</b> %s Found, %s." % (head, text, cap(res, wid))))

# ---------------------------------------------------------------- 9 limits
s += [P("9. Limits and open questions", h1)]
s += bullets([
    "<b>The marks are one person's judgement.</b> Faint and overlapping shots are uncertain. Ten set 3 "
    "windows are set aside to be labelled a second time, which will show how repeatable the marks are "
    "and so how high any detector's score can be expected to go.",
    "<b>The test sets are small.</b> In set 3, three DC_15 windows hold 65 of the network's 115 misses.",
    "<b>Fast shooting in a loud background (X04)</b> is mostly missed by both new detectors.",
    "<b>Shape-based detector: long multi-pulse shots.</b> It can miss the shot, mark a later pulse "
    "instead of the first, or make two detections in one shot. Its train rule also rejects some real "
    "faint shots when too few neighbours are detected.",
    "<b>Network: confidence on faint shots.</b> At a very low threshold it reaches 97% of the shots in "
    "set 2, but precision falls to 0.69. It notices most of the shots it misses and scores them too low.",
    "<b>Network: practical differences.</b> It needs Python and the trained model files, its results "
    "depend on which model version is used, it gives no durations, and a script to run it over a whole "
    "deployment has not been written yet.",
    "<b>One deployment per site.</b> Training and test windows come from different days of the same "
    "deployments, so performance on a new deployment or region has not been measured.",
    "<b>GC_15:</b> the recording level drops by 55 to 70 dB from mid-October 2022. Windows after that "
    "date were dropped as bad data and none are in set 3.",
    "<b>Instrument noise is not handled</b> by any of the three, and run time on full deployments has "
    "not been measured for the two new detectors.",
])

s += [KeepTogether([P("10. Status and next steps", h1), table([
    ["Item", "Status"],
    ["Bug fixes to the original", "Pull request #2 is open and awaiting review. Nothing has been merged."],
    ["Repository layout", "One folder per detector (xcorr_original, xcorr_shape, cnn) and a shared "
     "evaluation folder, on a local branch. Not committed and not yet proposed for the main repository."],
    ["Shape-based detector", "Tuning stopped. Ready to run on whole deployments (run_xcorr_shape.m)."],
    ["Network", "A third labelling round (train 3, 97 windows, half at DC_15) is prepared, aimed at faint "
     "shots and fast shooting. A whole-deployment run script is still to be written."],
    ["Second pass on set 3", "Ten windows to be labelled again to measure label repeatability."],
], [2.2, 4.8])])]

doc = SimpleDocTemplate(OUT, pagesize=letter, leftMargin=0.75 * inch, rightMargin=0.75 * inch,
                        topMargin=0.7 * inch, bottomMargin=0.8 * inch, title=TITLE, author="Vanessa ZoBell")
doc.build(s, onFirstPage=footer, onLaterPages=footer)
print("wrote", OUT)
for name, t in (("set1", T1), ("set2", T2), ("set3", T3)):
    print(name, t["shots"], {d: (t[d]["found"], t[d]["merged"], t[d]["extra"], t[d]["false"]) for d in DETS})
