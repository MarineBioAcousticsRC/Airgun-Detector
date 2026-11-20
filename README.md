:ocean: Scripps Acoustic Ecology Lab

# Airgun-Detector

Tools for detecting seismic airgun shots in decimated HARP recordings (2 kHz sample rate, `df100` xwavs).

There are three detectors. All three give **one detection per airgun shot** and write their results in the same layout: a folder `<deployment>_AirgunDetections` with one `.mat` per xwav.

| Folder | Detector | How it decides | Language |
|---|---|---|---|
| [`xcorr_original/`](xcorr_original/) | Original cross-correlation detector | Cross-correlates the signal envelope with an airgun template and thresholds how *large* the result is | MATLAB |
| [`xcorr_shape/`](xcorr_shape/) | Shape-based cross-correlation detector | Same template, but scores how well the *shape* matches regardless of loudness, then checks for a regular shot train | MATLAB |
| [`cnn/`](cnn/) | Spectrogram network (CNN) | A small neural network reads the spectrogram and gives the probability of a shot start every 0.1 s | Python (PyTorch) |

**The network is recommended.** It finds the most shots at every site tested, at the same precision as the others.

## How well they work

The detectors were tested on two deployments at each of three Gulf of Mexico sites. From each deployment, four random five-minute windows per month were chosen, and every airgun shot in them was marked by eye. Every detection that had no mark was then checked again by eye.

A shot counts as *found* if a detection falls within 2 s of it. Recall is shots found / shots marked. Precision is detections on a shot / all detections.

| Site | Deployments | Shots | Original | Shape-based | Network |
|---|---|---|---|---|---|
| GC | GC_17, GC_18 | 3,589 | recall 0.53, precision 1.00 | recall 0.72, precision 0.99 | recall 0.84, precision 0.99 |
| MC | MC_19, MC_20 | 1,790 | recall 0.27, precision 0.98 | recall 0.75, precision 0.99 | recall 0.80, precision 0.98 |
| DC | DC_17, DC_18 | 1,867 | recall 0.03, precision 1.00 | recall 0.62, precision 0.99 | recall 0.83, precision 0.99 |

Things to know:

- **The original detector misses most shots at quiet sites.** Its threshold is fixed in raw recording units, and background levels differ by 15 dB or more between sites.
- **Most shots the network misses are very faint.**
- **Remove data with instrument noise or dropouts before running any detector**, or remove those periods from the results afterwards. None of the detectors was built to reject every kind of instrument noise.
- **The network can count a later arrival (echo) of a shot** as a second shot when the two are more than 2.5 s apart.

## Original (`xcorr_original/`)

`xcorr_airgun_df100.m` processes each 75 s segment of an xwav:

1. Low-pass at 150 Hz and take the envelope.
2. Cross-correlate with the envelope of a 6 s airgun template (`air_template_df100.mat`) and square the result.
3. A detection is any stretch above the higher of (a) 15% above the segment median and (b) a fixed floor, `minCorrThresh = 7e15`.
4. Five tests on level and duration against the noise before and after; a detection is rejected when two or more fail.

`merge_airguns_one_file.m` groups detections into bouts.

**To run:** edit `deplList` and `outDisk` at the top of `xcorr_airgun_df100.m`, then run it in MATLAB with Triton's `utFindFiles` on the path (`Triton\Remoras\Logger`).

## Shape-based cross-correlation (`xcorr_shape/`)

`xcorr_airgun_shape.m` keeps the filter, envelope, template and level tests of the original, and changes how shots are found:

- **Shape score.** A correlation coefficient (-1 to 1) between the smoothed envelope and the template, so loudness does not matter. Each peak, at least 2.5 s from the next, is a candidate shot.
- **Acceptance.** A candidate is accepted if its shape score is 0.65 or more, or it is loud, or it stands 8 dB above the noise before it. Weaker candidates (down to 0.40) are kept only if at least three other detections line up with them at a steady spacing of 4 to 40 s.
- **Several arrivals of one shot.** Candidates within 3.5 s are joined into one event.
- **No fixed floor**, so it works across sites with different background levels.

**To run** on a whole deployment, in MATLAB with Triton's `utFindFiles` on the path:

```matlab
addpath('xcorr_shape')
run_xcorr_shape('GOM_MC_17','E:\AirgunShape')                                    % a deployment under \\frosty.ucsd.edu\GofMX_Decimated_8
run_xcorr_shape('//frosty.ucsd.edu/GofMX_Decimated_12/GOM_GC_18','E:\AirgunShape') % or the full path of a deployment on any share
```

Output has the same variables as the original detector, plus `allTier` and `allCand`. Files that already have a result are skipped, so an interrupted run can be started again.

## Spectrogram network (`cnn/`)

- **Input.** A spectrogram from 5 to 250 Hz with one column every 0.1 s. Each frequency has its median level removed, so the network sees how far a sound stands above the background.
- **Network.** About 208,000 learned values: four 2-D convolution blocks that learn what a shot looks like, then seven 1-D convolution blocks along time that give each output about 25 s of context on either side (enough to use the regular rhythm of a survey).
- **Output.** The probability that a shot starts, every 0.1 s. Peaks above 0.40 and at least 2.5 s apart are the detections. It gives **start times only**, not durations.
- **Trained network.** `models/network_airgun.pt` holds four copies trained from different random starts, whose probabilities are averaged. It was trained on 318 hand-labelled five-minute windows (6,551 shots) from eight HARP deployments (GC_15, MC_17, DC_15, DT_15, EB_01, SL_03, Y2B_01, CE_03). None of the test deployments above were used in training.

| File | What it is |
|---|---|
| `run_network.py` | Runs the network over whole deployments |
| `model.py` | Spectrogram, network layers, loading the trained network |
| `models/network_airgun.pt` | The trained network and its threshold |
| `requirements.txt` | Python packages needed |

**To run** (Python 3; install the packages once with `pip install -r cnn/requirements.txt`):

```bash
cd cnn
python run_network.py GOM_GC_15 GOM_MC_17 --out E:/AirgunDetections_network                                  # under //frosty.ucsd.edu/GofMX_Decimated_8
python run_network.py GOM_GC_18 --root //frosty.ucsd.edu/GofMX_Decimated_12 --out E:/AirgunDetections_network  # on another share
```

- It uses a GPU if there is one. On a desktop processor it runs at roughly 1,000 times real time.
- Files that already have a result are skipped, so an interrupted run can be started again.
- `--threshold` changes the threshold (default 0.40). Lower finds more faint shots but gives more false detections.

**Output.** The `.mat` files hold the same variables as the original detector (`allSmpPts`, `allExp`, `bt`, `rawStart`, `rawDur`, `parm`, ...).
- The level and duration measures the network does not make are filled with NaN.
- Each detection's end is set equal to its start.
- `allProb` holds the network's probability.

Three files are added per deployment:
- `<deployment>_detections.csv`: every detection with its time and probability.
- `<deployment>_minutes_per_day.csv`: per day, the minutes recorded and the minutes with airguns.
- `run_info.txt`: the model, threshold and settings used.
