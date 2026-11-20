"""The airgun network: spectrogram, network layers, loading the trained copies
and picking detections. Used by run_network.py.

The network reads a spectrogram and gives, for every 0.1 s, the probability
that an airgun shot starts there. Two parts:
  1. 2-D convolutions over (frequency, time) that learn what a shot looks like
     across frequency. Frequency is pooled away; time resolution is kept.
  2. 1-D convolutions over time with growing gaps (dilations), so each output
     can use about 25 s of context on either side, e.g. the rhythm of a survey.

The trained network (models/network_airgun.pt) holds four copies trained from
different random starts; their probabilities are averaged.
"""
import os

import numpy as np
import torch
import torch.nn as nn
from scipy import signal
from scipy.signal import find_peaks

FS = 2000                 # sample rate of the df100 data (Hz)
HOP_S = 0.1               # time step of the spectrogram and of the network output (s)
NPERSEG = 800             # 0.4 s analysis window -> 2.5 Hz resolution
FMIN, FMAX = 5.0, 250.0   # frequency range given to the network (Hz)
MODELS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "models")


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


class Block2d(nn.Module):
    def __init__(self, cin, cout, pool_f):
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv2d(cin, cout, 3, padding=1), nn.BatchNorm2d(cout), nn.ReLU(inplace=True),
            nn.MaxPool2d((pool_f, 1)))

    def forward(self, x):
        return self.net(x)


class Block1d(nn.Module):
    def __init__(self, ch, dilation):
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv1d(ch, ch, 5, padding=2 * dilation, dilation=dilation), nn.BatchNorm1d(ch), nn.ReLU(inplace=True))

    def forward(self, x):
        return x + self.net(x)


class ShotNet(nn.Module):
    def __init__(self, n_freq):
        super().__init__()
        self.front = nn.Sequential(Block2d(1, 16, 2), Block2d(16, 32, 2), Block2d(32, 48, 2), Block2d(48, 64, 3))
        f_out = n_freq // 2 // 2 // 2 // 3
        self.squeeze = nn.Sequential(nn.Conv1d(64 * f_out, 64, 1), nn.BatchNorm1d(64), nn.ReLU(inplace=True))
        self.context = nn.Sequential(*[Block1d(64, d) for d in (1, 2, 4, 8, 16, 32, 64)])
        self.head = nn.Conv1d(64, 1, 1)

    def forward(self, x):
        """x: (batch, n_freq, n_time) -> logits (batch, n_time)"""
        h = self.front(x.unsqueeze(1))
        h = h.flatten(1, 2)
        h = self.context(self.squeeze(h))
        return self.head(h).squeeze(1)


def load(name="network_airgun"):
    """The copies of a trained network (models/<name>.pt or a full path) and its threshold."""
    path = name if os.path.exists(name) else os.path.join(MODELS, name + ".pt")
    ck = torch.load(path, map_location="cpu")
    models = []
    for state in ck["states"]:
        m = ShotNet(ck["n_freq"])
        m.load_state_dict(state)
        m.eval()
        models.append(m)
    return models, float(ck["threshold"])


def pick(prob, threshold, min_sep_s=2.5):
    """Detection times (s) from a probability curve: peaks above threshold, at least min_sep_s apart."""
    pk, _ = find_peaks(prob, height=threshold, distance=max(1, int(round(min_sep_s / HOP_S))))
    return pk * HOP_S
