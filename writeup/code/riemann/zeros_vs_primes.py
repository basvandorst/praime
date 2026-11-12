"""The primes check published zeros: residual of the band formula, and what it says about the zeros.

For a band of zeros (Odlyzko's tables) and many x, compare
  zero side   S(x) = sum_gamma w(gamma) e^{i gamma ln x}
  prime side  P(x) = -(1/2 pi) sum_{m < 20x + 1000} Lambda(m) m^{-1/2} e^{i c ln(x/m)} What(ln(x/m))
(band_decoder.py, closed-form kernel). An error e_j in zero j changes S by i w_j e^{i gamma_j ln x} e_j ln x,
so random errors of rms eps give |S - P| about eps ln x sqrt(sum w_j^2) and a common shift d gives
|S| d ln x. Calibration: the same residual with every zero moved by a random eps or by a shift d.

usage: python3 zeros_vs_primes.py FILE [FILE ...]
       FILE: an Odlyzko table, or band:T1:W:zeros.f64 (own zeros in (T1, T1 + W], from zeros.js)
"""
import math
import os
import sys

import mpmath as mp
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import band_decoder as bd  # noqa: E402

XS = np.arange(500, 5001, 100)


def prime_side(b, x, ph):
    M = 20 * x + 1000
    tot = 0j
    for p in PRIMES[PRIMES < M]:
        q = int(p)
        lp = math.log(q)
        while q < M:
            v = math.log(x) - math.log(q)
            tot += -1 / (2 * np.pi) * lp / math.sqrt(q) * ph[q] * bd.tukey_ft(v, b.A, b.L)
            q *= int(p)
    return tot


def zero_side(b, x, g):
    return np.sum(b.w * np.exp(1j * g * math.log(x)))


s = np.ones(101001, bool); s[:2] = False
for p in range(2, 318):
    if s[p]:
        s[p * p::p] = False
PRIMES = np.nonzero(s)[0]

for path in sys.argv[1:]:
    b = bd.parse_spec(path)[0]
    two = 2 * mp.pi
    # e^{i c ln(x/m)} = e^{i c ln x} e^{-i c ln m}; S is computed demodulated, so compare S~ with P e^{-i c ln x}
    ph = b.phases(20 * int(XS[-1]) + 1000)
    res, mag = [], []
    rng = np.random.default_rng(1)
    pert = {'rms 1e-7': rng.normal(0, 1e-7, len(b.g)), 'rms 1e-6': rng.normal(0, 1e-6, len(b.g)),
            'shift 1e-7': np.full(len(b.g), 1e-7)}
    rp = {k: [] for k in pert}
    for x in XS:
        P = prime_side(b, int(x), ph)            # demodulated prime side (phase e^{-i c ln m})
        S = zero_side(b, int(x), b.g)            # demodulated zero side
        res.append(abs(S - P)); mag.append(abs(P))
        for k, e in pert.items():
            rp[k].append(abs(zero_side(b, int(x), b.g + e) - P))
    res, mag = np.array(res), np.array(mag)
    L = np.log(XS)
    eps_rms = np.sqrt(np.mean((res / (L * np.sqrt(np.sum(b.w ** 2)))) ** 2))
    shift = np.max(res / (mag * L))
    print(f"{os.path.basename(path) if ':' not in path else path.split('/')[0]}: {len(b.g)} zeros, height {float(b.c):.4e}; x = 500..5000 step 100")
    print(f"  residual |S - P|: rms {np.sqrt(np.mean(res ** 2)):.1e}, max {res.max():.1e} (|P| rms {np.sqrt(np.mean(mag ** 2)):.1f})")
    print(f"  as zero errors: rms per zero <= {eps_rms:.1e}; common shift <= {shift:.1e}")
    for k in pert:
        print(f"  zeros moved by {k}: residual rms {np.sqrt(np.mean(np.array(rp[k]) ** 2)):.1e}")
