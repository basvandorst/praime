"""Own analysis: lines of the flowed zeros at integers n, against ln(n) kappa_t(n) n^(-1/2) damping. Phase: rigid shift t pi/8."""
import sys, numpy as np
from sympy import factorint
d = np.load(sys.argv[1]); ts = d['t']; X = d['x']
lo, hi = int(sys.argv[2]), int(sys.argv[3]); NMAX = 200
def kappa_lnn(t):
    b = np.exp(t / 4 * np.log(np.arange(0, NMAX + 1).clip(1)) ** 2); ln = np.log(np.arange(0, NMAX + 1).clip(1))
    K = np.zeros(NMAX + 1)                     # K[n] = kappa_n ln n
    for n in range(2, NMAX + 1):
        s = b[n] * ln[n]
        for dd in range(2, n):
            if n % dd == 0: s -= K[dd] * b[n // dd]
        K[n] = s
    return K
kind = {}
for n in range(2, NMAX + 1):
    f = factorint(n); kind[n] = 'prime' if (len(f) == 1 and list(f.values())[0] == 1) else ('ppower' if len(f) == 1 else 'composite')
for ti in [int(a) for a in sys.argv[4].split(',')]:
    t = ts[ti]; g = X[ti, lo:hi]; g0, g1 = g[0], g[-1]; K = kappa_lnn(t)
    gg = np.linspace(g0, g1, 20001)
    meas = {}; pred = {}; noise = []
    for n in range(2, NMAX + 1):
        F = np.sum(np.exp(1j * (g + t * np.pi / 8) * np.log(n)))
        P = -K[n] / (2 * np.pi * np.sqrt(n)) * np.trapezoid((gg / (2 * np.pi)) ** (-(t / 4) * np.log(n)), gg)
        meas[n] = F; pred[n] = P
        noise.append(np.sum(np.exp(1j * (g + t * np.pi / 8) * np.log(n + 0.5))))
    noise = np.array(noise); sig = np.sqrt(np.mean(noise.real ** 2))
    print(f"t = {t:.2f}, zeros {lo}..{hi} (gamma {g0:.0f}..{g1:.0f}); noise per line (from half-integers) {sig:.1f}; mean imaginary part of noise {np.sqrt(np.mean(noise.imag**2)):.1f}")
    for n in (2, 3, 5, 7, 11, 13, 4, 9, 8, 6, 10, 12, 14, 15, 21, 30):
        print(f"   n = {n:3d} ({kind[n]:9s}): measured {meas[n].real:9.1f} {meas[n].imag:+8.1f}i   predicted {pred[n]:9.1f}   ratio {meas[n].real/pred[n]:6.3f} +- {abs(sig/pred[n]):5.3f}")
    for kd in ('prime', 'ppower', 'composite'):
        ns = [n for n in range(2, NMAX + 1) if kind[n] == kd]
        p = np.array([pred[n] for n in ns]); m = np.array([meas[n].real for n in ns])
        slope = np.sum(p * m) / np.sum(p * p); err = sig / np.sqrt(np.sum(p * p))
        print(f"   pooled slope, {kd:9s} (n <= {NMAX}, {len(ns)} lines): {slope:6.3f} +- {err:5.3f}   (from zero: {slope/err:5.1f} errors)")
