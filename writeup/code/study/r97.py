"""Round 97 (preregister/round97.txt): the mean square over ln X of the off-diagonal part of R(X, d) is (4/d^2) Q,
Q = sum_{j != k} |c_j|^2 |c_k|^2/(4 + (g_j - g_k)^2) (RH, linear independence). Compared with Q_smooth (partners from the
smooth density) minus the GUE sum-rule deficit sum_j |c_j|^4 K(L_j). usage: r97.py"""
import os, json, math
import numpy as np
from scipy import integrate

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', '..', 'data', 'study')
Z = np.array([float(l) for l in open(os.path.join(HERE, '..', '..', 'data', 'riemann', 'zeros_10000_flint.txt')) if l.strip()])
G = np.concatenate([Z, -Z]); RHO = 0.5 + 1j * G


def K(L):
    f = lambda u: (math.sin(math.pi * u) / (math.pi * u)) ** 2 / (4 + (2 * math.pi * u / L) ** 2) if u else 0.25
    tot = 0.0
    for a in range(0, 400):
        tot += integrate.quad(f, a, a + 1, limit=50)[0]
    tot += 1 / (math.pi ** 2) * 0.5 * (1 / 400) / (4 + 0) * 0        # tail beyond u = 400 is ~ 1/(2 pi^2 u) x kernel: negligible
    return 2 * tot


def q_actual(c2):
    Q = 0.0
    for i0 in range(0, len(G), 2000):
        dg = G[i0:i0 + 2000, None] - G[None, :]
        M = c2[i0:i0 + 2000, None] * c2[None, :] / (4 + dg * dg)
        idx = np.arange(i0, min(i0 + 2000, len(G)))
        M[idx - i0, idx] = 0.0
        Q += float(M.sum())
    return Q


def q_smooth(c2, d):
    # partner density g(t) = |c(t)|^2 N0'(|t|) on a grid, convolved with 1/(4 + x^2) by FFT, read at the zeros
    h = 0.02; Tm = 14000.0
    t = np.arange(-Tm, Tm, h)
    rho = 0.5 + 1j * t
    cc = np.abs(((1 + d) ** rho - 1) / rho) ** 2
    dens = np.where(np.abs(t) > 2 * math.pi, np.log(np.abs(t) / (2 * math.pi)) / (2 * math.pi), 0.0)
    g = cc * dens
    x = np.arange(-len(t) // 2, len(t) // 2) * h
    k = 1 / (4 + x * x)
    n = 1 << int(math.ceil(math.log2(2 * len(t))))
    conv = np.fft.irfft(np.fft.rfft(g, n) * np.fft.rfft(np.fft.ifftshift(np.concatenate([k, np.zeros(0)])), n), n)[:len(t)] * h
    vals = np.interp(G, t, conv)
    return float(np.sum(c2 * vals))


if __name__ == '__main__':
    res = {'zeros': len(Z), 'rows': []}
    Ks = {}
    for d in (0.1, 0.01):
        c2 = np.abs(((1 + d) ** RHO - 1) / RHO) ** 2
        Qa = q_actual(c2); Qs = q_smooth(c2, d)
        L = np.log(np.abs(G) / (2 * math.pi))
        Lr = np.round(L, 2)
        kk = np.array([Ks.setdefault(v, K(v)) for v in Lr])
        pred = float(np.sum(c2 * c2 * kk))
        row = {'delta': d, 'Q_actual': Qa, 'Q_smooth': Qs, 'deficit_measured': Qs - Qa, 'deficit_predicted': pred,
               'ratio': (Qs - Qa) / pred, 'rms_actual': 2 / d * math.sqrt(Qa), 'rms_poisson': 2 / d * math.sqrt(Qs)}
        res['rows'].append(row)
        print({a: (round(b, 6) if isinstance(b, float) else b) for a, b in row.items()}, flush=True)
    r96 = json.load(open(os.path.join(OUT, 'r96.json')))
    res['primes_std_over_round96_window'] = {k: float(np.std(v['R_primes'])) for k, v in r96['rows'].items()}
    print(res['primes_std_over_round96_window'])
    json.dump(res, open(os.path.join(OUT, 'r97.json'), 'w'), indent=1)
