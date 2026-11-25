"""chartools.py: per-character statistics and the two surgeries, shared by the development analysis and freeze_test.py
(freeze_test.py carries its own copy of these definitions, so it is self-contained; this module is for development).

Per character (values |L(1/2+it, chi_D)| on t_k = 100 + k/20, k = 0..18000):
  M, lnZ_b (b = 0.5, 1, 1.5, 2, 3), mean_ln : means over the 143 intervals of length 2 pi (obslib.py);
  varM      : variance (ddof 1) of M over the 143 intervals;
  lnEZ_1, lnEZ_2 : ln of the mean over the intervals of Z_1, Z_2 (annealed: ln of the mean square / fourth power).
Surgery for a prime p (E_p(t) = 1 - chi_D(p) p^(-1/2 - it), the Euler factor's inverse):
  removal  (p good): statistic of |L| * |E_p|  minus statistic of |L|            ('p missing minus p present')
  addition (p bad) : statistic of |L|  minus statistic of |L| / |1 - e p^(-1/2-it)|, averaged over e = +1, -1.
"""
import numpy as np
import obslib

PR = (2, 3, 5, 7)
KEYS = ['M', 'lnZ_0.5', 'lnZ_1', 'lnZ_1.5', 'lnZ_2', 'lnZ_3', 'mean_ln', 'varM', 'lnEZ_1', 'lnEZ_2']


def kron(D, p):
    if D % p == 0:
        return 0
    if p == 2:
        return 1 if D % 8 in (1, 7) else -1
    return 1 if pow(D % p, (p - 1) // 2, p) == 1 else -1


def stats(v):
    o = obslib.observables(v)
    s = {k: float(o[k].mean()) for k in ['M', 'lnZ_0.5', 'lnZ_1', 'lnZ_1.5', 'lnZ_2', 'lnZ_3', 'mean_ln']}
    s['varM'] = float(o['M'].var(ddof=1))
    s['lnEZ_1'] = float(np.log(np.mean(np.exp(o['lnZ_1']))))
    s['lnEZ_2'] = float(np.log(np.mean(np.exp(o['lnZ_2']))))
    return s


def surgery(D, v, s0=None):
    """returns {p: {'kind': 'rem'|'add', key: contrast}}"""
    g = obslib.grid()
    if s0 is None:
        s0 = stats(v)
    out = {}
    for p in PR:
        c = kron(D, p)
        if c != 0:
            s = stats(v * np.abs(1 - c * p ** -0.5 * np.exp(-1j * g * np.log(p))))
            out[p] = {'kind': 'rem'}
            out[p].update({k: s[k] - s0[k] for k in KEYS})
        else:
            acc = {k: 0.0 for k in KEYS}
            for e in (1, -1):
                s = stats(v / np.abs(1 - e * p ** -0.5 * np.exp(-1j * g * np.log(p))))
                for k in KEYS:
                    acc[k] += 0.5 * (s0[k] - s[k])
            out[p] = {'kind': 'add'}
            out[p].update(acc)
    return out
