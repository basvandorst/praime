"""Independent checks for the band results (RIEMANN.md section 12).

A. The formula. For Odlyzko's band near zero #10^22, at a few x: the zero side
   sum_gamma w(gamma) e^{i gamma ln x} in mpmath from the decimal strings of the table, against the
   prime side -(1/2 pi) sum_{m < 20x + 1000} Lambda(m) m^{-1/2} K(ln x - ln m) with
   K(v) = int w(t) e^{i t v} dt integrated numerically (no closed form), and the primes from sympy.
B. The decoded lists written by band_decoder.py --out: every number below the first wrong
   decision is prime (sympy.isprime) and every prime there except the wheel primes is in the list.

usage: python3 band_check.py formula ODLYZKO_FILE
       python3 band_check.py list FILE [FILE ...]
"""
import os
import re
import sys

import mpmath as mp
import numpy as np
from sympy import isprime, primerange, primepi

mp.mp.dps = 40


def formula(path):
    txt = open(path).read().splitlines()
    base = int(re.search(r'gamma - (\d+)', ' '.join(txt[:6])).group(1))
    offs = [s.strip() for s in txt if re.fullmatch(r'\s*\d+\.\d+\s*', s)]
    d = np.array([float(s) for s in offs])
    lo, hi = d[0], d[-1]
    W, L = hi - lo, (hi - lo) / 7

    def w(t):  # t: offset (float)
        if t < lo + L:
            return 0.5 * (1 - np.cos(np.pi * (t - lo) / L))
        if t > hi - L:
            return 0.5 * (1 + np.cos(np.pi * (t - (hi - L)) / L))
        return 1.0
    wz = [w(t) for t in d]
    gam = [mp.mpf(base) + mp.mpf(s) for s in offs]
    # K(v) = int_lo^hi w(t) e^{i (base + t) v} dt, by Simpson on a fine grid of the offset
    tg = np.linspace(lo, hi, 40001)
    wg = np.array([w(t) for t in tg])
    sw = np.ones_like(tg); sw[1:-1:2] = 4; sw[2:-1:2] = 2
    sw *= (tg[1] - tg[0]) / 3

    worst = 0
    for x in [101, 997, 1000, 2003, 4999]:
        lx = mp.log(x)
        zs = mp.fsum(mp.mpf(wz[j]) * mp.expj(gam[j] * lx) for j in range(len(gam)))
        zs = complex(zs)
        qs, ls = [], []
        for p in primerange(2, 20 * x + 1000):
            q = p
            while q < 20 * x + 1000:
                qs.append(q); ls.append(float(mp.log(p)))
                q *= p
        vs = [lx - mp.log(q) for q in qs]
        ph0 = np.array([complex(mp.expj(mp.mpf(base) * v)) for v in vs])
        vf = np.array([float(v) for v in vs])
        ks = np.empty(len(vf), complex)
        for a in range(0, len(vf), 200):
            ks[a:a + 200] = np.exp(1j * np.outer(vf[a:a + 200], tg)) @ (sw * wg)
        ps = complex(np.sum(-1 / (2 * np.pi) * np.array(ls) / np.sqrt(np.array(qs, float)) * ph0 * ks))
        worst = max(worst, abs(zs - ps))
        print(f"x = {x}: zero side {zs:.6f}, prime side {ps:.6f}, difference {abs(zs - ps):.1e}")
    print(f"largest difference {worst:.1e}")


def check_list(path):
    lines = open(path).read().splitlines()
    head = lines[0]
    wrong = [int(s) for s in re.findall(r'\d+', head.split('first wrong')[1])]
    cold = re.search(r'cold=(\d+)', head)
    lo = int(cold.group(1)) if cold else 3
    hi = wrong[0] if wrong else int(re.search(r'N=(\d+)', head).group(1)) + 1
    got = [int(s) for s in lines[1:] if s.strip()]
    got = [p for p in got if lo <= p < hi]
    bad = [p for p in got if not isprime(p)]
    wm = re.search(r'wheel=([\d,]+)', head)
    wheel = [int(q) for q in wm.group(1).split(',')] if wm else [2]
    want = int(primepi(hi - 1) - primepi(lo - 1)) - sum(1 for q in wheel if lo <= q < hi)
    print(f"{os.path.basename(path)}: {len(got)} primes in [{lo}, {hi}), all prime: {not bad}; "
          f"primes in range other than the wheel primes: {want}; complete: {len(got) == want and not bad}")


if __name__ == '__main__':
    if sys.argv[1] == 'formula':
        formula(sys.argv[2])
    else:
        for f in sys.argv[2:]:
            check_list(f)
