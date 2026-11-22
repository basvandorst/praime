"""Round 119 (preregister/round119.txt): Var/mean of prime counts in consecutive H-intervals per block of 1e8 (rounds 10
and 34), from the zeros of zeta: counts - H/ln(mid) ~ -2 Re sum [li((x+H)^rho) - li(x^rho)], plus the diagonal tail.
usage: r119.py"""
import os, json, math
import numpy as np
from scipy import integrate
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, '..', '..', 'data', 'study')
Z = np.array([float(l) for l in open(os.path.join(OUT, '..', 'riemann', 'zeros_10000_flint.txt')) if l.strip()])
T = float(Z[-1]); RHO = 0.5 + 1j * Z


def li_rho(lnx):                       # li(x^rho) for all zeros at the points with log lnx (array), asymptotic Ei
    w = np.outer(lnx, RHO)
    return np.exp(w) / w * (1 + 1 / w + 2 / w ** 2 + 6 / w ** 3 + 24 / w ** 4)


def tail(x, H):
    d = H / x; u = math.log1p(d); s = math.sqrt(1 + d)
    f = lambda t: ((1 + d) + 1 - 2 * s * math.cos(t * u)) / (0.25 + t * t) * math.log(t / (2 * math.pi)) / (2 * math.pi)
    per = 2 * math.pi / u; lo = T; tot = 0.0
    while lo < T + 40 * per:
        tot += integrate.quad(f, lo, lo + per, limit=100)[0]; lo += per
    tot += (2 + d) * (math.log(lo / (2 * math.pi)) + 1) / (2 * math.pi * lo)
    return 2 * x * tot / math.log(x) ** 2           # variance of the count contribution (both signs of gamma)


res = {}
for L0, fname in ((10 ** 9, 'r10_var.json'), (3 * 10 ** 9, 'r10_var_3000000000.json'), (5 * 10 ** 9, 'r10_var_5000000000.json')):
    meas = json.load(open(os.path.join(OUT, fname)))
    for H in (100000, 1000000):
        vals = []
        for b in range(10):
            Lb = L0 + b * 10 ** 7 * 10
            x = Lb + np.arange(10 ** 8 // H) * H
            S = np.zeros(len(x))
            for i0 in range(0, len(Z), 500):
                sl = slice(i0, i0 + 500); r = RHO[sl]
                a = np.exp(np.outer(np.log(x + H), r)) / np.outer(np.log(x + H), r); a *= 1 + 1 / np.outer(np.log(x + H), r)
                c = np.exp(np.outer(np.log(x), r)) / np.outer(np.log(x), r); c *= 1 + 1 / np.outer(np.log(x), r)
                S += -2 * np.real((a - c).sum(axis=1))
            mid = x + H / 2; e = H / np.log(mid)
            v_low = float(np.var(S, ddof=1)); v_tail = float(np.mean([tail(float(xx), H) for xx in (Lb, Lb + 5 * 10 ** 7, Lb + 10 ** 8)]))
            vals.append((v_low + v_tail) / e.mean())
        m = float(np.mean(meas[str(H)]['meas'])); th = float(np.mean(meas[str(H)]['theory'])); z = float(np.mean(vals))
        res['%d_%d' % (L0, H)] = {'measured': m, 'theory_HL': th, 'zeros': z, 'rel_err': z / m - 1, 'zeros_nearer': abs(z - m) < abs(th - m)}
        print('window %.0e H %.0e: measured %.4f  HL %.4f  zeros %.4f  (rel %+.3f)' % (L0, H, m, th, z, z / m - 1), flush=True)
json.dump(res, open(os.path.join(OUT, 'r119.json'), 'w'), indent=1)
