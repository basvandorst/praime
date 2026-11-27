"""Attack 1: teeth of the Rankin-Selberg constraint. Raw smallest eigenvalue of the truncated Weil matrix of L(f x g)
(no division by k) with true data, flipped d, perturbed a_p(g), shifted r_g, random fake g."""
import sys
sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review3")
import numpy as np
from eng3 import Eng
from sdp3 import cheb_b
from common import load_key, bdict

key = load_key()
F = {round(k["r"], 4): k for k in key}


def rs_min(E, f, g, d=None, apg=None, rg=None, apf=None):
    r, rg = f["r"], (g["r"] if rg is None else rg)
    d = (f["eps"] + g["eps"]) % 2 if d is None else d
    apg = g["ap"] if apg is None else apg
    apf = f["ap"] if apf is None else apf
    kap = [d + 1j * (r + rg), d - 1j * (r + rg), d + 1j * (r - rg), d - 1j * (r - rg)]
    M = E.arch(kap)
    for p, m in E.pp:
        M = M + E.prime(p, m) * cheb_b(apf[p], m)[m] * cheb_b(apg[p], m)[m]
    return min(np.linalg.eigvalsh(M[np.ix_(ix, ix)])[0] for ix in E.blocks)


pairs = [(16.1381, 14.3585), (16.1381, 13.7798), (16.6443, 16.1381), (19.4847, 19.4235), (18.1809, 17.7386),
         (12.1730, 9.5337)]
for N in (24, 48):
    E = Eng(16.0, N)
    print("=== N = %d (Legendre degree < N), X = 16" % N)
    for rf, rg in pairs:
        f, g = F[rf], F[rg]
        base = rs_min(E, f, g)
        flip = rs_min(E, f, g, d=1 - (f["eps"] + g["eps"]) % 2)
        print("f %.4f (%s) x g %.4f (%s): true %+.3e   d flipped %+.3e" % (rf, "odd" if f["eps"] else "even", rg,
              "odd" if g["eps"] else "even", base, flip))
        s = []
        for p in (2, 3, 5, 7, 11, 13):
            for dlt in (+0.3, -0.3):
                ap = dict(g["ap"]); ap[p] = ap[p] + dlt
                if abs(ap[p]) > 2:
                    s.append("p%d%+.1f:out" % (p, dlt)); continue
                s.append("p%d%+.1f:%+.2e" % (p, dlt, rs_min(E, f, g, apg=ap)))
        print("    a_p(g) +- 0.3:", " ".join(s))
        s = []
        for dr in (-0.05, -0.01, +0.01, +0.05):
            s.append("%+.2f:%+.2e" % (dr, rs_min(E, f, g, rg=g["r"] + dr)))
        print("    r_g shifted:", " ".join(s))
    # random fake g with the r_g of a true form, a_p uniform on the Sato-Tate law
    rng = np.random.default_rng(7)
    for rf, rg in pairs[:2]:
        f, g = F[rf], F[rg]
        vals = []
        for trial in range(200):
            th = np.array([rng.choice(np.linspace(0, np.pi, 2001), p=None) for _ in range(6)])
            # Sato-Tate sample: rejection
            ap = {}
            for p in (2, 3, 5, 7, 11, 13):
                while True:
                    t = rng.uniform(0, np.pi)
                    if rng.uniform(0, 1) < np.sin(t) ** 2:
                        break
                ap[p] = 2 * np.cos(t)
            vals.append(rs_min(E, f, g, apg=ap))
        vals = np.array(vals)
        print("random Sato-Tate g (r_g = %.4f) against true f %.4f: fraction with lambda_min < -0.004: %.3f; "
              "quantiles 10/50/90%%: %+.3f %+.3f %+.3f" % (rg, rf, np.mean(vals < -0.004),
                                                            *np.quantile(vals, [0.1, 0.5, 0.9])))
