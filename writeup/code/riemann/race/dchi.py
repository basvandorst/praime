"""Mean (over ln X) of the variance ratio R_chi(X, d) of psi(x, chi) increments over (x, (1+d)x], from the
Guinand-Weil formula with h(r) = |(1+d)^(1/2+ir) - 1|^2/(1/4 + r^2) (Fourier support |x| <= ln(1+d), empty prime sum):
  D_chi(d) = (1/(2 pi d)) int h(r) [ln(q/pi) + Re psi((1/2 + a)/2 + i r/2)] dr   (+ 2 ln(1+d) for zeta, the pole)
Checked on zeta against the closed form of round 95. usage: dchi.py"""
import math, json, os
import numpy as np
import mpmath as mp
from scipy import integrate

G0 = float(mp.euler)


def D(d, q, a, pole=False):
    u = math.log1p(d); s = math.sqrt(1 + d)
    h = lambda r: ((1 + d) + 1 - 2 * s * math.cos(r * u)) / (0.25 + r * r)
    f = lambda r: h(r) * (math.log(q / math.pi) + float(mp.re(mp.digamma((0.5 + a) / 2 + 0.5j * r))))
    per = 2 * math.pi / u; tot = integrate.quad(f, 0, per, limit=200)[0]; lo = per
    while lo < 4000 * per:
        tot += integrate.quad(f, lo, lo + per, limit=100)[0]; lo += per
    # tail: h ~ (2+d)/r^2 on average, Re psi ~ ln(r/2): integral of (2+d)(ln(q/pi) + ln(r/2))/r^2 from lo
    tot += (2 + d) * (math.log(q / math.pi) + math.log(lo / 2) + 1) / lo
    val = 2 * tot / (2 * math.pi * d)
    return val + (2 * math.log1p(d) if pole else 0.0)


closed_zeta = lambda d: math.log(1 / d) - G0 - math.log(2 * math.pi) + 2 * math.log1p(d) + ((2 + d) / d) * math.log1p(d / 2)
out = {}
for d in (0.1, 0.01):
    z = D(d, 1, 0, pole=True)
    c4 = D(d, 4, 1); c3 = D(d, 3, 1)
    out[str(d)] = {'zeta_weil': z, 'zeta_closed': closed_zeta(d), 'chi4': c4, 'chi3': c3,
                   'chi4_minus_zeta_without_pole_plus_lnq': c4 - (closed_zeta(d) - 2 * math.log1p(d) + math.log(4))}
    print(d, {k: round(v, 6) for k, v in out[str(d)].items()}, flush=True)
json.dump(out, open('dchi.json', 'w'), indent=1)
