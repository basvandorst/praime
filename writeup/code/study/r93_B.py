"""Round 93: the diagonal with the smooth zero density, after the twenty-ninth review (mpmath quadosc was unreliable;
the lower limit 2pi was arbitrary). D0(d) = (2/d) int_0^inf N0'(t) |(1+d)^(1/2+it) - 1|^2/(1/4 + t^2) dt with
N0'(t) = ln(t/2pi)/(2pi), integrated period by period with scipy. Compared with Hardy-Littlewood's exact finite-d value
ln(1/d) + B + (1+d) ln(1+d)/d - 1, and with the diagonal over the actual zeros (r93.json). usage: r93_B.py"""
import json, math, os
import numpy as np
from scipy import integrate

G0 = 0.5772156649015329; B = 1 - G0 - math.log(2 * math.pi)
out = {'B': B, 'rows': []}
for d in (1e-1, 1e-2, 1e-3, 1e-4):
    u = math.log1p(d); s = math.sqrt(1 + d)
    f = lambda t: math.log(t / (2 * math.pi)) / (2 * math.pi) * ((1 + d) + 1 - 2 * s * math.cos(t * u)) / (0.25 + t * t)
    per = 2 * math.pi / u; tot = 0.0; a = 0.0
    # first period with care near 0 (log singularity is integrable), then period by period, then an analytic tail
    edges = [0.0, min(per, 50.0)]
    while edges[-1] < 4000 * per and edges[-1] < 2e7:
        edges.append(edges[-1] + per)
    for lo, hi in zip(edges[:-1], edges[1:]):
        tot += integrate.quad(f, lo, hi, limit=200)[0]
    Tt = edges[-1]
    tot += (2 + d) / (2 * math.pi) * (math.log(Tt / (2 * math.pi)) + 1) / Tt          # tail, cosine averaged out
    D0 = 2 / d * tot
    hl = math.log(1 / d) + B + (1 + d) * math.log1p(d) / d - 1
    out['rows'].append({'delta': d, 'D0_minus_ln_minus_B': D0 - math.log(1 / d) - B, 'HL_finite_delta_term': hl - math.log(1 / d) - B})
    print('d %.0e  D0 - ln(1/d) - B = %.6f   HL finite-d term = %.6f' % (d, D0 - math.log(1 / d) - B, hl - math.log(1 / d) - B), flush=True)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study', 'r93_B.json'), 'w'), indent=1)
