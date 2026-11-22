"""Round 95 (preregister/round95.txt): D(delta), the diagonal over the actual zeros (10,000 zeros from python-flint,
writeup/data/riemann/zeros_10000_flint.txt, plus the smooth-density tail), against Hardy-Littlewood's exact value; the
difference over delta against 7/4 = 2 x 7/8 (the constant of the Riemann-von Mangoldt formula). Post hoc (after the run): the exact weight sqrt(1+d) ln^2(1+d)/d of the
zeros with 1 << gamma << 1/d refines (7/4) d to (7/4) sqrt(1+d) ln^2(1+d)/d. usage: r95.py"""
import json, math, os
import numpy as np
from scipy import integrate

HERE = os.path.dirname(os.path.abspath(__file__))
Z = np.array([float(l) for l in open(os.path.join(HERE, '..', '..', 'data', 'riemann', 'zeros_10000_flint.txt')) if l.strip()])
T = float(Z[-1]); G0 = 0.5772156649015329; B = 1 - G0 - math.log(2 * math.pi)
out = {'zeros': len(Z), 'T': T, 'rows': []}
for d in (0.1, 0.05, 0.02, 0.01, 0.005, 0.002):
    u = math.log1p(d); s = math.sqrt(1 + d)
    w = lambda g: ((1 + d) + 1 - 2 * s * np.cos(g * u)) / (0.25 + g * g)
    zsum = 2 * float(np.sum(w(Z))) / d
    f = lambda t: math.log(t / (2 * math.pi)) / (2 * math.pi) * float(w(t))
    # tail beyond T (after the thirty-first review): int_T^inf w dN = int_T^inf w N0' dt - w(T) S(T+), with
    # S(T+) = n - 1 - theta(T)/pi (n = number of zeros <= T); the remaining -int w' S is far smaller
    import mpmath
    S_T = len(Z) - 1 - float(mpmath.siegeltheta(T)) / math.pi
    per = 2 * math.pi / u; tail = 0.0; a = T
    while a < T + 3000 * per and a < 3e7:
        tail += integrate.quad(f, a, a + per, limit=100)[0]; a += per
    tail += (2 + d) / (2 * math.pi) * (math.log(a / (2 * math.pi)) + 1) / a
    tail -= float(w(T)) * S_T
    D = zsum + 2 * tail / d
    hl = math.log(1 / d) + B + (1 + d) * math.log1p(d) / d - 1
    refined = 1.75 * math.sqrt(1 + d) * math.log1p(d) ** 2 / d          # post hoc approximation (superseded)
    closed = math.log(1 / d) + B + 2 * math.log1p(d) + ((2 + d) / d) * math.log1p(d / 2) - 1   # exact under RH (Weil; review 30)
    row = {'delta': d, 'D': D, 'HL_exact': hl, 'ratio_to_delta': (D - hl) / d, 'predicted': 1.75, 'refined_post_hoc': refined, 'D_minus_HL_minus_refined': D - hl - refined, 'closed_form': closed, 'D_minus_closed': D - closed}
    out['rows'].append(row)
    print('d %.3f  D %.5f  HL %.5f  (D - HL)/d = %.4f  (7/4 = 1.75)' % (d, D, hl, (D - hl) / d), flush=True)
json.dump(out, open(os.path.join(HERE, '..', '..', 'data', 'study', 'r95.json'), 'w'), indent=1)
