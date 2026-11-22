"""Round 93: B from the zeros' counting function. The ln X-average of R(X, d) is the diagonal sum
D(d) = (1/d) sum_rho |(1 + d)^rho - 1|^2/|rho|^2. With the Riemann-von Mangoldt density N'(t) = ln(t/2pi)/(2pi) (both
signs of t), D(d) = (2/d) int (1/2pi) ln(t/2pi) |(1 + d)^(1/2+it) - 1|^2/(1/4 + t^2) dt, and for small d
|(1 + d)^(1/2+it) - 1|^2 ~ 2(1 - cos(t d)) gives, with int_0^inf ln t (1 - cos at)/t^2 dt = (pi a/2)(1 - gamma - ln a),
D(d) = ln(1/d) + 1 - gamma - ln 2pi = ln(1/d) + B. Checked here numerically. usage: r93_B.py"""
import json, math, os
import mpmath as mp

mp.mp.dps = 30
G0 = float(mp.euler); B = 1 - G0 - math.log(2 * math.pi)
out = {'B': B, 'identity': {}, 'D_minus_ln': {}}
for a in (0.5, 1.0, 3.0):
    lhs = mp.quadosc(lambda t: mp.log(t) * (1 - mp.cos(a * t)) / t ** 2, [0, mp.inf], omega=a)
    out['identity'][str(a)] = {'integral': float(lhs), 'formula': float(mp.pi * a / 2 * (1 - mp.euler - mp.log(a)))}
for d in (1e-1, 1e-2, 1e-3, 1e-4):
    u = mp.log(1 + d); s = mp.sqrt(1 + d)
    f = lambda t: mp.log(t / (2 * mp.pi)) / (2 * mp.pi) * ((1 + d) + 1 - 2 * s * mp.cos(t * u)) / (mp.mpf(1) / 4 + t ** 2)
    val = 2 / d * mp.quadosc(f, [2 * mp.pi, mp.inf], omega=u)     # density used from t = 2pi (where it is 0)
    out['D_minus_ln'][str(d)] = float(val - mp.log(1 / d))
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study', 'r93_B.json'), 'w'), indent=1)
