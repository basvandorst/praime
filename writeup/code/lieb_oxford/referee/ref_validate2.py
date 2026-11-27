"""(78) in exact rational arithmetic vs referee shell quadrature (float) at random points, incl. large a / small b."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from fractions import Fraction as Fr
from ref_core import Psi_ball
rng = np.random.default_rng(7)
def psi78_exact(a, b):
    a = Fr(a); b = Fr(b)
    x1 = max(a + b - a*b, 0); x2 = max(abs(a - b) - a*b, 0)
    return (x1**4/160*(a*a*b*b - 5*a*a - 5*b*b + 4*a*b*b + 4*b*a*a + 20*a*b)
            - x2**4/160*(a*a*b*b - 5*a*a - 5*b*b + 4*abs(a - b)*a*b - 20*a*b))
worst = 0
for i in range(20000):
    a = float(np.exp(rng.uniform(np.log(1e-3), np.log(1e4))))
    b = float(np.exp(rng.uniform(np.log(1e-3), np.log(1e1))))
    e = float(psi78_exact(a, b)); m = float(Psi_ball(a, b))
    rel = abs(e - m) / (a**3 * b**3)
    if rel > worst: worst = rel; wp = (a, b, e, m)
print("max |shell - (78) exact| / a^3 b^3 over 20000 points:", worst, wp)
