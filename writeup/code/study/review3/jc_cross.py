"""Round 19 check: crossings of gap 30 over gap 6 from the exact HL inclusion-exclusion sums U_j(h) (own enumeration)."""
import json, math
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
E = json.load(open('/home/bas/scratch/study-review3/enum.json'))
U = {int(h): np.array(v) for h, v in E['U'].items()}
def share(h, lam):
    return sum((-1) ** j * U[h][j] / lam ** (j + 1) for j in range(len(U[h])))
diff = lambda u: share(30, u) - share(6, u)
loc = brentq(diff, 40, 200)
print('local crossing log10 x = %.4f' % (loc / math.log(10)))
def Ndiff(Uend, U0):
    # (N(30) - N(6)) * e^-Uend: integral of e^(u-Uend) diff(u)/u du
    f = lambda u: math.exp(u - Uend) * diff(u) / u
    return quad(f, U0, Uend, limit=400, epsabs=0, epsrel=1e-12)[0]
for U0x in (1e3, 1e8, 1e20):
    cu = brentq(lambda Ue: Ndiff(Ue, math.log(U0x)), 60, 120, xtol=1e-12)
    x = math.exp(cu)
    print('cumulative crossing from %.0e: x = %.6e, log10 %.5f, vs 1.74274e35: %+.4f %%' % (U0x, x, cu / math.log(10), (x / 1.74274357e35 - 1) * 100))
# ratio at 1e8 from 1e3 (exact law) vs data 222847/768752
U8 = math.log(1e8)
N = {h: quad(lambda u: math.exp(u - U8) * share(h, u) / u, math.log(1e3), U8, limit=400)[0] for h in (6, 30)}
print('N30/N6 at 1e8 (exact sums, from 1e3): %.4f, data 0.2899' % (N[30] / N[6]))
print('share6, share30 at lam=81:', share(6, 81.), share(30, 81.))
