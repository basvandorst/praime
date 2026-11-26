"""Second and fourth moments of ghat_d = Xi^2 B_d for the class-A optima: -G''(0) = (1/2pi) int t^2 Xi^2 B,
G''''(0) = (1/2pi) int t^4 Xi^2 B; and B_d(i/2) (= 2 x pole term). NUMERICAL."""
import sys, json
import os
CERTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'convex', 'certs')
sys.path.insert(0, CERTS)
import mpmath as mp
from core import Basis
mp.mp.dps = 110
J = 26
m = [x[0] for x in Basis(J, 0).eval(0)]
print('m_0..m_3 =', [mp.nstr(x, 12) for x in m[:4]])
for d in range(2, 25, 2):
    P = json.load(open(CERTS + '/polish_A_d%d.json' % d))
    b = [mp.mpf(x) for x in P['b']]
    g0 = mp.fsum(b[j] * m[j] for j in range(d + 1))
    M2 = mp.fsum(b[j] * m[j + 1] for j in range(d + 1))
    M4 = mp.fsum(b[j] * m[j + 2] for j in range(d + 1))
    Bi = mp.fsum(b[j] * (-mp.mpf(1) / 4) ** j for j in range(d + 1))
    print('d=%2d G(0)=%s  -G\'\'(0)=%s  G\'\'\'\'(0)=%s  B(i/2)=%s pole=%s' % (d, mp.nstr(g0, 8), mp.nstr(M2, 8), mp.nstr(M4, 8), mp.nstr(Bi, 8), mp.nstr(Bi / 2, 8)))
