import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'convex', 'certs'))
import mpmath as mp
from core import Basis
mp.mp.dps = 50
B = Basis(4, 0)
for u in [2, 3, 4, 5]:
    v = B.eval(u)
    print(u, [mp.nstr(v[j][0] / v[0][0] / (-4 * mp.pi ** 2 * mp.exp(2 * u)) ** j, 8) for j in range(1, 5)])
