"""Integrated version of the gap model: N(d, x) ~ S(d) * integral_2^x exp(-d / ln t) / ln^2 t dt (each prime t has a gap d
with chance about S(d) exp(-d/ln t)/ln t), so the ratio N(30)/N(6) is S(30)/S(6) times a ratio of integrals.  The
single-x model of champ.py uses ln x for every prime; here each prime uses its own ln t.  usage: champ_model.py out.txt"""
import sys, math
from scipy.integrate import quad
out = open(sys.argv[1], 'w')
def I(d, x):  # integral in u = ln t
    return quad(lambda u: math.exp(u - d / u) / u ** 2, math.log(3), math.log(x), limit=400)[0]
for e in range(3, 9):
    x = 10 ** e
    r62 = 2 * I(6, x) / I(2, x); r306 = (4 / 3) * I(30, x) / I(6, x)
    s = '%6.0e  N(6)/N(2) integrated model %.3f   N(30)/N(6) integrated model %.3f' % (x, r62, r306); print(s); out.write(s + '\n')
