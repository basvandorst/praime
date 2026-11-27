# Turns an optimiser output (c3opt_d<deg>.json or c3ref_d<deg>.json) into a certificate file cert_d<deg>.json:
# b_k rounded to exact decimals with DIG significant digits; if B is not > 0 on [0, oo) in s = t^2 (Sturm), eps is added to b_0.
# Usage: c3mkcert.py SRC DEG [DIG]
import sys; sys.dont_write_bytecode = True
import json
from fractions import Fraction
from decimal import Decimal, getcontext
import sympy as sp
getcontext().prec = 200
src, deg = sys.argv[1], int(sys.argv[2]); DIG = int(sys.argv[3]) if len(sys.argv) > 3 else 40
d = json.load(open(src))
b = [Fraction(Decimal(x)) for x in d['b']]
scale = b[0]; b = [x/scale for x in b]                       # normalise b_0 = 1 (scale-free problem)
def rnd(x):
    if x == 0: return Fraction(0)
    s = format(Decimal(x.numerator)/Decimal(x.denominator), '.%de' % (DIG - 1)); return Fraction(Decimal(s))   # no float
b = [rnd(x) for x in b]
s = sp.symbols('s')
def npos(bb):
    q = sp.Poly([sp.Rational(x.numerator, x.denominator) for x in reversed(bb)], s)
    return q.count_roots(0, None) if q.degree() > 0 else 0
assert b[-1] > 0, 'leading coefficient of B is not positive: B < 0 for large t'
eps = Fraction(0); it = 0
while npos(b) > 0 or b[0] <= 0:
    it += 1; assert it < 14, 'B could not be made positive by a small eps'
    eps = Fraction(1, 10**12) if eps == 0 else eps*10
    b[0] += eps
print('deg %d: %d coefficients, signs %s, eps added %s' % (deg, len(b), ''.join('+' if x > 0 else '-' for x in b), eps))
out = {'deg': deg, 'source': src.split('/')[-1], 'b_exact': [str(x) for x in b], 'eps_added_to_b0': str(eps), 'NC': 400, 'prec': 600,
       'claim_ratio': d.get('claim_ratio', '1'), 'claim_bound': d.get('claim_bound', '0')}
json.dump(out, open('/home/bas/praime/writeup/code/study/bootstrap/maass/theory/cert_d%d.json' % deg, 'w'), indent=1)
