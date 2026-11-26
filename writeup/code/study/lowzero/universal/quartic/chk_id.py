import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
from local4 import T, F, l, X, coeffs
for n in T:
    mu, lam = coeffs(n)
    for lv in (5, 7):
        x = sp.Rational(1, 3)  # generic X
        lhs = sum(F(n, v).subs(l, lv) * x**v for v in range(-3, 22))
        rhs = sp.Rational(lv, lv - 1) * sum(mu[k] * (1 / (lv * x))**k for k in range(4)) * sum(lam[k] * x**k for k in range(26))
        print(n, lv, float(lhs - rhs), end='; ')
    print()
