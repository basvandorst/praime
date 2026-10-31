# Analytic model for Puzzle 24: n prime, M_b = (base-b digits of n read in base 10) prime for every b in B.
# Correlations between bases enter only through the last digit (n mod b decides M_b mod 2 and mod 5) and,
# for 3 | b-1, through M_b = s_b(n) = n (mod 3). Everything else is treated as independent (Mertens factor 1).
import numpy as np, math
from math import gcd
def joint(B):
    L = 1
    for b in B: L = L * b // gcd(L, b)
    L = L * 30 // gcd(L, 30)
    rs = [r for r in range(L) if gcd(r, 30) == 1]               # n prime > 5: unit mod 30, uniform over units mod L
    rs = [r for r in rs if gcd(r, L) == 1]
    tot = 0.0
    for r in rs:
        w = 1.0
        for b in B:
            if b == 10: continue
            d0 = r % b
            if b > 10 and d0 >= 10: w = 0; break
            if d0 % 2 == 0 or d0 % 5 == 0: w = 0; break
            w /= 0.4                                             # survive 2 and 5: 1 / ((1-1/2)(1-1/5))
            for q in (3, 7, 11, 13):                             # q | 10 - b: M_b = n (mod q), never 0
                if (10 - b) % q == 0: w *= q / (q - 1)
        tot += w
    return tot / len(rs)
def count(A, Bx, B, n=400):
    J = joint(B); xs = np.geomspace(A, Bx, n); f = np.zeros(n)
    for i, x in enumerate(xs):
        lx = math.log(x); p = J / lx
        for b in B:
            if b == 10: continue
            p /= lx * math.log(10) / math.log(b)                 # ln M_b ~ ln n * ln10 / ln b
            if b > 10: p *= (10 / b) ** (lx / math.log(b) - 1)   # all non-last digits < 10
        f[i] = p
    return np.trapezoid(f, xs), J
if __name__ == "__main__":
  for B, Xs in ((list(range(2, 11)), (5.0e10, 2.8e11, 3.1e12)), (list(range(2, 12)), (1e14, 1.814e14, 3e14)),
              (list(range(2, 13)), (6e15, 1e17, 1e18, 1e19, 1e20, 1e21))):
      for X in Xs:
          N, J = count(1e6, X, B)
          print(f"bases 2..{B[-1]}: joint factor {J:.2f}, expected N(<{X:.3g}) = {N:.3g}")
