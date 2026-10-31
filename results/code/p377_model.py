# Hardy-Littlewood count for Puzzle 377 variations #2 (n + 2*10^i, digits <= 7) and #3 (n - 2*10^i, digits >= 2):
# E(k) = #(k-digit strings over the digit set) * P(last digit fits) * prod_q local(q) / prod_i ln N_i.
import math
from sympy import primerange, n_order
QS = list(primerange(7, 200000))
def local_other(k, sign):
    lg = 0.0
    for q in QS:
        nu = len({(-sign * 2 * pow(10, i, q)) % q for i in range(k)})
        lg += math.log1p(-nu / q) - k * math.log1p(-1 / q)
    return math.exp(lg)
def E(k, var):
    sign = 1 if var == 2 else -1
    ndig = 8                                            # digits 0..7 or 2..9
    cand = (7 if var == 2 else 8) * ndig ** (k - 1)     # leading digit nonzero
    lastok = 2 / ndig                                   # last digit 1,7 (var 2) or 3,9 (var 3), for k >= 2
    f25 = (1 / 0.4) ** k                                # all k numbers coprime to 10 (forced by the last digit)
    f3 = (2 / 3) / (2 / 3) ** k                         # all forms equal mod 3: one class of n mod 3 is lost
    N = (4.0 if var == 2 else 5.5) * 10 ** (k - 1)
    return cand * lastok * f25 * f3 * local_other(k, sign) / math.log(N) ** k
if __name__ == "__main__":
    for var in (2, 3):
        print(f"variation #{var}:")
        tot = 0
        for k in range(2, 31):
            e = E(k, var); tot += e if k >= 14 else 0
            if k <= 16 or k % 5 == 0: print(f"  k={k:2d}: expected {e:.3g}")
        print(f"  expected number with 14 or more digits: {tot:.2e}")
