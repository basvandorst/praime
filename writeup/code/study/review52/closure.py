"""The 'prime line alone' comparator (wave closure): zeros with p good = zeros with p missing, warped by the
explicit-formula wave of p; ln|L'| gains the Euler factor A_p and the log of the local density 1 + eps_p.
Contrast (missing - good) = -Var_w(A_p + ln(1 + eps_p)), w = 1 + eps_p (zero weighting), theta uniform.
First order: -Li2(1/p)/2 + (ln p/L) l (2 + l)."""
import numpy as np
from mpmath import polylog

def closure(p, L, K=200, n=1 << 14):
    th = 2 * np.pi * np.arange(n) / n
    k = np.arange(1, K + 1)[:, None]
    c = p ** (-k / 2.0) * np.cos(k * th)
    A = (c / k).sum(0)
    eps = -(2 * np.log(p) / L) * c.sum(0)
    w = 1 + eps
    X = A + np.log(w)
    m = np.mean(w * X)
    full = -(np.mean(w * X * X) - m * m)
    # same without ln(1+eps) linearised: X1 = A + eps
    X1 = A + eps; m1 = np.mean(w * X1)
    lin = -(np.mean(w * X1 * X1) - m1 * m1)
    l = np.log(p / (p - 1.0))
    first = -0.5 * float(polylog(2, 1.0 / p)) + np.log(p) / L * l * (2 + l)
    return full, lin, first, eps.min()

if __name__ == '__main__':
    for L in (9.0, 10.2, 10.72, 11.4):
        print('L =', L)
        for p in (2, 3, 5, 7, 11, 13):
            f, li, fo, em = closure(p, L)
            print('  p=%2d closure full %+.5f  (eps linear %+.5f)  first order %+.5f  min(1+eps) %.3f' % (p, f, li, fo, 1 + em))
