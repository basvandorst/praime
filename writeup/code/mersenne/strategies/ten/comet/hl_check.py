"""Check of the Hardy-Littlewood integral I(n) = int_2^{n-2} dt / (ln t ln(n-t)): the 4001-point trapezoid (used before)
against scipy quad and against a trapezoid in v = ln t on [ln 2, ln(n/2)] (symmetry), 4001 points.  usage: hl_check.py"""
import math
import numpy as np
from scipy.integrate import quad
C2 = 0.6601618158468696
def trap_t(n):
    t = np.linspace(2, n - 2, 4001); f = 1 / (np.log(t) * np.log(n - t)); return float(np.sum((f[1:] + f[:-1]) / 2) * (t[1] - t[0]))
def trap_v(n, K=4000):
    v = np.linspace(math.log(2), math.log(n / 2), K + 1); t = np.exp(v); f = t / (np.log(t) * np.log(n - t)); return 2 * float(np.sum((f[1:] + f[:-1]) / 2) * (v[1] - v[0]))
def q(n):
    g = lambda t: 1 / (math.log(t) * math.log(n - t)); pts = [2, 3, 10, 100, n / 2]
    return 2 * sum(quad(g, a, b, limit=200)[0] for a, b in zip(pts[:-1], pts[1:]) if b > a)
for n, r, S in ((100000, 1620, 4 / 3), (1000000, 10804, 4 / 3)):
    a, b, c = trap_t(n), trap_v(n), q(n)
    print('n = %d: I by old trapezoid %.3f, by trapezoid in ln t %.3f, by quad %.3f; HL = 2 C2 S(n) I: old %.2f, new %.2f, quad %.2f; r = %d, r/HL(quad) = %.5f' % (n, a, b, c, 2*C2*S*a, 2*C2*S*b, 2*C2*S*c, r, r / (2*C2*S*c)))
