"""Classical estimate of the noise level of prime11.html: the stations a/q carry the mass (mu(q)^2 / phi(q)) * sum 1/ln^2 n
(major arcs), the total mass is N; noise = (1 - sum_{q <= Q} mu^2/phi(q) * sum_{n<=X} 1/ln^2 n / N) / (share of bins kept).
Compared with the measured levels of spec_1e6.txt / spec2_1e6.txt.   usage: noise_formula.py"""
import math
import numpy as np
from sympy import totient, mobius, primepi
X = 10 ** 6; N = int(primepi(X))
n = np.arange(2, X + 1, dtype=float); s2 = float(np.sum(1 / np.log(n) ** 2))
H = 2 ** 21 // 2 + 1
for Q, kept, meas in ((50, 1032709, 0.5968), (150, 1017499, 0.5234), (400, 912949, 0.4704)):
    m = sum(int(mobius(q)) ** 2 / int(totient(q)) for q in range(1, Q + 1))
    print('Q = %d: station mass share %.4f, predicted noise %.4f, measured %.4f' % (Q, m * s2 / N, (1 - m * s2 / N) / (kept / H), meas))
