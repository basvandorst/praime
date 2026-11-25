"""Shape of the r = 1 line: G(u) = sum e^{2i theta(g)} e^{-i u g} on (0, T], predicted -eps c (1/2pi) int_0^T e^{-iut} dt
(the law's own line at r = 1 seen through the sharp window). Slope b(u) per u (registered estimator), genus 2 (100)."""
import sys, os, math; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from g2data import load, slope
rows = load(); T = 150.0
for k in (-3, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3):
    u = k * 2 * math.pi / T
    kern = T if u == 0 else (np.exp(-1j * u * T) - 1) / (-1j * u)
    Gs = [complex((r['ph'] * np.exp(-1j * u * r['gg'])).sum()) for r in rows]; Ps = [-r['eps'] * r['c'] * kern / (2 * math.pi) for r in rows]
    b, se = slope(Gs, Ps)
    print('u = %+5.2f (2pi/T) |kernel|/T %.3f  slope %.3f +- %.3f' % (k, abs(kern) / T, b, se))
