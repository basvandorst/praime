"""Small extra readings: slope without the discovery sample (round 151), without the non-USp(4) curve, and their bootstrap."""
import sys, os, math; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from g2data import load, slope
rng = np.random.default_rng(1)
R = load()
for nm, S in (('rounds 152 + rep (85)', [r for r in R if r['rnd'] != '151']), ('without R_N2187_28 (99)', [r for r in R if r['key'] != 'R_N2187_28']),
              ('eps+1 without R_N2187_28', [r for r in R if r['key'] != 'R_N2187_28' and r['eps'] > 0])):
    G = np.array([r['G'] for r in S]); P = np.array([-r['eps'] * r['c'] * 150 / (2 * math.pi) for r in S], complex)
    b, se = slope(G, P); bs = []
    for _ in range(4000):
        i = rng.integers(0, len(S), len(S)); bs.append((G[i] * P[i].conj()).real.sum() / (np.abs(P[i]) ** 2).sum())
    bs = np.array(bs); print('%-26s n=%3d slope %.4f +- %.4f (bootstrap sd %.4f, P(b >= 1) %.4f)' % (nm, len(S), b, se, bs.std(), np.mean(bs >= 1)))
