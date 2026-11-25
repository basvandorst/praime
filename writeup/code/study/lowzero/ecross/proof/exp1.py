"""Experiment 1: zeros vs primes for one curve, several weights (sharp, Gaussian-edged, Hann)."""
import os, sys, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, cypari2
import gl2lib as G

cid = sys.argv[1] if len(sys.argv) > 1 else 'N307_0_0_1_m8_m9'
C = G.load_curves()[cid]
N, w, ainv, bad = C['N'], C['w'], C['ainv'], C['bad']
gam = G.load_zeros(cid)
pari = cypari2.Pari(); pari.allocatemem(600 * 10 ** 6)
X = int(N * 402.0 ** 2 / (4 * math.pi ** 2)) + 10
t0 = time.time()
n, Lam, _ = G.prime_power_data(pari, ainv, bad, X)
print(cid, 'N', N, 'w', w, 'X', X, 'prime powers', len(n), '%.1fs' % (time.time() - t0))
wins = [(40, 100), (100, 160), (160, 220), (220, 280), (280, 340), (340, 395)]
for kind, D in [('box', 0.0), ('box', 2.0), ('box', 6.0), ('hann', 0.0)]:
    for (a, b) in wins:
        H = G.Weight(kind, a, b, D)
        Z = G.zero_side(N, gam, H)
        P1 = G.prime_side(N, n, Lam, H, mode='exact')
        P2 = G.prime_side(N, n, Lam, H, mode='exact', second_order=True) if D > 0 or kind == 'hann' else P1
        P3 = G.prime_side(N, n, Lam, H, mode='simple')
        print('%s D=%4.1f [%3d,%3d]: zeros %8.3f%+8.3fi | primes exact %8.3f%+8.3fi (2nd %8.3f%+8.3fi) simple %8.3f%+8.3fi | |Z-P| %.3f %.3f %.3f' % (
            kind, D, a, b, Z.real, Z.imag, P1.real, P1.imag, P2.real, P2.imag, P3.real, P3.imag, abs(Z - P1), abs(Z - P2), abs(Z - P3)))
