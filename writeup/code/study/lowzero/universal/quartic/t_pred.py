import sys, os, time; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, mpmath as mp, pred4 as P
mp.mp.dps = 20
# zeta check
for s in (0.6 + 3.7j, 0.9 + 400.2j, 1.2 + 2500j, 1.7 + 30j):
    print('zeta_em', s, complex(P.zeta_em(np.array([s]))[0]), complex(mp.zeta(s)))
print('A1 =', P.A1_exact(), ' Q(1/2) =', P.Q_real(0.5), ' Q(1/3) =', P.Q_real(1/3))
w = np.array([0.3 + 5j, 0.3 + 200j]); l = 3.0; X = np.exp(-w * np.log(l))
print('closed P_l vs sum over types: ', P.P_l(l, X), P.E_local(l, X) * (1 - X**2) * (1 - X**3) * (1 - X**4) * (1 + l**-2 - l**-3 - l**-4))
t0 = time.time(); m = P.quartic_model(T=200.0); print('model T=200: %.1f s' % (time.time() - t0))
print('residue of D at 1/2 (claim 1/2):', m.poles[0][1], '  at 1/3:', m.poles[1][1])
for r2, mus in P.SIG_MU.items(): print('r2=%d tails' % r2, P.tail_coeffs(m, mus))
