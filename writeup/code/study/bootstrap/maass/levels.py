# Smallest admissible level for holomorphic newforms of weight kappa, squarefree level, trivial character:
# q*(kappa) = exp(-t*), conductor q^k for Sym^k, good primes: a measure on the Satake angle; bad primes: Steinberg columns.
# Usage: levels.py <X> <kappa,kappa,...> [bad=1] [N]
import sys, time, numpy as np
from maass import Problem
X = float(sys.argv[1]); bad = (sys.argv[3] != '0') if len(sys.argv) > 3 else True; N = int(sys.argv[4]) if len(sys.argv) > 4 else None
for kap in [float(v) for v in sys.argv[2].split(',')]:
    t0 = time.time(); P = Problem(kap, 0, X, hol=True, bad=bad, N=N, J=121); lam, t, rd = P.solve()
    print('weight %g X %g N %d bad %d: smallest level %.5f (t* = %.6f, %d rounds, %.0fs); forms %s; a_p %s' % (kap, X, P.N, bad, np.exp(-t), t, rd, time.time() - t0,
          ' '.join('%.1e' % v for v in P.mineigs(P.w)), P.summary()), flush=True)
