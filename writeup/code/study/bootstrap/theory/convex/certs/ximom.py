"""Xi-side moments for the archimedean check:
  m_j = (1/2pi) int Xi(t)^2 t^{2j} dt,  a_j = (1/2pi) int Xi(t)^2 t^{2j} (Re psi(1/4+it/2) - log pi) dt,  j=0..J.
Gauss-Legendre, 24 nodes per unit panel on [0, T]. Output xi_mom.pkl.  Usage: python ximom.py [J] [dps] [T] [outfile] [nodes per unit panel]"""
import sys, pickle, time
import mpmath as mp
from core import xi_moments
J = int(sys.argv[1]) if len(sys.argv) > 1 else 16
dps = int(sys.argv[2]) if len(sys.argv) > 2 else 50
T = int(sys.argv[3]) if len(sys.argv) > 3 else 200
mp.mp.dps = dps
t0 = time.time()
NODES = int(sys.argv[5]) if len(sys.argv) > 5 else 24
m, a = xi_moments(J, T=T, panel=1.0, nodes=NODES)
pickle.dump(dict(J=J, dps=dps, T=T, nodes=NODES, m=[mp.nstr(x, dps - 3) for x in m], a=[mp.nstr(x, dps - 3) for x in a]),
            open(sys.argv[4] if len(sys.argv) > 4 else "xi_mom.pkl", "wb"))
print("time", time.time() - t0)
for j in range(J + 1):
    print(j, mp.nstr(m[j], 25), mp.nstr(a[j], 25))
