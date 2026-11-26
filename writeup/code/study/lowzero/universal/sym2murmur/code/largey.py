import sys
sys.path.insert(0, sys.argv[1])
import numpy as np
from model import M_nat, nu, nu_local
print("mean nu over m<=1e5:", np.mean([nu(m) for m in range(1, 100001)]))
print("nu(8) =", nu(8), " (route C average 1.3848)")
for k in (2, 4, 6):
    ys = np.array([1.0, 3.0, 10.0, 30.0, 100.0, 300.0]) + 0.137
    mp = M_nat(ys, k, hyperbolic=False)
    print("k=%d  y:" % k, ys, "\n   M^P/sqrt(y):", np.round(mp / np.sqrt(ys), 3), " 6/(k-1) = %.3f" % (6 / (k - 1)),
          "\n   M^nat (with hyperbolic):", np.round(M_nat(ys, k), 3))
