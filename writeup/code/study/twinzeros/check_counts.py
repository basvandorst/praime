"""Independent check of tz_counts.py on a small range: recompute every statistic with sympy and plain integers."""
import sys, numpy as np, sympy, math
run = sys.argv[1]
z = np.load(run + ".npz")
X = z["X"]; PA, PB, QA, QB = (int(v) for v in z["meta"])
def cum(a): return np.cumsum(a)
# choose a few grid points and recompute everything up to them
pts = [1, len(X)//3, len(X)//2, len(X)-1]
primes = list(sympy.primerange(2, 2*PB + 40))
pset = set(primes)
for jj in pts:
    x = int(X[jj])
    th = sum(math.log(p) for p in primes if QA < p <= min(x, QB))
    assert abs(cum(z["th"])[jj] - th) < 1e-6*max(1,th), (jj, cum(z["th"])[jj], th)
    for q in (3, 4, 5, 21, 29):
        for r in range(q):
            v = sum(math.log(p) for p in primes if QA < p <= min(x, QB) and p % q == r)
            assert abs(cum(z["c%d" % q][:, r])[jj] - v) < 1e-6*max(1, v)
    for h in (2, 4, 6, 8, 10, 12, 30):
        v = sum(math.log(p)*math.log(p+h) for p in primes if PA < p <= min(x, PB) and p+h in pset)
        c = sum(1 for p in primes if PA < p <= min(x, PB) and p+h in pset)
        assert abs(cum(z["P%d" % h])[jj] - v) < 1e-6*max(1, v) and cum(z["N%d" % h])[jj] == c, (h, jj)
    for name, a, b in (("T026", 2, 6), ("T046", 4, 6)):
        v = sum(math.log(p)*math.log(p+a)*math.log(p+b) for p in primes if PA < p <= min(x, PB) and p+a in pset and p+b in pset)
        assert abs(cum(z[name])[jj] - v) < 1e-6*max(1, v)
    v = sum(math.log(p)*math.log(2*p+1) for p in primes if PA < p <= min(x, PB) and 2*p+1 in pset)
    assert abs(cum(z["SG"])[jj] - v) < 1e-6*max(1, v), ("SG", jj, cum(z["SG"])[jj], v)
    print("grid point", jj, "x =", x, "all statistics agree")
# doubled grid
for jj in pts:
    x2 = 2*int(X[jj])
    th = sum(math.log(p) for p in primes if QA < p <= min(x2, QB))
    assert abs(cum(z["th_2"])[jj] - th) < 1e-6*max(1, th), ("th_2", jj)
    for q in (3, 4, 12, 20):
        for r in range(q):
            v = sum(math.log(p) for p in primes if QA < p <= min(x2, QB) and p % q == r)
            assert abs(cum(z["c%d_2" % q][:, r])[jj] - v) < 1e-6*max(1, v)
    print("doubled grid point", jj, "agrees")
