"""Turtle closed 4-loops: runs (a,b,a,b) close a rectangle, i.e. 5 consecutive primes
p, p+a, p+a+b, p+2a+b, p+2a+2b.  Real counts by admissibility of that 5-tuple vs the
block-shuffle expectation (sum over blocks of 500 gaps of npos * f(a)^2 f(b)^2)."""
import numpy as np
from common import series, N

S = series(N)
g = np.diff(S["primes"])
def admissible(a, b):
    t = [0, a, a + b, 2 * a + b, 2 * a + 2 * b]
    return all(len({x % q for x in t}) < q for q in (2, 3, 5, 7))
real = {}
for k in range(g.size - 3):
    if g[k] == g[k + 2] and g[k + 1] == g[k + 3]:
        key = (int(g[k]), int(g[k + 1])); real[key] = real.get(key, 0) + 1
exp = {}
B = 500
core = g[1:]  # the shuffle kept gap 1 (3-2) in place
for s in range(0, core.size, B):
    blk = core[s:s + B]
    v, c = np.unique(blk, return_counts=True); f = dict(zip(v.tolist(), (c / blk.size).tolist()))
    for a in f:
        for b in f:
            exp[(a, b)] = exp.get((a, b), 0) + blk.size * f[a] ** 2 * f[b] ** 2
tot_r = sum(real.values()); tot_e = sum(exp.values())
adm_r = sum(v for k, v in real.items() if admissible(*k)); adm_e = sum(v for k, v in exp.items() if admissible(*k))
print("total real", tot_r, "shuffle-expected", round(tot_e, 1))
print("admissible patterns: real", adm_r, "expected", round(adm_e, 1))
print("inadmissible patterns: real", tot_r - adm_r, "expected", round(tot_e - adm_e, 1))
print("aaaa loops (5 consecutive primes in AP): real", sum(v for k, v in real.items() if k[0] == k[1]),
      "expected", round(sum(v for k, v in exp.items() if k[0] == k[1]), 1))
top = sorted(exp.items(), key=lambda kv: -kv[1])[:10]
for k, v in top:
    print(k, "adm" if admissible(*k) else "INADM", "expected", round(v, 1), "real", real.get(k, 0))
