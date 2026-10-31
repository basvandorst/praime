# Exotic look-back gap rules: g = f(g1, g2, g3, ...), backtested near 1e9.
# Train: primes in [X-T, X); test: primes in [X, X+T). Candidate = first number >= p+g
# with no prime factor <= y (y = 0: no sieve, candidate = p+g).
import numpy as np, itertools, sys
X, T = 10**9, 3*10**6
lo, hi = X - T - 2000, X + T + 2000
base = np.ones(40000, bool); base[:2] = False
for i in range(2, 200):
    if base[i]: base[i*i::i] = False
bp = np.nonzero(base)[0]
L = hi - lo + 1
isp = np.ones(L, bool)
for q in bp:
    if q*q > hi: break
    s = (-lo) % q; isp[s::q] = False
P = lo + np.nonzero(isp)[0]
def surv_next(y):
    s = np.ones(L, bool)
    for q in bp:
        if q > y: break
        st = (-lo) % q; s[st::q] = False
    s |= isp & (np.arange(L) + lo <= y)   # small primes themselves survive
    idx = np.where(s, np.arange(L), L)
    return np.minimum.accumulate(idx[::-1])[::-1]   # first survivor at or after
G = np.diff(P)                     # G[i] = P[i+1]-P[i]
K = 8                              # need g1..g8
i0 = np.arange(K, len(P) - 1)      # predict P[i+1] from P[i]
p = P[i0]; nxt = P[i0 + 1]
g = {k: G[i0 - k] for k in range(1, K + 1)}   # g[1] = last gap (P[i]-P[i-1])
train = p < X; test = ~train
def ev(gv):
    return np.maximum(2, 2*np.round(gv/2)).astype(np.int64)
rules = {}
# linear family a1 g1 + a2 g2 + a3 g3 + c
A = [-2, -1, -0.5, 0, 0.5, 1, 2]
for a1, a2, a3 in itertools.product(A, A, A):
    for c in (-4, -2, 0, 2, 4):
        rules[f"{a1:g}·g1 {a2:+g}·g2 {a3:+g}·g3 {c:+d}"] = (a1, a2, a3, c)
g1, g2, g3 = g[1].astype(float), g[2].astype(float), g[3].astype(float)
exo = {
 "g1 − 2(g1 − g2)  [= 2g2 − g1]": 2*g2 - g1,
 "g1 + (g1 − g2)  [rechte lijn]": 2*g1 - g2,
 "g1 + 2(g1 − g2)": 3*g1 - 2*g2,
 "3g1 − 3g2 + g3  [parabool]": 3*g1 - 3*g2 + g3,
 "g1 + g2 − g3": g1 + g2 - g3,
 "g2 + g3 − g1": g2 + g3 - g1,
 "g1 + g2  [Fibonacci]": g1 + g2,
 "|g1 − g2|": np.abs(g1 - g2),
 "|g1 − g2| + 2": np.abs(g1 - g2) + 2,
 "max(g1, g2)": np.maximum(g1, g2),
 "min(g1, g2, g3)": np.minimum(np.minimum(g1, g2), g3),
 "(g1 + g2)/2": (g1 + g2)/2,
 "g1·g2 mod 30": (g1*g2) % 30,
 "g1² mod 30": (g1**2) % 30,
 "(g1 + g2 + g3) mod 12": (g1 + g2 + g3) % 12,
 "30 − (g1 + g2 + g3) mod 30": 30 - (g1 + g2 + g3) % 30,
 "6 − g1 mod 6": 6 - g1 % 6,
 "g2  [periode 2]": g2, "g3  [periode 3]": g3, "g4": g[4].astype(float), "g6": g[6].astype(float), "g8": g[8].astype(float),
 "geometrisch g1²/g2": g1**2/g2,
 "harmonisch 2g1g2/(g1+g2)": 2*g1*g2/(g1+g2),
 "als g1 = g2 dan 6 anders g1": np.where(g1 == g2, 6, g1),
 "als g1 < g2 dan g2 anders g1 − 2": np.where(g1 < g2, g2, g1 - 2),
 "g1 + g2 − 2·g3 + g4": g1 + g2 - 2*g3 + g[4],
 "mediaan(g1..g8)": np.median(np.stack([g[k] for k in range(1, 9)]), axis=0),
}
out = {}
for y in (0, 7, 13, 100, 1000):
    nx = surv_next(y) if y else None
    def score(gv):
        c = p + ev(gv)
        if y:
            j = nx[np.minimum(c - lo, L - 1)]; c = np.where(j < L, lo + j, -1)
        h = c == nxt
        return 100*h[train].mean(), 100*h[test].mean()
    res = []
    for name, (a1, a2, a3, c) in rules.items():
        res.append((name,) + score(a1*g1 + a2*g2 + a3*g3 + c))
    for name, gv in exo.items():
        res.append((name,) + score(gv))
    res.append(("ZEEF: eerste overblijver",) + score(np.full(len(p), 2.0)))
    if y == 0:
        for yy in (3, 5):
            n2 = surv_next(yy); j = n2[np.minimum(p + 1 - lo, L - 1)]; print(f"  ter vergelijking: wiel tot {yy}: test {100*(lo + j == nxt)[test].mean():.1f}")
    out[y] = res
    best_tr = max(res[:-1], key=lambda r: r[1])
    zeef = res[-1]
    print(f"\n=== zicht {y}: {len(res)-1} regels; zeef {zeef[2]:.1f}% (test)")
    print(f"  beste op training: {best_tr[0]}  train {best_tr[1]:.1f}  test {best_tr[2]:.1f}")
    for r in sorted(res[:-1], key=lambda r: -r[2])[:6]: print(f"  top test: {r[0]:40s} {r[2]:.1f}")
    for r in res:
        if r[0] in list(exo)[:4]: print(f"  {r[0]:40s} test {r[2]:.1f}")
    print("  regels boven de zeef (test):", sum(r[2] > zeef[2] for r in res[:-1]))
import json; json.dump({str(k): v for k, v in out.items()}, open("exo.json", "w"))
