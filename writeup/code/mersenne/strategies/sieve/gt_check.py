"""Independent recheck of gt_main.py: pure Python, sympy (primerange, nextprime, isprime), random module.
Q3 histogram prediction recomputed by explicit inclusion-exclusion over pairs S subset T (not the binomial collapse)."""
import math, random, sympy
from itertools import combinations

N = 500_000
EXPS = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253,
        4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091]
PR = list(sympy.primerange(2, 600_001))
P = [p for p in PR if p <= N]
print("primes <= N", len(P), "last", P[-1], "next", sympy.nextprime(P[-1]))

print("== Q1")
cnt = {}
for a, b in zip(P, P[1:]):
    cnt[b - a] = cnt.get(b - a, 0) + 1
print("gap counts", dict(sorted(cnt.items())))
print("champion", max(cnt, key=cnt.get), cnt[max(cnt, key=cnt.get)], " twins", cnt[2], " gaps", sum(cnt.values()))
best = 0; rec = []
for a, b in zip(P, P[1:]):
    if b - a > best:
        best = b - a; rec.append((b - a, a))
print("maximal gaps", rec)
print("mean gap", (P[-1] - 2) / (len(P) - 1))

print("== Q2")
pos = {p: i for i, p in enumerate(PR)}
def r(i): return (PR[i + 1] - PR[i]) / math.log(PR[i])
for e in EXPS:
    assert sympy.nextprime(e) == PR[pos[e] + 1]
def test(exps, B=100_000, W=50, seed=12345):
    S = sum(r(pos[e]) for e in exps)
    wins = []
    mu = var = 0.0
    for e in exps:
        i = pos[e]
        w = [r(j) for j in range(max(i - W, 0), i + W + 1)]
        m = sum(w) / len(w); mu += m; var += sum((x - m) ** 2 for x in w) / len(w)
        wins.append(w)
    rnd = random.Random(seed)
    ge = le = 0
    for _ in range(B):
        t = sum(rnd.choice(w) for w in wins)
        ge += t >= S; le += t <= S
    p2 = min(1.0, 2 * min(ge, le) / B)
    print(f"  n={len(exps)} S={S:.3f} null mean {mu:.3f} sd {math.sqrt(var):.3f} z={(S-mu)/math.sqrt(var):+.2f} two-sided p={p2:.3f} (B={B})")
test(EXPS)
test([e for e in EXPS if e >= 100])

print("== Q3")
pset = set(PR)
starts = []; n = 100_000
while True:
    L = round(math.log(n))
    if n + L - 1 > N: break
    starts.append((n, L)); n += L
counts = [sum(1 for m in range(s, s + h) if m in pset) for s, h in starts]
J = len(counts)
mean = sum(counts) / J; var = sum((c - mean) ** 2 for c in counts) / J
hist = {}
for c in counts: hist[c] = hist.get(c, 0) + 1
print("intervals", J, "mean", round(mean, 4), "var", round(var, 4), "hist", dict(sorted(hist.items())))
for c in counts:
    pass
# singular series by product over primes < 2e5 (pure Python)
SMALL = list(sympy.primerange(2, 200_000))
def sing(T):
    t = len(T)
    if t <= 1: return 1.0
    f = 1.0
    for p in SMALL:
        nu = len({x % p for x in T})
        if nu == p: return 0.0
        f *= (1 - nu / p) / (1 - 1 / p) ** t
    return f
coef = {}
for h in (12, 13):
    adm = {}
    for t in range(0, 7):
        for T in combinations(range(h), t):
            s_ = sing(T)
            if s_ > 0: adm[T] = s_
    # c[j][t] = sum over S subset T, |S|=j, |T|=t of (-1)^(t-j) S(T)
    c = [[0.0] * 8 for _ in range(8)]
    for T, sT in adm.items():
        t = len(T)
        for j in range(t + 1):
            for S in combinations(T, j):
                c[j][t] += (-1) ** (t - j) * sT
    coef[h] = c
    print(f"h={h}: admissible sets {len(adm)}, largest {max(len(T) for T in adm)}")
pred = [0.0] * 8
for s, h in starts:
    Lg = math.log(s + (h - 1) / 2)
    for j in range(8):
        pred[j] += sum(coef[h][j][t] / Lg ** t for t in range(8))
pois = [J * math.exp(-1) / math.factorial(j) for j in range(8)]
for j in range(7):
    print(f"  {j}: obs {hist.get(j, 0)}  HL {pred[j]:.1f}  Poisson(1) {pois[j]:.1f}")
obs5 = [hist.get(j, 0) for j in range(4)] + [sum(v for k, v in hist.items() if k >= 4)]
hl5 = pred[:4] + [J - sum(pred[:4])]
po5 = pois[:4] + [J - sum(pois[:4])]
chi = lambda o, e: sum((a - b) ** 2 / b for a, b in zip(o, e))
print(f"  chi-square, bins 0,1,2,3,>=4 (4 df): vs HL {chi(obs5, hl5):.2f}; vs Poisson(1) {chi(obs5, po5):.1f}")
lam = [h / math.log(s + (h - 1) / 2) for s, h in starts]
print("  HL variance from histogram coefficients:",
      round(sum(sum(j * j * sum(coef[h][j][t] / math.log(s + (h - 1) / 2) ** t for t in range(8)) for j in range(8)) for s, h in starts) / J
            - (sum(lam) / J) ** 2, 4))
