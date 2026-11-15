"""gaptowers, main computation (numpy).
Q1 gap distribution up to N, champion, twins, maximal gaps.
Q2 gaps after the Mersenne exponents vs random primes of similar size (simulation).
Q3 Gallagher: prime counts in consecutive intervals of length round(ln n) in [1e5, 5e5] vs Poisson(1),
   plus the Hardy-Littlewood prediction of the variance.
Rechecked by gt_check.py (pure Python, sympy)."""
import numpy as np, math

N = 500_000
EXPS = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253,
        4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091]
M = 600_000
s = np.ones(M + 1, dtype=bool); s[:2] = False
for i in range(2, int(M ** 0.5) + 1):
    if s[i]: s[i * i::i] = False
PR = np.nonzero(s)[0]
P = PR[PR <= N]
nP = len(P)
nxt = PR[nP]   # first prime > N
g = np.diff(P)  # gaps between consecutive primes both <= N
print("primes <= N:", nP, " gaps:", len(g), " last prime", P[-1], "next prime", nxt,
      " last column height (truncated at N):", N - P[-1] + 1, " true gap", nxt - P[-1])

# ---------------- Q1
print("\n== Q1 gap distribution")
vals, cnts = np.unique(g, return_counts=True)
print("gap counts:", {int(v): int(c) for v, c in zip(vals, cnts)})
champ = vals[np.argmax(cnts)]
print("jumping champion:", int(champ), "count", int(cnts.max()), " second:", sorted(zip(cnts.tolist(), vals.tolist()))[-2])
print("twin pairs (gap 2):", int(np.sum(g == 2)))
print("mean gap:", float(g.mean()), " ln N:", math.log(N))
rec = []; best = 0
for i, x in enumerate(g):
    if x > best:
        best = x; rec.append((int(x), int(P[i])))
print("maximal gaps (gap, start):", rec)
w = 1415
print(f"w={w}: bands = ceil({nP}/{w}) = {math.ceil(nP / w)}; numbers drawn = {N - 1} (2..N)")
with open("gt_gapcounts.csv", "w") as fh:
    fh.write("gap,count\n")
    for v, c in zip(vals, cnts): fh.write(f"{v},{c}\n")

# ---------------- Q2 Mersenne exponents
print("\n== Q2 gaps after Mersenne exponents")
idx = {int(p): i for i, p in enumerate(PR)}
gall = np.diff(PR).astype(np.float64)       # gap after PR[i]
lnp = np.log(PR[:-1].astype(np.float64))
ratio = gall / lnp
rows = []
for e in EXPS:
    i = idx[e]
    rows.append((e, int(PR[i + 1]), int(gall[i]), float(ratio[i]), int(gall[i - 1]) if i > 0 else None))
for r in rows:
    print("  p=%d next=%d gap=%d gap/ln p=%.3f (gap before %s)" % r)
rng = np.random.default_rng(20261115)
W = 50
def test(exps, label, B=200_000):
    ii = np.array([idx[e] for e in exps])
    S_obs = float(ratio[ii].sum())
    lo = np.maximum(ii - W, 0); hi = ii + W   # window of 101 primes by index (fewer at the start)
    sims = np.zeros(B)
    mu = 0.0; var = 0.0
    for a, b in zip(lo, hi):
        win = ratio[a:b + 1]
        mu += win.mean(); var += win.var()
        sims += win[rng.integers(0, len(win), B)]
    pr_hi = float(np.mean(sims >= S_obs)); pr_lo = float(np.mean(sims <= S_obs))
    p2 = min(1.0, 2 * min(pr_hi, pr_lo))
    # rank test: percentile of each exponent's gap in its window (mid-rank for ties)
    pct = []
    for i, a, b in zip(ii, lo, hi):
        win = gall[a:b + 1]
        pct.append((np.sum(win < gall[i]) + 0.5 * np.sum(win == gall[i])) / len(win))
    print(f"  {label}: n={len(exps)} S=sum g/ln p = {S_obs:.3f}; null mean {mu:.3f} sd {math.sqrt(var):.3f};"
          f" z={(S_obs-mu)/math.sqrt(var):+.2f}; simulated two-sided p = {p2:.3f} (B={B}); mean gap percentile {np.mean(pct):.3f}")
    return S_obs, mu, math.sqrt(var), p2
test(EXPS, "all 31 exponents <= N")
test([e for e in EXPS if e >= 100], "21 exponents >= 100")
# all primes: mean of g/ln p over p in [100, N]
b = (PR[:-1] >= 100) & (PR[:-1] <= N)
print("  mean g/ln p over all primes 100..N:", float(ratio[b].mean()))

# ---------------- Q3 Gallagher
print("\n== Q3 Gallagher intervals")
cum = np.concatenate([[0], np.cumsum(s[:M + 1])])   # cum[x] = #primes < x
starts, lens = [], []
n = 100_000
while True:
    L = int(round(math.log(n)))
    if n + L - 1 > N: break
    starts.append(n); lens.append(L); n += L
starts = np.array(starts); lens = np.array(lens)
cnt = cum[starts + lens] - cum[starts]
J = len(cnt)
print("intervals:", J, " lengths used:", {int(v): int(c) for v, c in zip(*np.unique(lens, return_counts=True))},
      " covered up to", int(starts[-1] + lens[-1] - 1))
mean = float(cnt.mean()); var = float(cnt.var())
print(f"primes counted {int(cnt.sum())}; mean {mean:.4f}; variance {var:.4f}; var/mean {var/mean:.4f}")
hist = np.bincount(cnt)
lam_i = lens / np.log(starts + (lens - 1) / 2.0)
print("expected mean sum h/ln n: %.4f" % lam_i.mean())
# Hardy-Littlewood variance: Var = lam - h/L^2 + (2/L^2) sum_{d<h} (h-d)(S(d)-1), S(d)=singular series of {0,d}
C2 = 0.0
sv = np.ones(10 ** 6, dtype=bool); sv[:2] = False
for i in range(2, 1001):
    if sv[i]: sv[i * i::i] = False
rr = np.nonzero(sv)[0][1:].astype(np.float64)
C2 = float(np.exp(np.sum(np.log1p(-1.0 / (rr - 1) ** 2))))
def Sd(d):
    if d % 2: return 0.0
    f = 2 * C2; m = d
    for r in range(3, d + 1, 2):
        if m % r == 0 and all(r % t for t in range(3, int(r ** 0.5) + 1, 2)):
            f *= (r - 1) / (r - 2)
    return f
hl_var_i = []
for st, h in zip(starts, lens):
    Lg = math.log(st + (h - 1) / 2.0)
    lam = h / Lg
    corr = sum((h - d) * (Sd(d) - 1) for d in range(1, h))
    hl_var_i.append(lam - h / Lg ** 2 + 2 * corr / Lg ** 2)
hl_var = float(np.mean(hl_var_i) + lam_i.var())
print(f"HL-predicted variance {hl_var:.4f} (ratio to mean {hl_var/lam_i.mean():.4f}); Poisson would give {lam_i.mean():.4f}")
print(f"S(d) for d=2,4,6,8,10,12: {[round(Sd(d),4) for d in (2,4,6,8,10,12)]}")
print(" j  observed   Poisson(1)   Poisson(mean)")
chi1 = chi2 = 0.0
pois1 = []; poism = []
for j in range(len(hist) + 1):
    o = int(hist[j]) if j < len(hist) else 0
    e1 = J * math.exp(-1) / math.factorial(j); e2 = J * math.exp(-mean) * mean ** j / math.factorial(j)
    pois1.append(e1); poism.append(e2)
    print(f" {j}  {o:7d}   {e1:9.1f}   {e2:9.1f}")
# chi-square with tail bin >= 5
def chisq(exp):
    obs = list(hist[:5]) + [int(hist[5:].sum())]
    ex = exp[:5] + [J - sum(exp[:5])]
    return sum((o - e) ** 2 / e for o, e in zip(obs, ex))
print(f"chi-square (bins 0..4, >=5) vs Poisson(1): {chisq(pois1):.1f}; vs Poisson(mean): {chisq(poism):.1f}  (5 df)")
np.savetxt("gt_gallagher_hist.txt", hist, fmt="%d")

# ---------------- Q3b full Hardy-Littlewood histogram
# P(count = j) = sum_t (-1)^(t-j) C(t,j) A_t / L^t, A_t = sum over t-subsets T of {0..h-1} of the singular series S(T)
from itertools import combinations
small_pr = np.nonzero(sv)[0]
def sing(T):
    t = len(T)
    if t <= 1: return 1.0
    lg = 0.0
    for p in small_pr:
        nu = len(set(x % p for x in T))
        if nu == p: return 0.0
        if p > 2 * max(T) and nu == t:
            break
        lg += math.log1p(-nu / p) - t * math.log1p(-1.0 / p)
    # tail p > bound: all residues distinct, factor (1-t/p)/(1-1/p)^t, summed in numpy
    big = small_pr[small_pr >= p].astype(np.float64)   # includes the prime at which the loop stopped
    lg += float(np.sum(np.log1p(-t / big) - t * np.log1p(-1.0 / big)))
    return math.exp(lg)
A = {}
for h in (12, 13):
    A[h] = [0.0] * (h + 1)
    for t in range(0, h + 1):
        tot = 0.0
        for T in combinations(range(h), t):
            if len(set(x % 2 for x in T)) > 1: continue   # inadmissible mod 2 (quick skip)
            tot += sing(T)
        A[h][t] = tot
    print(f"h={h}: A_t =", [round(a, 4) for a in A[h]])
pred = np.zeros(8)
for st, h in zip(starts, lens):
    Lg = math.log(st + (h - 1) / 2.0)
    for j in range(8):
        pred[j] += sum((-1) ** (t - j) * math.comb(t, j) * A[h][t] / Lg ** t for t in range(j, h + 1))
print(" j  observed   HL-predicted   Poisson(1)")
for j in range(7):
    o = int(hist[j]) if j < len(hist) else 0
    print(f" {j}  {o:7d}   {pred[j]:10.1f}   {pois1[j]:9.1f}")
obs6 = list(hist[:5]) + [int(hist[5:].sum())]
ex6 = list(pred[:5]) + [J - float(pred[:5].sum())]
obs5 = list(hist[:4]) + [int(hist[4:].sum())]
ex5 = list(pred[:4]) + [J - float(pred[:4].sum())]
po5 = pois1[:4] + [J - sum(pois1[:4])]
print("HL predicted total", float(pred.sum()), " chi-square (bins 0,1,2,3,>=4; 4 df) vs HL:",
      round(sum((o - e) ** 2 / e for o, e in zip(obs5, ex5)), 2), " vs Poisson(1):", round(sum((o - e) ** 2 / e for o, e in zip(obs5, po5)), 1))
print("HL variance from the histogram:", round(float(sum(j * j * pred[j] for j in range(8)) / J - (lam_i.mean()) ** 2), 4),
      " (total variance over all intervals; must equal the pair-sum value", round(hl_var, 4), ")")
