"""mersieve analysis (numpy + gmpy2), reads ms_first.txt and ms_pairs.txt from ms_main.py.

Q1 row counts, Q2 mod-4 law, Q3 Euler's rule, Q4 survivor models, Q5 classes/exponents,
Q6 operation counts.  Rechecked by ms_analyse_check.py (pure Python + sympy).
"""
import numpy as np, gmpy2, math, json

N = 500_000
KS = [16, 64, 256, 1000]
EXPS = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217,
        4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049,
        216091, 756839, 859433, 1257787, 1398269, 2976221, 3021377, 6972593, 13466917,
        20996011, 24036583, 25964951, 30402457, 32582657, 37156667, 42643801, 43112609,
        57885161, 74207281, 77232917, 82589933]
assert len(EXPS) == 51 and len(set(EXPS)) == 51
out = {}

P, F = np.loadtxt("ms_first.txt", dtype=np.int64, unpack=True)
pairs = np.loadtxt("ms_pairs.txt", dtype=np.int64)
pp, kk = pairs[:, 0], pairs[:, 1]
n = len(P)
print("primes:", n, " pairs (all k<=1000 factors):", len(pairs))
idx = {int(p): i for i, p in enumerate(P)}
EX = [e for e in EXPS if e <= N]
print("exponents <= N:", len(EX))

# ---------------- Q1 row counts
print("\n== Q1 row counts")
rowsK = {}
for K in KS:
    f = np.where((F > 0) & (F <= K), F, 0)
    cnt = np.bincount(f, minlength=K + 1)
    surv = int(cnt[0])
    rowsK[K] = cnt
    nonempty = int(np.count_nonzero(cnt[1:]))
    print(f"K={K}: survivors {surv} = {surv/n:.4f}; drawn in rows {n-surv} ({(n-surv)/n:.4f}); nonempty rows {nonempty}")
    out[f"K{K}"] = dict(survivors=surv, frac=surv / n, nonempty=nonempty)
    ex_surv = all((F[idx[e]] == 0) or (F[idx[e]] > K) for e in EX)
    print("   all exponents <= N are survivors:", ex_surv)
cnt = rowsK[1000]
print("row counts k=1..64:", {k: int(cnt[k]) for k in range(1, 65) if cnt[k]})
print("rows k=65..256 total:", int(cnt[65:257].sum()), " rows 257..1000 total:", int(cnt[257:1001].sum()))
top = np.argsort(-cnt[1:])[:12] + 1
print("largest rows (k, count):", [(int(k), int(cnt[k])) for k in top])
with open("ms_rows.csv", "w") as fh:
    fh.write("k,count_smallest_k\n")
    for k in range(1, 1001):
        fh.write(f"{k},{int(cnt[k])}\n")
# cumulative share by row group for K=64
print("K=64: share of drawn primes in rows 1,3,4,5:", {k: round(int(rowsK[64][k]) / (n - int(rowsK[64][0])), 4) for k in (1, 3, 4, 5)})

# ---------------- Q2 mod-4 law, on ALL factors q=2kp+1 (k<=1000), not only the smallest
print("\n== Q2 mod-4 law")
q = 2 * kk * pp + 1
bad8 = int(np.sum((q % 8 != 1) & (q % 8 != 7)))
bad2 = int(np.sum(kk % 4 == 2))
bad1 = int(np.sum((kk % 4 == 1) & (pp % 4 != 3)))
bad3 = int(np.sum((kk % 4 == 3) & (pp % 4 != 1)))
print(f"factors checked {len(q)}: q not +-1 mod 8: {bad8}; k=2 mod 4: {bad2}; k=1 mod 4 with p!=3 mod 4: {bad1}; k=3 mod 4 with p!=1 mod 4: {bad3}")
print("factor counts by k mod 4:", {r: int(np.sum(kk % 4 == r)) for r in range(4)})
print("factor counts by q mod 8:", {r: int(np.sum(q % 8 == r)) for r in range(8)})

# ---------------- Q3 Euler
print("\n== Q3 Euler's rule")
row1 = set(int(p) for p in P[F == 1])
rule = set(int(p) for p in P if p > 3 and p % 4 == 3 and gmpy2.is_prime(2 * int(p) + 1))
print("row 1 size:", len(row1), " rule set size:", len(rule), " equal:", row1 == rule)
print("p=3: 2p+1 = 7 = M_3 (not proper), so 3 is excluded from row 1")

# ---------------- Q4 models
print("\n== Q4 survivor models")
# twin prime constant, product over odd primes < 10^6 (tail < 1e-7)
s = np.ones(10 ** 6, dtype=bool); s[:2] = False
for i in range(2, 1001):
    if s[i]: s[i * i::i] = False
rr = np.nonzero(s)[0][1:].astype(np.float64)
C2 = float(np.exp(np.sum(np.log1p(-1.0 / (rr - 1) ** 2))))
print("C2 =", C2)
small = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]
# per-k HL factor prod_{odd r | k} (r-1)/(r-2)
def hlk(k):
    f = 1.0; m = k
    while m % 2 == 0: m //= 2
    r = 3
    while r * r <= m:
        if m % r == 0:
            f *= (r - 1) / (r - 2)
            while m % r == 0: m //= r
        r += 2
    if m > 1: f *= (m - 1) / (m - 2)
    return f
Pf = P.astype(np.float64)
logS_exact = {K: np.zeros(n) for K in KS}   # exact-q model: prod (1-1/k) over q prime, q=+-1 mod 8, q<M_p
dead_exact = {K: np.zeros(n, dtype=bool) for K in KS}
logS_hl = {K: np.zeros(n) for K in KS}      # Hardy-Littlewood averaged model
adm_count = np.zeros(1001, dtype=np.int64)  # number of admissible (p,k) per k
adm_hits = np.zeros(1001, dtype=np.int64)   # of which q | M_p
adm_expect = np.zeros(1001)                 # expected hits = count/k
hitset = set(zip(pp.tolist(), kk.tolist()))
Mbits = P  # q < 2^p - 1 test only matters for p < 64
for k in range(1, 1001):
    qv = 2 * k * P + 1
    m8 = (qv % 8 == 1) | (qv % 8 == 7)
    proper = np.ones(n, dtype=bool)
    sm = P < 64
    proper[sm] = [int(x) < (1 << int(p)) - 1 for x, p in zip(qv[sm], P[sm])]
    cand = m8 & proper
    for r in small:
        cand &= ~((qv % r == 0) & (qv != r))
    ci = np.nonzero(cand)[0]
    isp = np.array([gmpy2.is_prime(int(x)) for x in qv[ci]], dtype=bool)
    ai = ci[isp]
    adm_count[k] = len(ai)
    adm_hits[k] = sum((int(P[i]), k) in hitset for i in ai)
    adm_expect[k] = len(ai) / k
    # HL averaged probability that q is prime given its forced residues, times 1/k
    piq = np.where(m8 & proper, 2 * C2 * hlk(k) / np.log(qv.astype(np.float64)), 0.0)
    # odd r | k equal to p: factor (p-1)/(p-2) should not be applied; p | k only for p <= 1000, negligible but fix it
    for K in KS:
        if k <= K:
            if k == 1:
                dead_exact[K][ai] = True
            else:
                logS_exact[K][ai] += math.log(1 - 1 / k)
            logS_hl[K] += np.log1p(-np.minimum(piq / k, 1.0 - 1e-12))
# direct test of the 1/k rule
print("1/k test: admissible = q prime, q = +-1 mod 8, q < M_p")
groups = [(1, 1), (3, 3), (4, 4), (5, 5), (7, 7), (8, 8), (9, 9), (11, 12), (13, 16), (17, 32), (33, 64), (65, 256), (257, 1000)]
tab = []
for a, b in groups:
    c = int(adm_count[a:b + 1].sum()); h = int(adm_hits[a:b + 1].sum()); e = float(adm_expect[a:b + 1].sum())
    sd = math.sqrt(sum(adm_count[k] * (1 / k) * (1 - 1 / k) for k in range(a, b + 1)))
    tab.append((a, b, c, h, e, sd))
    print(f"  k {a}-{b}: admissible {c}, divides {h}, expected sum 1/k = {e:.1f} (sd {sd:.1f})")
tot_h = int(adm_hits[1:].sum()); tot_c = int(adm_count[1:].sum())
print("  all pairs dividing:", len(pairs), " of which admissible:", tot_h, "(must be all)")
bands = [(2, 10 ** 4), (10 ** 4, 10 ** 5), (10 ** 5, N + 1)]
print("survivors observed vs models (exact-q, HL-averaged, crude ln(2p)/ln(2Kp))")
q4 = []
for K in KS:
    obs = (F == 0) | (F > K)
    pe = np.where(dead_exact[K], 0.0, np.exp(logS_exact[K]))
    ph = np.exp(logS_hl[K])
    pc = np.log(2 * Pf) / np.log(2 * K * Pf)
    for lo, hi in bands + [(2, N + 1)]:
        b = (P >= lo) & (P < hi)
        nb = int(b.sum()); o = int(obs[b].sum())
        e1 = float(pe[b].sum()); sd1 = float(np.sqrt((pe[b] * (1 - pe[b])).sum()))
        e2 = float(ph[b].sum()); e3 = float(pc[b].sum())
        q4.append((K, lo, hi, nb, o, e1, sd1, e2, e3))
        print(f"  K={K:4d} p in [{lo},{hi}): n={nb:5d} obs {o:5d} ({o/nb:.4f})  exact-q {e1:8.1f} ({e1/nb:.4f}, z={(o-e1)/sd1:+.2f})  HL {e2:8.1f} ({e2/nb:.4f})  crude {e3:8.1f} ({e3/nb:.4f})")

# ---------------- Q5 classes and exponents
print("\n== Q5 classes")
for K in (64, 1000):
    obs = (F == 0) | (F > K)
    pe = np.where(dead_exact[K], 0.0, np.exp(logS_exact[K]))
    for r in (1, 3):
        b = (P % 4 == r)
        print(f"  K={K} p=={r} mod 4: n={int(b.sum())} survivors {int(obs[b].sum())} frac {obs[b].mean():.4f}  exact-q model {pe[b].mean():.4f}")
    for r in (1, 3, 5, 7):
        b = (P % 8 == r)
        print(f"  K={K} p=={r} mod 8: n={int(b.sum())} survivor frac {obs[b].mean():.4f}")
odd = [e for e in EX if e > 2]
print("odd exponents <= N:", len(odd), " =1 mod 4:", sum(e % 4 == 1 for e in odd), " =3 mod 4:", sum(e % 4 == 3 for e in odd))
print("all 50 odd known exponents: =1 mod 4:", sum(e % 4 == 1 for e in EXPS if e > 2), " =3 mod 4:", sum(e % 4 == 3 for e in EXPS if e > 2))
def binom_tail(n, x):  # P(X >= x), X ~ Bin(n, 1/2)
    return sum(math.comb(n, i) for i in range(x, n + 1)) / 2 ** n
a1 = sum(e % 4 == 1 for e in odd); a2 = sum(e % 4 == 1 for e in EXPS if e > 2)
print(f"binomial P(X >= {a1} | n={len(odd)}, 1/2) = {binom_tail(len(odd), a1):.4f};  P(X >= {a2} | n=50, 1/2) = {binom_tail(50, a2):.4f}")
print("odd exponents mod 8:", {r: sum(e % 8 == r for e in EXPS if e > 2) for r in (1, 3, 5, 7)})
# exponents with p = 3 mod 4: is 2p+1 composite (required)?
print("exponents p=3 mod 4, p>3, with 2p+1 prime (must be none):", [e for e in EXPS if e > 3 and e % 4 == 3 and gmpy2.is_prime(2 * e + 1)])
# neighbours of the exponents, rows for K=64
K = 64
nb_rows = []
for e in odd:
    i = idx[e]
    a = int(F[i - 1]) if i > 0 else None; b = int(F[i + 1])
    a = 0 if (a is not None and a > K) else a; b = 0 if b > K else b
    nb_rows.append((e, int(P[i - 1]), a, int(P[i + 1]), b))
print("exponent: (prev prime, row), (next prime, row)  [K=64, row 0 = survivor]")
for t in nb_rows:
    print("  ", t)
neigh = [t[2] for t in nb_rows] + [t[4] for t in nb_rows]
obs64 = ((F == 0) | (F > 64))
print(f"neighbours survivor frac {sum(x == 0 for x in neigh)}/{len(neigh)} vs all primes {obs64.mean():.4f}")
# exponents' survival probability under exact-q model with K=1000 compared with other primes of similar size
pe = np.where(dead_exact[1000], 0.0, np.exp(logS_exact[1000]))
print("exact-q survival prob (K=1000) of odd exponents vs mean over survivors in +-10% window:")
ratios = []
for e in odd:
    if e < 100: continue
    i = idx[e]
    w = (P > 0.9 * e) & (P < 1.1 * e) & obs64 & ((F == 0) | (F > 1000))
    ratios.append(pe[i] / pe[w].mean())
print("  mean ratio over", len(ratios), "exponents >= 100:", round(float(np.mean(ratios)), 3))

# ---------------- Q6 operation counts for the straightforward loop
print("\n== Q6 operation counts")
for K in (16, 64, 94, 128, 256, 1000):
    tested = np.where((F > 0) & (F <= K), F, K)
    tot = int(tested.sum())
    # only k with q = +-1 mod 8 (k=0 mod 4, or k=1 mod 4 & p=3 mod 4, or k=3 mod 4 & p=1 mod 4)
    adm = 0
    for p, t in zip(P.tolist(), tested.tolist()):
        if p == 2:
            adm += t; continue
        c = t // 4  # multiples of 4 up to t
        c += sum(1 for j in range(1, min(t, 4) + 1) if 0) if False else 0
        # count k<=t with k%4==0, plus k%4==(1 if p%4==3 else 3)
        want = 1 if p % 4 == 3 else 3
        c += (t - want) // 4 + 1 if t >= want else 0
        adm += c
    print(f"  K={K}: modpows straightforward {tot} ({tot/n:.1f} per prime); with mod-8 filter {adm}")
json.dump(out, open("ms_analyse.json", "w"), indent=1)
