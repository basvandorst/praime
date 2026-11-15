"""Independent recheck of ms_analyse.py: pure Python ints, sympy.isprime / primerange, no numpy, no gmpy2.
Reads chk_first.txt / chk_pairs.txt written by ms_check.py (pure-Python pow loop, p outer, k inner)."""
import math, sympy

N = 500_000
KS = [16, 64, 256, 1000]
EXPS = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217,
        4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049,
        216091, 756839, 859433, 1257787, 1398269, 2976221, 3021377, 6972593, 13466917,
        20996011, 24036583, 25964951, 30402457, 32582657, 37156667, 42643801, 43112609,
        57885161, 74207281, 77232917, 82589933]
first = {}
for line in open("chk_first.txt"):
    p, k = map(int, line.split()); first[p] = k
pairs = []
for line in open("chk_pairs.txt"):
    p, k = map(int, line.split()); pairs.append((p, k))
plist = sorted(first)
assert plist == list(sympy.primerange(2, N + 1))
n = len(plist)
print("primes", n, "pairs", len(pairs))

print("== Q1")
for K in KS:
    rows = {}
    surv = 0
    for p in plist:
        k = first[p]
        if k == 0 or k > K: surv += 1
        else: rows[k] = rows.get(k, 0) + 1
    print(f"K={K}: survivors {surv} ({surv/n:.4f}), nonempty rows {len(rows)}, exponents all survivors:",
          all(first[e] == 0 or first[e] > K for e in EXPS if e <= N))
    if K == 1000:
        print("rows 1..64:", {k: rows[k] for k in sorted(rows) if k <= 64})
        print("rows 65..256 total", sum(v for k, v in rows.items() if 65 <= k <= 256),
              " rows 257..1000 total", sum(v for k, v in rows.items() if k >= 257))

print("== Q2")
viol = 0
bykmod = [0] * 4; byq8 = [0] * 8
for p, k in pairs:
    q = 2 * k * p + 1
    bykmod[k % 4] += 1; byq8[q % 8] += 1
    if q % 8 not in (1, 7): viol += 1
    if k % 4 == 2: viol += 1
    if k % 4 == 1 and p % 4 != 3: viol += 1
    if k % 4 == 3 and p % 4 != 1: viol += 1
print("violations:", viol, " by k mod 4:", bykmod, " by q mod 8:", byq8)

print("== Q3")
row1 = {p for p in plist if first[p] == 1}
rule = {p for p in plist if p > 3 and p % 4 == 3 and sympy.isprime(2 * p + 1)}
print("row1", len(row1), "rule", len(rule), "equal", row1 == rule)

print("== Q4")
C2 = 1.0
for r in sympy.primerange(3, 10 ** 6):
    C2 *= 1 - 1 / (r - 1) ** 2
print("C2", C2)
def hl_factor(k):
    f = 1.0
    for r in sympy.primefactors(k):
        if r > 2: f *= (r - 1) / (r - 2)
    return f
HLF = [0.0] + [hl_factor(k) for k in range(1, 1001)]
hitset = set(pairs)
grp = [(1, 1), (3, 3), (4, 4), (5, 5), (7, 7), (8, 8), (9, 9), (11, 12), (13, 16), (17, 32), (33, 64), (65, 256), (257, 1000)]
gc = {g: [0, 0, 0.0, 0.0] for g in grp}   # admissible, hits, expected, variance
def gkey(k):
    for a, b in grp:
        if a <= k <= b: return (a, b)
    return None
bands = [(2, 10 ** 4), (10 ** 4, 10 ** 5), (10 ** 5, N + 1), (2, N + 1)]
acc = {(K, b): [0, 0, 0.0, 0.0, 0.0, 0.0] for K in KS for b in bands}  # n, obs, exact, var, hl, crude
cls = {(K, r): [0, 0, 0.0] for K in (64, 1000) for r in (1, 3)}
for p in plist:
    surv_e = {K: 1.0 for K in KS}
    surv_h = {K: 1.0 for K in KS}
    M = (1 << p) - 1
    for k in range(1, 1001):
        q = 2 * k * p + 1
        if q % 8 not in (1, 7) or q >= M:
            continue
        x = min(2 * C2 * HLF[k] / math.log(q) / k, 1 - 1e-12)
        prime = sympy.isprime(q)
        if prime:
            g = gkey(k)
            gc[g][0] += 1; gc[g][1] += (p, k) in hitset; gc[g][2] += 1 / k; gc[g][3] += (1 / k) * (1 - 1 / k)
        for K in KS:
            if k <= K:
                surv_h[K] *= 1 - x
                if prime: surv_e[K] *= 1 - 1 / k
    for K in KS:
        o = 1 if (first[p] == 0 or first[p] > K) else 0
        crude = math.log(2 * p) / math.log(2 * K * p)
        for b in bands:
            if b[0] <= p < b[1]:
                a = acc[(K, b)]
                a[0] += 1; a[1] += o; a[2] += surv_e[K]; a[3] += surv_e[K] * (1 - surv_e[K]); a[4] += surv_h[K]; a[5] += crude
        if K in (64, 1000) and p > 2:
            c = cls[(K, p % 4)]; c[0] += 1; c[1] += o; c[2] += surv_e[K]
for g in grp:
    a, h, e, v = gc[g]
    print(f"  k {g[0]}-{g[1]}: admissible {a} divides {h} expected {e:.1f} sd {math.sqrt(v):.1f}")
for K in KS:
    for b in bands:
        nb, o, e, v, h, c = acc[(K, b)]
        print(f"  K={K} [{b[0]},{b[1]}): n={nb} obs {o} ({o/nb:.4f}) exact-q {e:.1f} ({e/nb:.4f}) z={(o-e)/math.sqrt(v):+.2f} HL {h:.1f} ({h/nb:.4f}) crude {c:.1f} ({c/nb:.4f})")

print("== Q5")
for (K, r), (m, o, e) in sorted(cls.items()):
    print(f"  K={K} p={r} mod 4: n={m} survivors {o} frac {o/m:.4f} model {e/m:.4f}")
odd = [e for e in EXPS if 2 < e <= N]
c1 = sum(e % 4 == 1 for e in odd); c1all = sum(e % 4 == 1 for e in EXPS if e > 2)
tail = lambda m, x: sum(math.comb(m, i) for i in range(x, m + 1)) / 2 ** m
print(f"  odd exponents <= N: {len(odd)}, =1 mod 4: {c1}, P(X>={c1}) = {tail(len(odd), c1):.4f}; all 50: {c1all}, P = {tail(50, c1all):.4f}")
print("  mod 8:", {r: sum(e % 8 == r for e in EXPS if e > 2) for r in (1, 3, 5, 7)})
nbs = 0; tot = 0
for e in odd:
    i = plist.index(e)
    for j in (i - 1, i + 1):
        tot += 1; nbs += (first[plist[j]] == 0 or first[plist[j]] > 64)
print(f"  neighbours survivors (K=64): {nbs}/{tot}")

print("== Q6")
for K in (16, 64, 94, 128, 256, 1000):
    tot = 0; adm = 0
    for p in plist:
        t = first[p] if 0 < first[p] <= K else K
        tot += t
        adm += sum(1 for k in range(1, t + 1) if p == 2 or (2 * k * p + 1) % 8 in (1, 7))
    print(f"  K={K}: straightforward {tot}, mod-8 filtered {adm}")
