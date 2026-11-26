"""Referee check (b): phi-weighted average over squarefree N of M^(N) vs Zubrilina's M_k (her Thm 1 constants)."""
import math
import numpy as np
import sympy as sp
PR = list(sp.primerange(2, 2 * 10 ** 6))
la = np.array(PR, dtype=float)
# her constants, verbatim from the source (murmurationsv2.tex, lines 266-268)
alpha = 2 * math.pi * math.exp(np.sum(np.log((1 - la - 2 * la ** 2 + la ** 4) / (la ** 4 - 2 * la ** 2 + la))))
beta = 2 * math.pi * math.exp(np.sum(np.log((-1 + la ** 2 + la ** 3) / (la * (-1 + la + la ** 2)))))
gamma = 12 * math.exp(np.sum(np.log(la * (1 + la) / (-1 + la + la ** 2))))
def nuz(r): return math.prod(1 + l * l / (l ** 4 - 2 * l * l - l + 1.0) for l in sp.factorint(r)) if r > 1 else 1.0
def U(k, x): return 1.0 if k == 2 else (4 * x * x - 1 if k == 4 else None)
def M_zub(y, k):
    s = beta * math.sqrt(y) / (k - 1) - (gamma * y if k == 2 else 0)
    r = 1
    while r * r < 4 * y:
        s += (-1) ** (k // 2 - 1) * alpha * nuz(r) * U(k, r / (2 * math.sqrt(y))) * math.sqrt(4 * y - r * r) / (k - 1)
        r += 1
    return s
# own per-level local data
def u_(l): return 1.0 / ((l - 1) ** 2 * (l + 1))
def e_(l, r): return l * l / (l * l - 1.0) if r % l == 0 else 1.0 - u_(l)
def avg_phi(f):
    """lim sum^sq phi(N) prod_{q|N} f(q) / sum^sq phi(N): Euler product over l < 2e6"""
    lg = 0.0
    for l in PR:
        c = f(l) * (1 - 1.0 / l)
        lg += math.log1p((c - 1) / (l + 1)) - math.log1p(-1.0 / (l * (l + 1)))
    return math.exp(lg)
A_pi = avg_phi(lambda q: 1 + 1.0 / q)
A_kap = avg_phi(lambda q: q / (q - 1.0))
C0 = math.exp(sum(math.log1p(-u_(l)) for l in PR))
def nup(r):
    v = C0
    for l in sp.factorint(r): v *= e_(l, 0) / (1 - u_(l))
    return v
A_r = {r: nup(r) * avg_phi(lambda q, r=r: q / ((q - 1.0) * e_(q, r))) for r in range(1, 4)}
def M_avg(y, k):
    s = 2 * math.pi * A_pi * math.sqrt(y) / (k - 1) - (12 * A_kap * y if k == 2 else 0)
    r = 1
    while r * r < 4 * y:
        s += (-1) ** (k // 2 - 1) * (12 / math.pi) * A_r[r] * U(k, r / (2 * math.sqrt(y))) * math.sqrt(4 * y - r * r) / (k - 1)
        r += 1
    return s
print(f'beta = {beta:.6f} vs 2 pi avg(pi(N)) = {2*math.pi*A_pi:.6f};  gamma = {gamma:.6f} vs 12 avg(kappa) = {12*A_kap:.6f}')
for r in range(1, 4):
    print(f'  r={r}: alpha nu(r) = {alpha*nuz(r):.6f} vs (12/pi) avg(kappa nu_N(r)) = {12/math.pi*A_r[r]:.6f}')
# brute force over squarefree N <= X
X = 400000
spf = np.zeros(X + 1, dtype=np.int64)
for p in sp.primerange(2, X + 1):
    m = spf[p::p]; m[m == 0] = p
ys = [0.1, 0.2, 0.5, 1.0, 2.0]
num = {k: np.zeros(len(ys)) for k in (2, 4)}; den = 0.0
for N in range(2, X + 1):
    S = []; n = N; sq = True
    while n > 1:
        p = int(spf[n]); n //= p
        if n % p == 0: sq = False; break
        S.append(p)
    if not sq: continue
    phi = math.prod(q - 1 for q in S)
    kap = N / phi; piN = math.prod(1 + 1.0 / q for q in S)
    nu = {r: nup(r) / math.prod(e_(q, r) for q in S) for r in (1, 2)}
    den += phi
    for k in (2, 4):
        for i, y in enumerate(ys):
            s = 2 * math.pi * piN * math.sqrt(y) / (k - 1) - (12 * kap * y if k == 2 else 0)
            r = 1
            while r * r < 4 * y:
                s += (-1) ** (k // 2 - 1) * (12 / math.pi) * kap * nu[r] * U(k, r / (2 * math.sqrt(y))) * math.sqrt(4 * y - r * r) / (k - 1)
                r += 1
            num[k][i] += phi * s
for k in (2, 4):
    print(f'k={k}:  y   | Zubrilina | Euler avg of M^(N) | brute force N<=4e5')
    for i, y in enumerate(ys):
        print(f'      {y:4.1f} | {M_zub(y,k):9.5f} | {M_avg(y,k):9.5f} | {num[k][i]/den:9.5f}')
