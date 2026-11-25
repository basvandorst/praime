"""Referee check 5: c_L = prod_p G_p(1/p) (conditionally convergent) against B_L(1)/L(1, Ad pi) (absolutely convergent B),
own Euler products; L(1, Ad) from PARI:
  Sym^2 11a (d = 3): Ad = Sym^4 E + Sym^2 E (unramified p); L(1,Ad) = L(Sym^4 E, 3) L(Sym^2 E, 2) (arithmetic s);
     at p = 11 B uses PARI's own local factors (lfuneuler), so the bad-prime convention cancels exactly.
  Delta (d = 2): Ad = Sym^2 Delta, built with lfuncreate from L(s,Sym^2 f) = zeta(2s) sum lambda(n^2) n^-s, mu = (1, 11, 12)."""
import math, cmath, time
import numpy as np
import cypari2
P = cypari2.Pari(); P.allocatemem(800 * 10**6); P.set_real_precision(30)

def primes(N):
    s = np.ones(N + 1, bool); s[:2] = False
    for i in range(2, int(N ** .5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]

PMAX = 1000000
pr = primes(PMAX)
# ---------------- Sym^2 11a
E = P.ellinit([0, -1, 1, -10, -20])
t0 = time.time()
an = np.array([int(x) for x in P.ellan(E, PMAX)], dtype=float)
L2 = P.lfunsympow(E, 2); L4 = P.lfunsympow(E, 4)
v2 = float(P.lfun(L2, 2)); v4 = float(P.lfun(L4, 3))
print(f'L(Sym^2 11a, 2) = {v2:.12f}, L(Sym^4 11a, 3) = {v4:.12f}  (arithmetic; = analytic s = 1); checkfeq {P.lfuncheckfeq(L2)}, {P.lfuncheckfeq(L4)}; {time.time()-t0:.1f}s')
eu2 = P.lfuneuler(L2, 11); eu4 = P.lfuneuler(L4, 11)
A11_2 = float(P.substpol(eu2, P('x'), P('11^-2'))) if False else float(P(f'subst({eu2}, x, 11^-2)'))
A11_4 = float(P(f'subst({eu4}, x, 11^-3)'))
print(f'PARI local factors at 11: Sym^2 {eu2}, Sym^4 {eu4}; values at analytic s=1: {A11_2:.12f}, {A11_4:.12f}')
logG = 0j; logB = 0.0; marks = {}
Bdev = []
for p in pr:
    p = int(p)
    if p == 11:
        # pi = Sym^2 11a at 11: lambda(11^k) = 11^-k, mu_L(11) = -1/11: F_11(X) = 1 - X/121
        F = 1 - (1 / 11) / 121
        Gp = F / (1 - 1 / 11)
        logG += math.log(Gp)
        logB += math.log(Gp * A11_2 * A11_4)
        continue
    x = an[p - 1] / math.sqrt(p)
    A = x * x - 1                       # lambda(p) (real), omega = 1
    l3 = A ** 3 - 2 * A * A + 1         # lambda(p^3) = h_3 with e1 = A, e2 = A, e3 = 1
    Gp = 1 - (A * A - 1) / p + l3 / p ** 2     # claim (iii), d = 3
    # Satake of E: beta = e^{i phi}, 2 cos phi = x
    c = max(-1.0, min(1.0, x / 2)); ph = math.acos(c)
    L2p = 1.0; L4p = 1.0
    for j in (-1, 0, 1):
        L2p *= 1 / abs(1 - cmath.exp(2j * j * ph) / p) if j == 0 else 1
    L2p = 1 / ((1 - 1 / p) * abs(1 - cmath.exp(2j * ph) / p) ** 2)
    L4p = 1 / ((1 - 1 / p) * abs(1 - cmath.exp(2j * ph) / p) ** 2 * abs(1 - cmath.exp(4j * ph) / p) ** 2)
    Bp = Gp * L2p * L4p
    Bdev.append((p, abs(Bp - 1) * p * p))
    logG += math.log(Gp); logB += math.log(Bp)
    if p in (997, 9973, 99991, 999983): marks[p] = (math.exp(logG.real), math.exp(logB))
cS = None
for p, (g, b) in sorted(marks.items()):
    cS = b / (v2 * v4)
    print(f'Sym2 11a P={p}: prod G_p(1/p) = {g:.6f};  B_(P)(1)/L(1,Ad) = {cS:.7f}')
print(f'   max_p |B_p - 1| p^2 = {max(v for _, v in Bdev):.2f} (B_p = 1 + O(p^-2) at unramified p)')
print(f'   note: c = 1.94470 (constants file), partial product 1.94469 at 10^6')

# ---------------- Delta
t0 = time.time()
PM = 300000
tau = P(f'vector({PM}, n, ramanujantau(n))')
tl = np.array([float(tau[i]) for i in range(PM)])  # tau(1..PM)
P('lamS(n) = sumdiv(n, m, if(issquare(m), my(k = n/m); ramanujantau(k^2)/k^11, 0))')
LS = P.lfuncreate(P('[n -> vector(n, j, lamS(j)), 0, [1, 11, 12], 1, 1, 1]'))
cf = P.lfuncheckfeq(LS)
vS = float(P.lfun(LS, 1))
print(f'L(1, Sym^2 Delta) = {vS:.12f}; lfuncheckfeq = {cf}; {time.time()-t0:.1f}s')
logG = 0.0; logB = 0.0; marks = {}
for p in primes(PM):
    p = int(p)
    lp = tl[p - 1] / p ** 5.5
    lp2 = lp * lp - 1
    Gp = 1 - lp2 / p                              # claim (iii), d = 2, omega = 1
    c = max(-1.0, min(1.0, lp / 2)); ph = math.acos(c)
    L2p = 1 / ((1 - 1 / p) * abs(1 - cmath.exp(2j * ph) / p) ** 2)
    logG += math.log(Gp); logB += math.log(Gp * L2p)
    if p in (997, 9973, 99991, 299993): marks[p] = (math.exp(logG), math.exp(logB))
for p, (g, b) in sorted(marks.items()):
    print(f'Delta P={p}: prod G_p(1/p) = {g:.6f};  B_(P)(1)/L(1,Sym^2) = {b / vS:.7f}')
print('   note: c = 2.36082 (constants file), partial product 2.36059 at 3*10^5')
