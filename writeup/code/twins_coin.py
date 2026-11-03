import sys
import numpy as np
from math import log, exp, sqrt

GAMMA = 0.5772156649015329
C2 = 0.6601618158468696
X = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**9
W = int(float(sys.argv[2])) if len(sys.argv) > 2 else 10**8
H = 1e-4
UMAX = 40.0


def buchstab(kappa, K):
    n = int(round((UMAX - 1.0) / H)) + 1
    u = 1.0 + H * np.arange(n)
    G = np.empty(n)
    d = int(round(1.0 / H))
    lo = u <= 2.0 + 1e-12
    G[lo] = K
    for i in range(lo.sum(), n):
        ui = u[i]
        f_prev = kappa * (ui - H) ** (kappa - 1) * G[i - 1 - d] / (ui - H - 1.0) ** kappa
        f_now = kappa * ui ** (kappa - 1) * G[i - d] / (ui - 1.0) ** kappa
        G[i] = G[i - 1] + 0.5 * H * (f_prev + f_now)
    return u, G / u ** kappa


def value_and_slope(arr, u):
    i = int(round((min(u, UMAX - 0.1) - 1.0) / H))
    return arr[i], (arr[i + 1] - arr[i - 1]) / (2 * H)


def small_primes(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.flatnonzero(s)


U, B1 = buchstab(1, exp(GAMMA))
_, B2naive = buchstab(2, exp(2 * GAMMA))
print(f"Buchstab limits: e^g omega -> {B1[-1]:.4f}; naive two-dimensional equation (u^2 B)' = 2u B(u-1) -> {B2naive[-1]:.4f} (fundamental lemma needs 1)")

top = int(sqrt(X + W + 2)) + 1
primes = small_primes(top)
comp = np.zeros(W + 2, dtype=bool)
for p in primes:
    comp[(-X) % p::p] = True
P = ~comp
lnX = log(X + W / 2)
nP = int(P[:W].sum())
nTw = int((P[:W] & P[2:W + 2]).sum())
print(f"X={X:.0e} W={W:.0e}: primes {nP} (ratio to W/ln X {nP*lnX/W:.4f}), twins {nTw} (ratio to 2C2 W/ln^2 X {nTw*lnX**2/(2*C2*W):.4f})")

tail = {}
t = 1.0
for p in primes[::-1]:
    if p > 2:
        t *= 1 - 1 / (p - 1) ** 2
    tail[p] = t
tail_top = 1 - 1 / (top * log(top))

ys = [7, 11, 23, 53, 101, 211, 401, 797, 1601]
y04 = int(round(X**0.4))
ys = sorted(set(y for y in ys + [y04, 2 * y04, 4 * y04] if y <= top)) + [top]
comp = np.zeros(W + 2, dtype=bool)
V1, V2, j = 1.0, 0.5, 0
print()
print("     y      u  surv/(W V1 B1)  twinsurv/(W V2 B1emp^2)   rho    rho2  rho2/rho^2  C2 tail   c/(ab)    sd")
for y in ys:
    while j < len(primes) and primes[j] <= y:
        p = primes[j]
        comp[(-X) % p::p] = True
        V1 *= 1 - 1 / p
        if p > 2:
            V2 *= 1 - 2 / p
        j += 1
    S = ~comp
    idx = np.flatnonzero(S[:W] & S[2:W + 2])
    nS = int(S[:W].sum())
    nT = idx.size
    u = lnX / log(y)
    b1, d1 = value_and_slope(B1, u)
    b1emp = nS / (W * V1)
    a = P[idx].mean()
    b = P[idx + 2].mean()
    both = P[idx] & P[idx + 2]
    c = both.mean()
    rho = nP / nS
    rho2 = nTw / nT
    sd = sqrt((1 - c) / (c * nT))
    nxt = primes[j] if j < len(primes) else None
    ct = (tail[nxt] if nxt is not None else 1.0) * tail_top
    print(f"{y:8d} {u:6.3f}     {nS/(W*V1*(b1+d1/log(y))):.4f}          {nT/(W*V2*b1emp**2):.5f}          {rho:.4f} {rho2:.4f}   {rho2/rho**2:.4f}    {ct:.4f}   {c/(a*b):.4f}  {sd:.4f}")
    if y in (y04, 2 * y04):
        pn, pm = P[idx], P[idx + 2]
        cells = [int((pn & pm).sum()), int((pn & ~pm).sum()), int((~pn & pm).sum()), int((~pn & ~pm).sum())]
        ea, eb = pn.mean(), pm.mean()
        pred = [ea * eb, ea * (1 - eb), (1 - ea) * eb, (1 - ea) * (1 - eb)]
        kind = "P2 exactly" if y**3 > X + W + 2 else "P2 or P3"
        print(f"         census at y={y} ({kind}), {nT} twin survivors: (p,p) (p,P2) (P2,p) (P2,P2) = {cells}; independent coins: {[round(q*nT) for q in pred]}")
