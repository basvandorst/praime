"""Trap lemma, numerical evaluation (NUMERICAL, not a proof).
For the class-A optimum G of half-degree dp (polish_A_d<dp>.json): Hermite cardinal functions h_k
(h(u_j)=h'(u_j)=0 for j!=k, h(u_k)=1, h'(u_k)=0, top coefficient 0), W(h_k) from the prime sum,
eta0 = lam/W(h_k) (condition (ii) reduces to eta W(h) > lam since w(0) <= 0 under (i)),
and the set S = {u>0 : eta h_k > G} on a grid, for eta = fac*eta0.
usage: python trap.py dp dps fac step umax_extra
"""
import sys, json
import os
CERTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'convex', 'certs')
sys.path.insert(0, CERTS)
import mpmath as mp
from core import basis_polys, divisor_count, von_mangoldt

dp = int(sys.argv[1]); mp.mp.dps = int(sys.argv[2])
fac = mp.mpf(sys.argv[3]) if len(sys.argv) > 3 else mp.mpf(2)
step = float(sys.argv[4]) if len(sys.argv) > 4 else 0.005
uext = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
P = json.load(open(CERTS + '/polish_A_d%d.json' % dp))
b = [mp.mpf(x) for x in P['b']]
us = [mp.mpf(x) for x in P['touch_points']]
mu = [mp.mpf(x) for x in P['mu']]
lam = mp.mpf(P['lam'])
K = len(us); n = dp + 1
fr = basis_polys(dp, 2)
tomp = lambda L: [mp.mpf(a.numerator) / a.denominator for a in L]
polys = [[(tomp(P0), tomp(P1)) for (P0, P1) in lst] for lst in fr]
maxdeg = max(max(len(P0), len(P1)) for lst in polys for (P0, P1) in lst)


def comb(a, r):
    """coefficient lists of sum_j a_j D^r g_j."""
    Q0 = [mp.mpf(0)] * maxdeg; Q1 = [mp.mpf(0)] * maxdeg
    for j in range(n):
        P0, P1 = polys[j][r]
        for i, c in enumerate(P0): Q0[i] += a[j] * c
        for i, c in enumerate(P1): Q1[i] += a[j] * c
    return Q0, Q1


def horner(Q, X):
    s = mp.mpf(0)
    for c in reversed(Q): s = s * X + c
    return s


def evalc(combs, u):
    """values of several combined functions at u (list of (Q0,Q1))."""
    u = mp.mpf(u); eu = mp.exp(u); X1 = mp.pi * eu
    lim = (mp.mp.dps + 20) * mp.log(10)
    out = [mp.mpf(0)] * len(combs)
    N = 1
    while True:
        X = mp.pi * N * eu
        if N >= 2 and 2 * (X - X1) - (maxdeg + 2) * mp.log(X / X1) - mp.log(N) > lim: break
        k0 = mp.besselk(0, 2 * X); k1 = mp.besselk(1, 2 * X)
        pref = 4 * divisor_count(N) / mp.sqrt(mp.pi * N) * mp.sqrt(X)
        for i, (Q0, Q1) in enumerate(combs):
            out[i] += pref * (horner(Q0, X) * k0 + horner(Q1, X) * k1)
        N += 1
    return out


# individual basis rows at the nodes (r = 0, 1, 2) and at 0
E = [[mp.mpf(1) if i == j else mp.mpf(0) for i in range(n)] for j in range(n)]
rows = {}
for k, u in enumerate(us):
    rows[k] = [evalc([comb(E[j], r) for j in range(n)], u) for r in range(3)]  # rows[k][r][j]
m0 = evalc([comb(E[j], 0) for j in range(n)], mp.mpf(0))
# Hermite cardinal functions
A = mp.matrix(n, n)
for k in range(K):
    for j in range(n):
        A[2 * k, j] = rows[k][0][j]; A[2 * k + 1, j] = rows[k][1][j]
A[n - 1, n - 1] = 1
H = []
for k in range(K):
    rhs = mp.matrix(n, 1); rhs[2 * k] = 1
    H.append([x for x in mp.lu_solve(A, rhs)])
G2 = [mp.fsum(b[j] * rows[k][2][j] for j in range(n)) for k in range(K)]
G0 = [mp.fsum(b[j] * rows[k][0][j] for j in range(n)) for k in range(K)]
G1 = [mp.fsum(b[j] * rows[k][1][j] for j in range(n)) for k in range(K)]
combG = comb(b, 0)
combH = [comb(h, 0) for h in H]
# prime-power sums
pp = [m for m in range(2, 90) if von_mangoldt(m) != 0]
Wv = [mp.mpf(0)] * (K + 1)
for m in pp:
    vals = evalc([combG] + combH, mp.log(m))
    w = 2 * von_mangoldt(m) / mp.sqrt(m)
    for i in range(K + 1): Wv[i] += w * vals[i]
g0val = mp.fsum(b[j] * m0[j] for j in range(n))
print('dp=%d K=%d dps=%d  G(0)=%s  W(G)=%s lam=%s' % (dp, K, mp.mp.dps, mp.nstr(g0val, 15), mp.nstr(Wv[0], 15), mp.nstr(lam, 15)))
print('node residuals max |G(u_k)|, |G\'(u_k)|:', mp.nstr(max(abs(x) for x in G0), 3), mp.nstr(max(abs(x) for x in G1), 3))
res = []
for k in range(K):
    h = H[k]
    h0 = mp.fsum(h[j] * m0[j] for j in range(n))
    Wh = Wv[k + 1]
    kkt = lam * h0 + mu[k]
    eta0 = lam / Wh
    eta = fac * eta0
    r = mp.sqrt(2 * eta / G2[k])
    # condition (i): B_G - eta B_h >= 0 on R: polynomial in s = t^2
    q = [b[j] - eta * h[j] for j in range(n)]
    roots = mp.polyroots(list(reversed(q)), maxsteps=400, extraprec=300) if n > 1 else []
    realpos = [rt for rt in roots if abs(mp.im(rt)) < mp.mpf(10) ** -30 and mp.re(rt) >= 0]
    hmax_coef = max(abs(x) for x in h)
    res.append(dict(k=k, eta=eta, r=r, h0=h0))
    print('k=%2d e^u=%-14s mu=%-12s W(h)=%-14s lam h0+mu=%-14s h(0)=%-11s eta0=%-10s G\'\'=%-10s r=%-10s r/u=%s  cond(i) real roots s>=0: %d  max|h_j|=%s' % (
        k + 1, mp.nstr(mp.exp(us[k]), 12), mp.nstr(mu[k], 10), mp.nstr(Wh, 12), mp.nstr(kkt, 12), mp.nstr(h0, 5),
        mp.nstr(eta0, 4), mp.nstr(G2[k], 4), mp.nstr(r, 4), mp.nstr(r, 3), len(realpos), mp.nstr(hmax_coef, 3)), flush=True)

# grid scan of S = {eta h_k > G}
umax = float(us[-1]) + uext
nstep = int(umax / step)
pos = {k: [] for k in range(K)}
minratio = {k: mp.inf for k in range(K)}
for i in range(1, nstep + 1):
    u = mp.mpf(i * step)
    vals = evalc([combG] + combH, u)
    Gv = vals[0]
    for k in range(K):
        w = res[k]['eta'] * vals[k + 1] - Gv
        if w > 0:
            pos[k].append(float(u))
        if vals[k + 1] > 0 and abs(u - us[k]) > 0.05:
            ratio = Gv / (res[k]['eta'] * vals[k + 1])
            if ratio < minratio[k]: minratio[k] = ratio
for k in range(K):
    print('k=%2d grid points of S away from u_k: %s ; min over grid (|u-u_k|>0.05) of G/(eta h) where h>0: %s' % (
        k + 1, (pos[k][:6], len(pos[k])), mp.nstr(minratio[k], 4)))
