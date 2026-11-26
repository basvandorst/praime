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

PPL = [2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32]
res = []
for k in range(K):
    h = H[k]
    h0 = mp.fsum(h[j] * m0[j] for j in range(n))
    Wh = Wv[k + 1]
    eta0 = lam / Wh
    eta = fac * eta0
    r = mp.sqrt(2 * eta / G2[k])
    q = [b[j] - eta * h[j] for j in range(n)]
    roots = mp.polyroots(list(reversed(q)), maxsteps=400, extraprec=300)
    realpos = [rt for rt in roots if abs(mp.im(rt)) < mp.mpf(10) ** -30 and mp.re(rt) >= 0]
    h2k = mp.fsum(h[j] * rows[k][2][j] for j in range(n))
    # eps_max from the other nodes: G''(u_j)/(-h''(u_j)) where h''(u_j) < 0
    epsn = mp.inf
    for jj in range(K):
        if jj == k: continue
        h2j = mp.fsum(h[j] * rows[jj][2][j] for j in range(n))
        if h2j < 0: epsn = min(epsn, G2[jj] / (-h2j))
    # B side: min over t of B_G/(-B_h) where B_h<0, t in [0,400] log/linear grid
    epsB = mp.inf
    for i2 in range(1, 4001):
        t2 = (mp.mpf(i2) / 10) ** 2
        BG = mp.fsum(b[j] * t2 ** j for j in range(n)); Bh = mp.fsum(h[j] * t2 ** j for j in range(n))
        if Bh < 0: epsB = min(epsB, BG / (-Bh))
    nk = PPL[k]
    res.append(dict(k=k, eta=eta, r=r, h0=h0, Wh=Wh, h2k=h2k, epsn=epsn, epsB=epsB, nk=nk))
    print('k=%2d e^u=%-14s n_k=%2d mu=%-11s 2L/sqrt n=%-9s W(h)=%-12s h(0)=%-10s eta0=%-9s G2=%-9s h2=%-9s r=%-9s trap e^u in [%s, %s] contains log n_k: %s  (u_k-log n_k=%s)  cond(i) roots: %d  eps_node=%s eps_B=%s' % (
        k + 1, mp.nstr(mp.exp(us[k]), 12), nk, mp.nstr(mu[k], 8), mp.nstr(2 * von_mangoldt(nk) / mp.sqrt(nk), 7),
        mp.nstr(Wh, 10), mp.nstr(h0, 5), mp.nstr(eta0, 4), mp.nstr(G2[k], 4), mp.nstr(h2k, 4), mp.nstr(r, 4),
        mp.nstr(mp.exp(us[k] - r), 8), mp.nstr(mp.exp(us[k] + r), 8), abs(us[k] - mp.log(nk)) < r,
        mp.nstr(us[k] - mp.log(nk), 4), len(realpos), mp.nstr(epsn, 4), mp.nstr(epsB, 4)), flush=True)

umax = float(us[-1]) + uext
nstep = int(umax / step)
pos = {k: [] for k in range(K)}
minratio = {k: mp.inf for k in range(K)}
epsgrid = {k: mp.inf for k in range(K)}
for i in range(1, nstep + 1):
    u = mp.mpf(i * step)
    vals = evalc([combG] + combH, u)
    Gv = vals[0]
    for k in range(K):
        w = res[k]['eta'] * vals[k + 1] - Gv
        if w > 0:
            pos[k].append(i * step)
        far = abs(u - us[k]) > 1.5 * res[k]['r'] + step
        if vals[k + 1] > 0 and far:
            ratio = Gv / (res[k]['eta'] * vals[k + 1])
            if ratio < minratio[k]: minratio[k] = ratio
        if vals[k + 1] < 0 and min(abs(u - x) for x in us) > step:
            epsgrid[k] = min(epsgrid[k], Gv / (-vals[k + 1]))


def runs(L):
    out = []
    for x in L:
        if out and abs(x - out[-1][1] - step) < 1e-9: out[-1][1] = x
        else: out.append([x, x])
    return [(round(a, 4), round(c, 4)) for a, c in out]


for k in range(K):
    R = res[k]
    epsmax = min(R['epsn'], R['epsB'], epsgrid[k])
    dl = abs(R['h2k']) * R['r'] ** 2 / 2
    lower = R['Wh'] * (1 - 1 / fac) / (1 + dl)
    upper = (R['Wh'] + lam / epsmax + lam * max(-R['h0'], 0)) / (1 - dl)
    print('k=%2d S on grid (components): %s ; trap (u) = [%s, %s] ; min_grid outside 1.5r of G/(eta h): %s ; eps_max~%s ; mass bounds for every d>=dp: %s <= mu_d(S_k), mu_d(I_k) <= %s' % (
        k + 1, runs(pos[k]), mp.nstr(us[k] - R['r'], 6), mp.nstr(us[k] + R['r'], 6), mp.nstr(minratio[k], 4),
        mp.nstr(epsmax, 4), mp.nstr(lower, 6), mp.nstr(upper, 6)), flush=True)
