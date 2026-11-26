"""Re-polish the class-A KKT optimum at half-degree d (deg B = 2d = 4K) at ~100 digits, and build the trap data.

Same equations as certs/polish.py (independent code: arb + own K0/K1 instead of mpmath.besselk):
  A b = e_0, rows of A: m (moments, so g(0) = 1), g(u_k), g'(u_k) (k = 1..K);  y = A^{-T} c;
  lambda = y_0, mu_k = y_{2k-1}, z_k = y_{2k};  Newton on z(u) = 0 (central differences, h = 2^-150).
Then for every node k: h_k = element of V_d with h(u_j) = h'(u_j) = 0 (j != k), h(u_k) = 1, h'(u_k) = 0,
top coefficient 0 (A with row 0 replaced by e_d);  W(h_k) = c.b_h, h_k(0) = m.b_h, h_k''(u_j), g''(u_k).
Start: touch points of certs/polish_A_d<d>.json (30 digits).
Output: hp_d<d>.json (strings, 100 digits).   Usage: python polish_hp.py d"""
import sys, json, time, os
from fam import *

d = int(sys.argv[1])
K = d // 2
n = d + 1
t0 = time.time()
BV = json.load(open("basevals.json"))
m = [arb(s) for s in BV["m"][:n]]
c = [arb(s) for s in BV["c"][:n]]
ST = json.load(open(os.path.join(common.CERTS, "polish_A_d%d.json" % d)))
us = [arb(s) for s in ST["touch_points"]]
F = Fam(d, 3)
_cache = {}


def rows_at(u, R=1):
    key = (u.mid().str(200, radius=False), R)
    if key not in _cache:
        r = F.rows(u.exp(), R=R)
        _cache[key] = r
    return _cache[key]


def system(us):
    A = [m]
    for u in us:
        r = rows_at(u)
        A.append([r[j][0] for j in range(n)])
        A.append([r[j][1] for j in range(n)])
    b = solve(A, [arb(1)] + [arb(0)] * (n - 1))
    AT = [[A[i][j] for i in range(n)] for j in range(n)]
    y = solve(AT, c)
    return A, b, y


H = arb(2) ** -150
for it in range(12):
    A, b, y = system(us)
    z = [y[2 + 2 * k] for k in range(K)]
    lam = y[0]
    zn = max(abs(float(x.mid().str(5, radius=False))) for x in z)
    print("it %d  lambda=%s  |z|=%.3e  (%.0fs)" % (it, lam.str(20, radius=False), zn, time.time() - t0), flush=True)
    if it > 0 and zn < 1e-95 * abs(float(lam.mid().str(5, radius=False))) + 1e-200:
        break
    Jm = []
    for k in range(K):
        up = list(us); up[k] = mid(up[k] + H)
        um = list(us); um[k] = mid(um[k] - H)
        zp = system(up)[2]; zm = system(um)[2]
        Jm.append([(zp[2 + 2 * i] - zm[2 + 2 * i]) / (2 * H) for i in range(K)])
    JT = [[Jm[k][i] for k in range(K)] for i in range(K)]       # J[i][k] = dz_i/du_k
    step = solve(JT, [-x for x in z])
    mx = max(abs(float(s.mid().str(5, radius=False))) for s in step)
    if mx > 0.02:
        step = [s * (0.02 / mx) for s in step]
    us = [mid(us[k] + step[k]) for k in range(K)]
    if mx < 1e-100:
        break
    _cache = {kk: v for kk, v in _cache.items() if any(kk[0] == u.mid().str(200, radius=False) for u in us)}

A, b, y = system(us)
lam = y[0]
mu = [y[1 + 2 * k] for k in range(K)]
z = [y[2 + 2 * k] for k in range(K)]
cost = dot(b, c)
r2 = [rows_at(u, R=3) for u in us]
gpp = [dot(b, [r2[k][j][2] for j in range(n)]) for k in range(K)]
gppp = [dot(b, [r2[k][j][3] for j in range(n)]) for k in range(K)]
gval = [dot(b, [r2[k][j][0] for j in range(n)]) for k in range(K)]
g1val = [dot(b, [r2[k][j][1] for j in range(n)]) for k in range(K)]
# h_k
H_out = []
for k in range(K):
    Ah = [[arb(0)] * d + [arb(1)]] + A[1:]
    rhs = [arb(0)] * n
    rhs[1 + 2 * k] = arb(1)
    bh = solve(Ah, rhs)
    W = dot(bh, c)
    h0 = dot(bh, m)
    hpp = [dot(bh, [r2[j][jj][2] for jj in range(n)]) for j in range(K)]
    hval = [dot(bh, [r2[j][jj][0] for jj in range(n)]) for j in range(K)]
    H_out.append(dict(k=k + 1, b=[astr(x, 100) for x in bh], W=astr(W, 60), h0=astr(h0, 60),
                      kkt_resid=astr(W - lam * h0 - mu[k], 5),
                      hpp=[astr(x, 30) for x in hpp], hval_check=[astr(x, 5) for x in hval]))
out = dict(d=d, deg=2 * d, K=K, lam=astr(lam, 100), cost=astr(cost, 100), mu=[astr(x, 100) for x in mu],
           z=[astr(x, 5) for x in z], u=[astr(x, 110) for x in us], b=[astr(x, 110) for x in b],
           gpp=[astr(x, 60) for x in gpp], gppp=[astr(x, 30) for x in gppp],
           g_at_u=[astr(x, 5) for x in gval], g1_at_u=[astr(x, 5) for x in g1val],
           lam_acc_bits=lam.rel_accuracy_bits(), b_acc_bits=min(x.rel_accuracy_bits() for x in b),
           stored_lam=ST["lam"], stored_mu=ST["mu"], h=H_out, time_s=time.time() - t0)
json.dump(out, open("hp_d%d.json" % d, "w"), indent=1)
print("done d=%d lam=%s  max|z|=%s  time %.0fs" % (d, lam.str(30, radius=False), max(abs(float(x.mid().str(5, radius=False))) for x in z), time.time() - t0))
