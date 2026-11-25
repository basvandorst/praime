"""Step 1: own pair covariances C_ij(w) for the Fn sets, own completeness check, own Kronecker table,
own diagonal-law predictions (T -> inf prime sum to 2e5, finite-range part by direct synthesis for n <= 300)."""
import sys, math, time
import numpy as np
import mpmath as mp
sys.path.insert(0, __import__('os').path.dirname(__file__))
from mylib import *

t0 = time.time()
# theta check against mpmath
for q, a, t in [(1001, 0, 100.0), (1003, 1, 333.3), (1300, 0, 599.9), (1007, 1, 101.7)]:
    ex = float(mp.mpf(t) / 2 * mp.log(mp.mpf(q) / mp.pi) + mp.im(mp.loggamma(mp.mpf(0.25) + mp.mpf(a) / 2 + 0.5j * mp.mpf(t))))
    print('theta check q=%d a=%d t=%.1f: diff %.2e' % (q, a, t, theta(np.array([t]), q, a)[0] - ex))

Z = load_zeros(); Ds = sorted(Z.keys())
print('sets', len(Ds))
# completeness: N(t) - (theta(t) - theta(0))/pi on 200 heights; plus nearest-integer offset; largest unfolded gap
bad = []
for D in Ds:
    g = Z[D]; q = abs(D); a = 0 if D > 0 else 1
    t = np.linspace(T1, T2, 200)
    th0 = float(mp.im(mp.loggamma(0.25 + a / 2)))  # = 0 (real argument)
    s = np.searchsorted(g, t) - (theta(t, q, a) - th0) / math.pi
    b = g[(g >= T1) & (g <= T2)]; u = b / (2 * math.pi) * np.log(q * b / (2 * math.pi * math.e))
    if np.abs(s).max() > 1.5 or np.diff(u).max() > 3.0:
        bad.append((D, round(float(np.abs(s).max()), 3), round(float(np.diff(u).max()), 3), round(float(s.mean()), 3)))
print('flagged by the round-138 rule:', bad)
# mean of s over sets (smooth-count constant)
Ds = [D for D in Ds if D not in [b[0] for b in bad]]
n = len(Ds); print('kept', n, 'pairs', n * (n - 1) // 2)

C = np.zeros((16, n, n)); MU = np.zeros((16, n)); V = np.zeros((16, n))
for k, w in enumerate(W):
    tg = grid(w)
    S = np.array([dn(Z[D], D, w, tg) for D in Ds])
    C[k] = S @ S.T / len(tg); MU[k] = S.mean(1); V[k] = S.var(1)
print('covariances done %.0fs' % (time.time() - t0))

P, X = kron_pari(Ds, 200000)
print('kronecker done %.0fs, primes %d' % (time.time() - t0, len(P)))
# check against the frozen file's own psi columns later
# per-prime A_p(w), B_p(w) (T -> inf) for all primes
lp = np.log(P.astype(float))
A = np.zeros((16, len(P))); B = np.zeros((16, len(P)))
for m in range(1, 60):
    c = 1.0 / (m * m * P.astype(float) ** m)
    if c.max() < 1e-16: break
    t = c[None, :] * np.sin(m * W[:, None] * lp[None, :] / 2) ** 2
    if m % 2: A += t
    else: B += t
A *= 2 / math.pi ** 2; B *= 2 / math.pi ** 2
Xf = X.astype(float); X2 = Xf * Xf
CINF = np.array([(Xf * A[k]) @ Xf.T + (X2 * B[k]) @ X2.T for k in range(16)])
REST = np.array([(Xf[:, 4:] * A[k, 4:]) @ Xf[:, 4:].T + (X2[:, 4:] * B[k, 4:]) @ X2[:, 4:].T for k in range(16)])
# finite-range part, direct synthesis for prime powers n <= 300
PW = []
for j, p in enumerate(P[P <= 300]):
    m = 1
    while p ** m <= 300: PW.append((int(p) ** m, m, j)); m += 1
lnn = np.array([math.log(a) for a, _, _ in PW]); mm = np.array([b for _, b, _ in PW]); jj = np.array([c for _, _, c in PW])
CH = Xf[:, jj] ** mm[None, :]
DEL = np.zeros_like(CINF); ATAB = np.zeros((16, 4)); BTAB = np.zeros((16, 4))
for k, w in enumerate(W):
    tg = grid(w)
    cn = -(2 / math.pi) / (mm * np.exp(lnn / 2)) * np.sin(w * lnn / 2)
    Ssyn = (CH * cn[None, :]) @ np.cos(np.outer(lnn, tg + w / 2))      # (n, M)
    Cm = Ssyn @ Ssyn.T / len(tg)
    DEL[k] = Cm - (CH * (cn ** 2 / 2)[None, :]) @ CH.T
    for a_, p in enumerate((2, 3, 5, 7)):
        Ap = Bp = 0.0; m = 1
        while 1 / (m * m * p ** m) > 1e-15:
            base = 2 / math.pi ** 2 / (m * m * p ** m) * math.sin(m * w * math.log(p) / 2) ** 2 * (1 + np.mean(np.cos(2 * m * (tg + w / 2) * math.log(p))))
            if m % 2: Ap += base
            else: Bp += base
            m += 1
        ATAB[k, a_] = Ap; BTAB[k, a_] = Bp
FULL = CINF + DEL
np.savez(OUT + '/s1.npz', Ds=np.array(Ds), C=C, MU=MU, V=V, CINF=CINF, REST=REST, DEL=DEL, FULL=FULL, ATAB=ATAB, BTAB=BTAB,
         X7=X[:, :4], A=A[:, :200], B=B[:, :200])
np.save(OUT + '/kron.npy', X); np.save(OUT + '/primes.npy', P)
print('done %.0fs' % (time.time() - t0))
