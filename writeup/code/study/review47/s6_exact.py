"""frozen J (Euler product to 3e4) vs exact: H(u) = |L(1 - iu, psi)|^2 A(u), L from PARI lfun of the primitive psi*
(Kronecker symbol of the fundamental discriminant attached to D1 D2) with the removed Euler factors at p | q1 q2,
A = prod_{p not | q1q2} (1 - y)(1 - 2/p + y)/(1 - 1/p)^2 to p < 2e6 (y = psi(p) p^{iu}/p); own quadrature (step 0.02, |u| <= 7)."""
import sys, math, time
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import numpy as np
from rlib import *
import cypari2

pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6)
P = np.load(FZ + 'pred_OO.npz'); Ds = [int(x) for x in P['Ds']]
rng = np.random.default_rng(11)
ks = list(rng.choice(len(P['I']), 5, replace=False))
J5 = P['J_0.5']; ks += [int(np.argmin(J5)), int(np.argmax(J5))]
# one pair with gcd > 1
for k in rng.permutation(len(P['I'])):
    if math.gcd(abs(Ds[P['I'][k]]), abs(Ds[P['J'][k]])) > 1: ks.append(int(k)); break
pr = primes_upto(400000)
u = np.arange(-7.0, 7.0 + 1e-9, 0.02)
t0 = time.time()
for k in ks:
    D1, D2 = Ds[P['I'][k]], Ds[P['J'][k]]
    q1, q2 = abs(D1), abs(D2)
    Dstar = int(pari('coredisc(%d)' % (D1 * D2)))
    pari('L = lfuninit(lfuncreate(%d), [1, 0, 7.5])' % Dstar)
    Lv = np.array([complex(pari('lfun(L, 1 - I*(%.6f))' % x)) for x in u])
    bad = [p for p in factor(q1 * q2) if Dstar % p != 0]
    for p in bad:
        Lv = Lv * (1 - kron(Dstar, p) * np.exp((-1 + 1j * u) * math.log(p)))
    keep = (q1 * q2) % pr != 0
    ps = pr[keep].astype(float)
    # Kronecker for all primes up to 2e6 via tables mod q1, q2
    t1, t2 = real_table(D1), real_table(D2)
    psi = t1[pr[keep] % q1] * t2[pr[keep] % q2]
    lnA = np.zeros(len(u), complex)
    for j in range(0, len(ps), 4000):
        pp = ps[j:j + 4000]; ss = psi[j:j + 4000]
        y = ss[None, :] * np.exp(1j * u[:, None] * np.log(pp)[None, :]) / pp[None, :]
        lnA += np.sum(np.log((1 - y) * (1 - 2 / pp[None, :] + y) / (1 - 1 / pp[None, :]) ** 2), axis=1)
    H = np.abs(Lv) ** 2 * np.exp(lnA)
    out = []
    for s in (0.5, 1.0):
        w = np.exp(-u ** 2 / (2 * s * s)) * 0.02; w[0] /= 2; w[-1] /= 2
        Jx = np.sum(w * H)
        out.append('sigma %.1f: exact %.6f%+.1ei frozen %.6f diff %.1e' % (s, Jx.real, Jx.imag, P['J_%.1f' % s][k], Jx.real - P['J_%.1f' % s][k]))
    print('(%d, %d) D* = %d: %s   (%.0fs)' % (D1, D2, Dstar, ' | '.join(out), time.time() - t0), flush=True)
