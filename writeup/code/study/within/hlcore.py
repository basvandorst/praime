"""Own Hardy-Littlewood machinery (no repo imports).
- S_tail(h, y): the pair singular series without the primes <= y,  S_{>y}(h) = C_{>y} prod_{p>y, p|h} (p-1)/(p-2),
  C_{>y} = prod_{p>y} (1 - 1/(p-1)^2).  (HL pair factor; mean 1 over every residue class mod prod_{p<=y} p.)
- neighbour laws per fine class u mod Q = prod_{p<=y} p in the order-1 conditioned expansion: given a prime at x = u
  (mod Q), each rough site x + t is prime with chance rho0 * S_{>y}(t), rho0 = (Q/phi(Q))/lam, independently.
"""
import numpy as np, math
def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]
_PR = primes_upto(10_000_000)
def C_tail(y):
    p = _PR[_PR > y].astype(float)
    s = np.sum(np.log1p(-1.0 / (p - 1) ** 2)); X = 1e7
    s -= 1.0 / (X * math.log(X))                  # tail sum_{p > 1e7} 1/p^2
    return math.exp(s)
def S_tail_table(hmax, y):
    """S_{>y}(h) for 0 <= h <= hmax (h = 0 entry set to nan)."""
    f = np.ones(hmax + 1)
    for p in _PR[(_PR > y) & (_PR <= hmax)].tolist():
        f[p::p] *= (p - 1) / (p - 2)
    f *= C_tail(y); f[0] = np.nan
    return f
def primorial(y):
    return int(np.prod([p for p in _PR[_PR <= y].tolist()]))
def rough_classes(Q):
    u = np.arange(Q); ok = np.ones(Q, bool)
    for p in _PR[_PR <= 50].tolist():
        if Q % p == 0: ok &= (u % p != 0)
    return np.nonzero(ok)[0]
def neighbour_class_laws(y, lam, qs, T=800, chunk=4096):
    """For each rough u mod Q: P_R[q][u_index, b] = P(next prime after x is = b mod q | x prime, x = u mod Q) and
    P_L likewise for the previous prime (b = residue mod q, full residue index 0..q-1).  Returns (Q, U, PR, PL, rho0)."""
    Q = primorial(y); U = rough_classes(Q); phiQ = len(U); rho0 = (Q / phiQ) / lam
    St = S_tail_table(T, y)[1:]                         # t = 1..T
    t = np.arange(1, T + 1)
    isr = np.zeros(Q, bool); isr[U] = True
    PR = {q: np.zeros((phiQ, q)) for q in qs}; PL = {q: np.zeros((phiQ, q)) for q in qs}
    for s0 in range(0, phiQ, chunk):
        uu = U[s0:s0 + chunk]
        for side, store in ((1, PR), (-1, PL)):
            pos = (uu[:, None] + side * t[None, :]) % Q
            pt = np.where(isr[pos], rho0 * St[None, :], 0.0)
            surv = np.cumprod(1 - pt, axis=1)
            surv = np.concatenate([np.ones((len(uu), 1)), surv[:, :-1]], axis=1)
            w = pt * surv                                # P(neighbour at offset t)
            for q in qs:
                cls = (uu[:, None] + side * t[None, :]) % q
                for b in range(q):
                    store[q][s0:s0 + len(uu), b] = (w * (cls == b)).sum(1)
    return Q, U, PR, PL, rho0
