"""Check of the restricted Montgomery-Soundararajan sums against the closed forms of derivation.md section 3.
T_q(c; H) = sum_{0<|h|<H, h = c mod q} (H - |h|) (S'_q(h) - 1), S'_q = singular series without the primes p | q.
Closed forms (derivation.md, eq. (3.6), (3.7)):
  c = 0 mod q : T = -(phi/q) H [ln(H/q) + S_q + gamma + ln 2pi - 1] + H
  c != 0      : T = -H (phi(q)/phi(q/g)) Lambda(q/g)/q,   g = gcd(c, q)
and the class average of the full S over h = c mod q: A_q(c) = prod_{p|q} (p/(p-1) if p|c else p(p-2)/(p-1)^2)."""
import math, sys
import numpy as np
from sympy import factorint, totient

HMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 2_000_000
N = 2 * HMAX + 10


def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, math.isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


PR = primes_upto(N)
# log of prod_{p | h, p > 2} (p-1)/(p-2), h = 0..N-1
logF = np.zeros(N)
for p in PR[1:]:
    p = int(p)
    logF[p::p] += math.log((p - 1) / (p - 2))
C2_trunc = float(np.exp(np.sum(np.log1p(-1.0 / (PR[1:].astype(float) - 1) ** 2))))   # twin constant from primes < N
C2 = 0.66016181584686957393                       # exact twin prime constant (the truncated product is 3e-8 too high,
# which adds 3e-8 H^2/q to T: visible as a uniform shift 0.03/q at H = 1e6)
h = np.arange(N)
S = np.where(h % 2 == 0, 2 * C2 * np.exp(logF), 0.0)                      # full singular series, h >= 1


def local(p, hh):
    return np.where(hh % p == 0, p / (p - 1), p * (p - 2) / (p - 1) ** 2)


def Sprime(q):
    """S'_q(h) = S(h) / prod_{p|q} local factor at p (for p = 2 and odd h: S' = C_q-type value with the 2 removed)."""
    out = 2 * C2 * np.exp(logF)                     # value of prod over p>2 times 2C2 (local factor at 2 for even h)
    out = out / 2.0                                 # remove the local factor at 2 (=2 for even h) -> product over p>2 only
    # now out = C2 * prod_{p|h,p>2}(p-1)/(p-2) = prod_{p>2} local(p,h)
    for p in factorint(q):
        if p == 2:
            continue
        out = out / local(p, h)
    if q % 2 == 1:                                  # keep the factor at 2
        out = out * np.where(h % 2 == 0, 2.0, 0.0)
    return out


def T(q, c, H, Sp):
    hp = np.arange(1, H)
    w = (H - hp).astype(float)
    tot = 0.0
    for sgn in (1, -1):
        sel = ((sgn * hp - c) % q) == 0
        tot += float(np.sum(w[sel] * (Sp[hp[sel]] - 1.0)))
    return tot


def T0_closed(q, H):
    phi = int(totient(q)); Sq = sum(math.log(p) / (p - 1) for p in factorint(q))
    return -(phi / q) * H * (math.log(H / q) + Sq + 0.5772156649015329 + math.log(2 * math.pi) - 1) + H


def Tc_closed(q, c, H):
    g = math.gcd(c, q); m = q // g
    f = factorint(m)
    Lam = math.log(list(f)[0]) if len(f) == 1 else 0.0
    return -H * (int(totient(q)) / int(totient(m))) * Lam / q


if __name__ == '__main__':
    print('C2 exact %.12f, from primes < %d: %.12f' % (C2, N, C2_trunc))
    # class averages of the full S
    for q in (3, 4, 10, 30):
        for c in range(q):
            hh = np.arange(1, 2 * HMAX)
            sel = hh % q == c
            A = 1.0
            for p in factorint(q):
                A *= (p / (p - 1)) if c % p == 0 else (p * (p - 2) / (p - 1) ** 2)
            if c < 4 or c == q - 1:
                print('avg S over h=%d mod %d: %.5f  closed %.5f' % (c, q, S[hh[sel]].mean(), A))
    for q in (2, 3, 4, 5, 6, 7, 9, 10, 12, 30, 210):
        Sp = Sprime(q)
        for H in (HMAX // 10, HMAX):
            row = []
            for c in range(q):
                if q > 12 and c not in (0, 2, 4, 6, 10, 12, 14, 30, 42, 70, 84, 105, 1):
                    continue
                t = T(q, c, H, Sp)
                cl = T0_closed(q, H) if c == 0 else Tc_closed(q, c, H)
                row.append('c=%d: %.4f/%.4f' % (c, t / H, cl / H))
            print('q=%d H=%d (T/H data/closed):' % (q, H), '  '.join(row), flush=True)
