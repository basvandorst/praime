# Item F: Goldbach strands against the prime race, mod 12 and mod 8, N = 4,000,000.
import numpy as np, math
from sympy import primerange, isprime

N = 4_000_000
s = np.ones(N + 1, bool); s[:2] = False
for i in range(2, int(N**0.5) + 1):
    if s[i]: s[i*i::i] = False
ar = np.arange(N + 1)
FL = 1 << 23

SPEC = {
    12: {0: [(1, 11, 1), (11, 1, 1), (5, 7, -1), (7, 5, -1)],
         2: [(1, 1, 1), (7, 7, -1)],
         6: [(1, 5, 1), (5, 1, 1), (7, 11, -1), (11, 7, -1)],
         10: [(5, 5, 1), (11, 11, -1)]},
    8: {0: [(1, 7, 1), (7, 1, 1), (3, 5, -1), (5, 3, -1)],
        2: [(1, 1, 1), (5, 5, -1)],
        4: [(1, 3, 1), (3, 1, 1), (5, 7, -1), (7, 5, -1)],
        6: [(1, 5, 1), (5, 1, 1), (3, 3, -1), (7, 7, -1)]},
}
S0 = {12: 5, 8: 3}
CLAIM = {
    12: {'corr': [0.993, 0.993, 0.989, 0.990], 'mad': [3.99, 2.27, 2.34, 2.29], 'madS': [19.41, 11.07, 10.49, None]},
    8: {'corr': [0.986, 0.987, 0.988, 0.989], 'mad': [1.82, 1.69, 2.33, 1.60], 'madS': [8.92, 8.35, 12.23, 7.02]},
}

def SQ(s0):
    """S(n), Q(n) for all n (even and odd), primes s >= s0."""
    S = np.zeros(N + 1, np.int32); Q = np.zeros(N + 1, np.int32)
    for p in primerange(s0, math.isqrt(N) + 1):
        sq = p * p
        S[sq + 1:] += s[1:N + 1 - sq]          # n - s^2 >= 1, s^2 < n
        Q[sq:] += 1
    return S, Q

edges = np.exp(np.linspace(math.log(N / 16), math.log(N), 41))
print('last edge = %r (n = N included: %s)' % (edges[-1], N < edges[-1]))

def direct_delta(n, spec, M):
    c = {}
    for p in primerange(2, n):
        if s[n - p]:
            c[(p % M, (n - p) % M)] = c.get((p % M, (n - p) % M), 0) + 1
    return sum(sg * c.get((a, b), 0) for a, b, sg in spec)

out = []
for M in (12, 8):
    classes = [a for a in range(M) if math.gcd(a, M) == 1]
    ind = {a: (s & (ar % M == a)).astype(np.float64) for a in classes}
    FT = {a: np.fft.rfft(ind[a], FL) for a in classes}
    pic = {a: np.cumsum(ind[a]).astype(np.int64) for a in classes}
    del ind
    S, Q = SQ(S0[M])
    for ci, (r, spec) in enumerate(SPEC[M].items()):
        n = np.arange(r if r > 0 else M, N + 1, M)
        n = n[n >= 4]
        D = np.zeros(len(n), np.int64)
        SIG = np.zeros(len(n), np.int64)
        for a, b, sg in spec:
            c = np.fft.irfft(FT[a] * FT[b], FL)[:N + 1]
            assert np.abs(c - np.rint(c)).max() < 0.01
            D += sg * np.rint(c).astype(np.int64)[n]
            SIG += sg * (pic[a][n] + pic[b][n])
        for nn in (n[len(n) // 3], n[-1]):
            assert direct_delta(int(nn), spec, M) == D[np.searchsorted(n, nn)]
        Sn = S[n - 2] if (M == 12 and r == 10) else S[n]
        Qn = Q[n]
        okQ = Qn > 0
        pred = np.where(okQ, Sn / np.maximum(Qn, 1) * SIG, 0.0)
        mD, mP, mS = [], [], []
        for k in range(40):
            sel = (n >= edges[k]) & (n < edges[k + 1]) & okQ
            mD.append(D[sel].mean()); mP.append(pred[sel].mean()); mS.append((-S[n][sel]).mean())
        mD, mP, mS = map(np.array, (mD, mP, mS))
        corr = np.corrcoef(mD, mP)[0, 1]
        mad = np.mean(np.abs(mD - mP))
        madS = np.mean(np.abs(mD - mS))
        cl = CLAIM[M]
        line = ('mod %2d n=%2d: corr %.4f (claim %.3f)  MAD pred %.3f (claim %.2f)  MAD -S %.3f (claim %s)' % (
            M, r, corr, cl['corr'][ci], mad, cl['mad'][ci], madS, cl['madS'][ci]))
        if M == 12 and r == 10:
            # variant: Q(n-2) together with S(n-2)
            pred2 = np.where(Q[n - 2] > 0, S[n - 2] / np.maximum(Q[n - 2], 1) * SIG, 0.0)
            mP2 = np.array([pred2[(n >= edges[k]) & (n < edges[k + 1]) & okQ].mean() for k in range(40)])
            line += '  [variant Q(n-2): corr %.4f MAD %.3f]' % (np.corrcoef(mD, mP2)[0, 1], np.mean(np.abs(mD - mP2)))
        print(line, flush=True)
        out.append(line)
    del FT
