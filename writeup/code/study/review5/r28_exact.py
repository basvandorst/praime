"""Own exact Hardy-Littlewood share of gap h (sum over admissible interior sets with signs), vectorised BFS.
share(h) = sum_{A subset interior, {0,h} u A admissible} (-1)^|A| S({0,h} u A) / lam^(|A|+1)."""
import sys, math, json
import numpy as np
from lib import primes_upto

PR = primes_upto(20_000_000)


def tailG(k, P):
    """prod_{p > P} (1-k/p)/(1-1/p)^k with primes to 2e7 and an analytic tail."""
    p = PR[PR > P].astype(np.float64)
    s = float(np.sum(np.log1p(-k / p) - k * np.log1p(-1 / p)))
    return math.exp(s - k * (k - 1) / 2 / (2e7 * math.log(2e7)))


def exact(h, lams):
    pts = list(range(2, h, 2))
    ODD = [int(p) for p in primes_upto(h) if p > 2]
    TR = [p for p in ODD if 2 * p <= h]        # primes where even points can collide; others give nu = k
    P1 = ODD[-1]
    sumlog1 = {p: math.log1p(-1 / p) for p in ODD}
    G = {k: tailG(k, P1) for k in range(2, h // 2 + 3)}
    bits = np.array([[1 << (x % p) for p in TR] for x in pts], dtype=np.uint64)
    occ = np.array([[(1 << 0) | (1 << (h % p)) for p in TR]], dtype=np.uint64)
    last = np.array([-1]); k = 2
    U = {}
    nsets = 0
    Parr = np.array(TR, dtype=np.int64)

    def Ssum(occ, k):
        nu = np.bitwise_count(occ).astype(np.float64)             # (n, len(TR))
        logS = np.sum(np.log1p(-nu / Parr), axis=1)
        logS += sum(math.log1p(-k / p) for p in ODD if 2 * p > h)
        logS -= k * sum(sumlog1.values())
        return float(np.sum(np.exp(logS))) * 2 ** (k - 1) * G[k]
    U[0] = Ssum(occ, 2); nsets = 1
    while len(last):
        new_occ, new_last = [], []
        for i in range(len(pts)):
            sel = last < i
            if not sel.any():
                continue
            o = occ[sel] | bits[i]
            ok = np.all(np.bitwise_count(o) < Parr, axis=1)
            if ok.any():
                new_occ.append(o[ok]); new_last.append(np.full(int(ok.sum()), i))
        if not new_occ:
            break
        occ = np.concatenate(new_occ); last = np.concatenate(new_last); k += 1
        nsets += len(last)
        U[k - 2] = Ssum(occ, k)
    res = {lam: sum((-1) ** j * U[j] / lam ** (j + 1) for j in U) for lam in lams}
    return res, nsets


out = {}
for h in [int(x) for x in sys.argv[1:]]:
    res, n = exact(h, (17.8, 23.0))
    out[h] = {'sets': n, 'exact_17.8': res[17.8], 'exact_23.0': res[23.0]}
    print(h, out[h], flush=True)
json.dump(out, open('/home/bas/scratch/study-review5/r28_exact_%s.json' % '_'.join(sys.argv[1:]), 'w'), indent=1)
