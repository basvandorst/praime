"""Round 28 (diagnostic): the expansion's bias for longer gaps. Exact Hardy-Littlewood share (round-22 enumeration)
against order 4 and order 4 + overlap correction, at lam = 17.8 and 23, for single gaps h = 74..90."""
import os, sys, json, math, time
import theory2 as T2
import numpy as np


def U_sums(h, batch=200000):
    """as r22.U_sums, but streaming: sets are flushed to S_many in batches, so memory stays bounded."""
    PS = [int(p) for p in T2.small_primes(h + 1) if p > 2]
    B = [0, h]
    occ = {p: [0] * p for p in PS}; cov = {p: 0 for p in PS}
    for x in B:
        for p in PS:
            if occ[p][x % p] == 0:
                cov[p] += 1
            occ[p][x % p] += 1
    pts = list(range(2, h, 2))
    U = np.zeros(h); buf = {}; cnt = [0]
    cur = []

    def flush(k):
        if buf.get(k):
            arr = np.array([B + list(a) for a in buf[k]], dtype=np.int64)
            U[k] += float(T2.S_many(arr, P=max(h, 3)).sum()); buf[k] = []

    def rec(start):
        k = len(cur); cnt[0] += 1
        buf.setdefault(k, []).append(tuple(cur))
        if len(buf[k]) >= batch:
            flush(k)
        for i in range(start, len(pts)):
            x = pts[i]
            if any(occ[p][x % p] == 0 and cov[p] == p - 1 for p in PS):
                continue
            for p in PS:
                if occ[p][x % p] == 0:
                    cov[p] += 1
                occ[p][x % p] += 1
            cur.append(x); rec(i + 1); cur.pop()
            for p in PS:
                occ[p][x % p] -= 1
                if occ[p][x % p] == 0:
                    cov[p] -= 1
    rec(0)
    for k in list(buf):
        flush(k)
    return np.trim_zeros(U, 'b'), cnt[0]

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
path = os.path.join(OUT, 'r28.json')
res = json.load(open(path)) if os.path.exists(path) else {}
for h in [int(x) for x in sys.argv[1:]]:
    if str(h) in res:
        continue
    t0 = time.time()
    U, n = U_sums(h)
    pt = T2.Pattern([0, h], 4)
    r = {'sets': n}
    for lam in (17.8, 23.0):
        ex = sum((-1) ** j * U[j] / lam ** (j + 1) for j in range(len(U)))
        L4 = pt.log_avoid(lam, 4); oc = T2.overlap_correction(pt, lam)
        sh = lambda la: pt.SB / lam * math.exp(la)
        r['lam%.1f' % lam] = {'exact': ex, 'order4_rel': sh(L4) / ex - 1, 'order4_overlap_rel': sh(L4 + oc) / ex - 1}
    r['secs'] = round(time.time() - t0)
    res[str(h)] = r
    json.dump(res, open(path, 'w'), indent=1)
    print(h, r, flush=True)
