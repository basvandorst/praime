"""f(n) = G(h(n)): G a prime-rich quadratic, h(n) = a*n(n-1)/2 + b*n + c integer-valued quadratic (or cubic with d*C(n,3)).
Longest run of consecutive n in [-NW, NW] with |f| prime and all values distinct. usage: comp.py DEG"""
import sys, numpy as np, gmpy2, itertools, time
DEG = int(sys.argv[1]) if len(sys.argv) > 1 else 2
BASES = {'euler': (1, 1, 41), 'ruby': (36, 18, -1801), 'fr47': (47, 9, -5209), 'sp103': (103, 31, -3391),
         'legendre2': (2, 0, 29), 'h4': (4, 4, 59)}
NW = 70; M = 400000
n = np.arange(-NW, NW + 1, dtype=np.int64)
t0 = time.time(); best_all = []
for name, (A, B, C) in BASES.items():
    m = np.arange(-M, M + 1, dtype=np.int64)
    G = A * m * m + B * m + C
    isp = np.fromiter((gmpy2.is_prime(abs(int(v))) for v in G), bool, len(G))
    best = (0,)
    arange = range(1, 13) if DEG == 2 else range(1, 7)
    for a in list(arange) + [-x for x in arange]:
        for d in ([0] if DEG == 2 else [1, 2, 3, 6, -1, -2, -3, -6]):
            for b in range(0, abs(a) if DEG == 2 else 40):
                h0 = a * n * (n - 1) // 2 + b * n + (d * n * (n - 1) * (n - 2) // 6)
                for c in range(-3000, 3001):
                    h = h0 + c
                    ok = (np.abs(h) <= M)
                    if not ok.all():
                        ok &= True
                    hv = np.clip(h, -M, M) + M
                    pr = isp[hv] & (np.abs(h) <= M)
                    # longest run
                    if pr.sum() < best[0]: continue
                    x = np.concatenate([[0], pr.astype(np.int8), [0]]); dd = np.diff(x)
                    st = np.nonzero(dd == 1)[0]; en = np.nonzero(dd == -1)[0]; L = en - st
                    for i in np.argsort(-L)[:3]:
                        if L[i] <= best[0]: break
                        vals = G[hv[st[i]:en[i]]]
                        # longest sub-run with distinct |values|
                        seen = {}; lo = 0; bl = 0; bs = 0
                        for j, v in enumerate(np.abs(vals).tolist()):
                            if v in seen and seen[v] >= lo: lo = seen[v] + 1
                            seen[v] = j
                            if j - lo + 1 > bl: bl = j - lo + 1; bs = lo
                        if bl > best[0]:
                            best = (bl, a, b, c, d, int(n[st[i] + bs]))
    print(name, 'best distinct run', best, round(time.time() - t0), 's', flush=True)
