# Independent numerical test: P = 3, k = 2, levels N = 3q, q prime in [2000, 4000] (X = 6000, beta = 2).
import sys, math, time; sys.path.insert(0, '/tmp/claude-1000/referee_martin/scripts')
from ref_core import *
t0 = time.time()
X = 6000; beta = 2
qs = [int(p) for p in pari.primes([2000, 4000])]
# eigenspace sizes from trace forms: tr W_Q = (-1)^omega(Q) a_Q(new trace form), k = 2
sizes = {}
for q in qs:
    N = 3*q
    tf = pari.mftraceform([N, 2], 0)
    d = int(pari.mfdim([N, 2], 0))
    w3 = -int(pari.mfcoef(tf, 3)); wq = -int(pari.mfcoef(tf, q)); wN = int(pari.mfcoef(tf, N))
    for e3 in (1, -1):
        for eq in (1, -1):
            num = d + e3*w3 + eq*wq + e3*eq*wN
            assert num % 4 == 0
            sizes[(q, e3, eq)] = num//4
print('sizes done', round(time.time()-t0, 1), 's;  sum|H^eps|/(sum dim/4):',
      [round(sum(sizes[(q,a,b)] for q in qs)/(sum(int(pari.mfdim([3*q,2],0)) for q in qs)/4), 3) for a in (1,-1) for b in (1,-1)])

def H_f(n): return float(H(n))
def traces(l, q):
    """t_Q^(3q)(l) for Q in {1, 3, q, 3q}; the fixed ones via the q-independent class numbers (own loop)."""
    N = 3*q
    out = {}
    out[q] = float(tQ(N, q, l)); out[N] = float(tQ(N, N, l))
    out[1] = float(tQ(N, 1, l)); out[3] = float(tQ(N, 3, l))
    return out

def predD(y):
    return D_dens(y, 1), D_dens(y, 3)   # D_q-type (Q0 = 1, M = 3), D_N-type (Q0 = 3, M = 1)

targets = [0.2, 0.5, 1.0, 1.5]
rows = []
for xt in targets:
    l0 = int(pari.nextprime(int(xt*X)))
    ls = [l0]
    while len(ls) < 6: ls.append(int(pari.nextprime(ls[-1]+1)))
    for l in ls:
        tr = {q: traces(l, q) for q in qs if q != l}
        for e3 in (1, -1):
            for eq in (1, -1):
                den = sum(sizes[(q, e3, eq)] for q in qs)
                numq = {q: (eq*tr[q][q] + e3*eq*tr[q][3*q]) / 4 for q in tr}   # q-terms
                numf = {q: (tr[q][1] + e3*tr[q][3]) / 4 for q in tr}          # fixed terms
                A = (sum(numq.values()) + sum(numf.values())) / den
                Afix = sum(numf.values()) / den
                # standard error by ratio estimator over levels
                r = [numq[q] + numf[q] - A*sizes[(q, e3, eq)] for q in tr]
                se = math.sqrt(sum(x*x for x in r)) / den
                Mlim = M_eps(l/X, beta, e3, eq, n=400)
                # finite-X model: prime-weighted average of the limit densities
                wsum = 0.0; wden = 0.0
                for q in qs:
                    if q == l: continue
                    a, b = predD(l/(3*q)); wsum += (q-1)*eq*(a + e3*b); wden += (q-1)
                Mq = wsum/wden
                # alternative: W_N-only (drop the D_q-type term) and the "other sign at 3"
                rows.append((xt, l, e3, eq, A, Afix, se, Mlim, Mq))
        print('l=%d done %.0fs' % (l, time.time()-t0), flush=True)
import json
json.dump(rows, open('/tmp/claude-1000/referee_martin/rows.json', 'w'))
print('%5s %6s %3s %3s %8s %8s %8s %8s %8s' % ('x', 'l', 'e3', 'eq', 'A', 'A_fixed', 'se', 'M_lim', 'M_qsum'))
for r in rows:
    print('%5.2f %6d %3d %3d %8.3f %8.3f %8.3f %8.3f %8.3f' % r)
