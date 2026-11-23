import numpy as np, cypari2
from common import *
pari = cypari2.Pari()
d = complex_zeros()
print('n chars', len(d))
ls = np.arange(0.5, 8.01, 0.5)
rows = []
bad = []
for (N, m), g in sorted(d.items()):
    G = pari.znstar(N, 1)
    cond = int(pari.znconreyconductor(G, m)) if True else None
    # znconreyconductor returns an integer if primitive? It returns conductor (int) or [cond, char]
    c = pari.znconreyconductor(G, m)
    try:
        cond = int(c)
    except Exception:
        cond = int(c[0])
    order = int(pari.charorder(G, m))
    # smallest Conrey label of order > 2 primitive?
    smallest = None
    for mm in range(1, N):
        if pari.gcd(mm, N) != 1: continue
        cc = pari.znconreyconductor(G, mm)
        try: ccond = int(cc)
        except Exception: ccond = int(cc[0])
        if ccond == N and int(pari.charorder(G, mm)) > 2:
            smallest = mm; break
    if cond != N or order <= 2 or smallest != m:
        bad.append((N, m, cond, order, smallest))
    gp, gn = g[g > 0], -g[g < 0][::-1]
    gn = np.sort(gn)
    sp, Lp, np_ = unfold_stats(gp, N, 100, 500, ls)
    sn, Ln, nn = unfold_stats(gn, N, 100, 500, ls)
    spt = unfold_stats(gp, N, 100, 500, ls, frac=0.1)[0]; snt = unfold_stats(gn, N, 100, 500, ls, frac=0.1)[0]
    rng = np.random.default_rng(N)
    spr = unfold_stats(gp, N, 100, 500, ls, mode='rand', rng=rng, nrand=2000)[0]; snr = unfold_stats(gn, N, 100, 500, ls, mode='rand', rng=rng, nrand=2000)[0]
    # pooled: weight by zero count (pool window counts) ~ approx mean of the two
    rows.append(dict(N=N, m=m, order=order, S=(sp + sn) / 2, Sp=sp, Sn=sn, St=(spt + snt) / 2, Sr=(spr + snr) / 2,
                     L=(Lp * np_ + Ln * nn) / (np_ + nn), nz=np_ + nn, ps=primes_of(N),
                     ming=g.min(), maxabs=np.abs(g).max(), minabs=np.abs(g).min()))
print('bad', bad)
Ns = [r['N'] for r in rows]
from sympy import isprime, primerange
pr = [r['N'] for r in rows if isprime(r['N'])]
print('prime N', len(pr), pr)
for k in (3, 4, 5, 7):
    print(k, 'N', sorted(r['N'] for r in rows if r['N'] % k == 0 and not isprime(r['N'])))
print('expected primes every 2nd 101..397', list(primerange(101, 398))[::2])
print('orders', sorted(set(r['order'] for r in rows)))
print('minabs gamma', min(r['minabs'] for r in rows))
il1 = 1
x = np.array([sum(dS(1, p, r['L']) for p in r['ps']) for r in rows])
for key in ['S', 'St', 'Sr', 'Sp', 'Sn']:
    y = np.array([r[key][il1] for r in rows])
    b, se = ols(np.column_stack([np.ones(len(y)), x]), y)
    print(key, 'beta %.3f +- %.3f, a %.4f' % (b[1], se[1], b[0]))
y = np.array([r['S'][il1] for r in rows])
rng = np.random.default_rng(3); bb = []
for _ in range(2000):
    i = rng.integers(0, len(y), len(y)); bb.append(ols(np.column_stack([np.ones(len(y)), x[i]]), y[i])[0][1])
print('bootstrap beta sd %.3f' % np.std(bb))
# group means
isp = np.array([isprime(r['N']) for r in rows])
print('prime N mean %.4f n=%d' % (y[isp].mean(), isp.sum()))
for k in (2, 3, 5, 7):
    s = np.array([r['N'] % k == 0 for r in rows]); print(k, '%.4f n=%d' % (y[s].mean(), s.sum()))
# per-prime fit with indicators, multi-l, to test integer/half-integer alternation
Lr = np.array([r['L'] for r in rows])
I = np.array([[1.0 if r['N'] % p == 0 else 0.0 for p in (2, 3, 5, 7)] for r in rows])
rest = lambda l: np.array([sum(dS(l, p, r['L']) for p in r['ps'] if p > 7) for r in rows])
X = np.column_stack([np.ones(len(rows)), I])
print('complex per-l, p=2..7 fitted/pred (offset for primes>7):')
for j, l in enumerate(ls):
    Sy = np.array([r['S'][j] for r in rows]) - rest(l)
    b, se = ols(X, Sy)
    P = [np.mean([dS(l, p, Lr[k]) for k in range(len(rows)) if I[k, i]]) for i, p in enumerate((2, 3, 5, 7))]
    print(l, np.round(b[1:], 4), np.round(P, 4), 'slope %.2f' % (np.dot(b[1:], P) / np.dot(P, P)))
np.savez('/home/bas/scratch/review-logic3/complex.npz', S=np.array([r['S'] for r in rows]), L=Lr, N=np.array(Ns))
