"""Seeded choice of the blind test set (no zero of any L-function is computed here).
Rule: for each class, enumerate all candidates in a fixed box, sort by a canonical key, then rng = numpy default_rng(SEED + class
index) permutes them; the first n_test are the test set, the next n_proxy are proxies (similar L-functions used only for
the PARI checks and timings, never tested). SEED = 20261125.
 T1  primitive Dirichlet characters chi mod q, 8 <= q <= 60, order >= 3, as conjugate pairs (Conrey m < m^-1 mod q): 20 pairs (40 L).
 T2  newforms of weight k in {4,6,8,10,12}, squarefree level 2 <= N <= 60, trivial character, rational coefficients: 30.
 T3  Sym^2 E, E semistable with 20 <= N <= 150, one curve per isogeny class (first a-invariants in the box order); classes of
     11a, 14a, 15a, 17a, 19a, 37a and the class of [1,0,1,-19,26] (probed) excluded: 30.
 T4a Sym^3 E, same curve box, 20 <= N <= 80: 15.     T4b genus-2 curves y^2 + h y = f, f in {-1,0,1}^7 (deg <= 6),
     h in the 16 0/1-polynomials of degree <= 3, conductor < 1000, not 169, 249, 277, one curve per L-function, mean
     |lambda_p|^2 over p <= 10^4 in [0.8, 1.2] (primitive): 15.
 T5  complex (hence non-normal) cubic fields with -500 < d < 0, d not -23, -31, from x^3 + a x^2 + b x + c, a in {-1,0,1},
     |b|, |c| <= 12, polredabs: 30.   Pairs: 10 disjoint pairs of the 138 round-144 curves (fresh_zeros.npz keys).
Curve box: a1, a3 in {0,1}, a2 in {-1,0,1}, |a4| <= 20, |a6| <= 40; isogeny class key = (N, a_p for p <= 97)."""
import sys, os, json, math, time
import numpy as np, cypari2
pari = cypari2.Pari(); pari.allocatemem(800 * 10 ** 6)
SEED = 20261125
out = dict(seed=SEED, rule=__doc__)
def pick(cands, cls_idx, n_test, n_proxy):
    rng = np.random.default_rng(SEED + cls_idx)
    perm = rng.permutation(len(cands))
    sel = [cands[i] for i in perm]
    return sel[:n_test], sel[n_test:n_test + n_proxy], len(cands)
t0 = time.time()
# T1
pairs = []
for q in range(8, 61):
    if q % 4 == 2: continue                       # no primitive characters
    G = pari('znstar(%d,1)' % q)
    for m in range(1, q):
        if math.gcd(m, q) != 1: continue
        mi = pow(m, -1, q)
        if mi <= m: continue
        chi = pari.znconreychar(G, m)
        if int(pari.zncharconductor(G, chi)) != q: continue
        if int(pari.charorder(G, chi)) < 3: continue
        pairs.append([q, m, mi])
pairs.sort()
t1, p1, n1 = pick(pairs, 1, 20, 3)
out['T1'] = dict(n_candidates=n1, test=t1, proxy=p1)
print('T1 candidates (pairs)', n1, 'secs', round(time.time() - t0, 1), flush=True)
# T2
forms = []
for k in (4, 6, 8, 10, 12):
    for N in range(2, 61):
        if not pari.issquarefree(N): continue
        mf = pari.mfinit([N, k], 0)
        if int(pari.mfdim(mf)) == 0: continue
        fl = pari.mffields(mf)
        for i in range(len(fl)):
            if int(pari.poldegree(fl[i])) == 1:
                forms.append([k, N, i])
forms.sort()
t2, p2, n2 = pick(forms, 2, 30, 3)
out['T2'] = dict(n_candidates=n2, test=t2, proxy=p2)
print('T2 candidates', n2, 'secs', round(time.time() - t0, 1), flush=True)
# curves
PR = [int(p) for p in pari.primes([2, 97])]
classes = {}
for a1 in (0, 1):
    for a2 in (-1, 0, 1):
        for a3 in (0, 1):
            for a4 in range(-20, 21):
                for a6 in range(-40, 41):
                    a = [a1, a2, a3, a4, a6]
                    try:
                        E = pari.ellinit(a)
                        if pari.ellglobalred(E) is None: continue
                    except Exception:
                        continue
                    if len(E) == 0: continue
                    N = int(pari.ellglobalred(E)[0])
                    if N > 150 or N < 11 or not pari.issquarefree(N): continue
                    key = (N,) + tuple(int(pari.ellap(E, p)) for p in PR)
                    if key not in classes:
                        classes[key] = a
print('curve classes (semistable, N <= 150) in the box', len(classes), 'secs', round(time.time() - t0, 1), flush=True)
excl_ainv = [[0, -1, 1, -10, -20], [1, 0, 1, 4, -6], [1, 1, 1, -10, -10], [1, -1, 1, -1, -14], [0, 1, 1, -9, -15], [0, 0, 1, -1, 0], [1, 0, 1, -19, 26]]
excl = set()
for a in excl_ainv:
    E = pari.ellinit(a); N = int(pari.ellglobalred(E)[0]); excl.add((N,) + tuple(int(pari.ellap(E, p)) for p in PR))
cl = sorted([[k[0], v] for k, v in classes.items() if k not in excl])
c3 = [c for c in cl if 20 <= c[0] <= 150]
t3, p3, n3 = pick(c3, 3, 30, 3)
c4 = [c for c in cl if 20 <= c[0] <= 80]
t4, p4, n4 = pick(c4, 4, 15, 2)
out['T3'] = dict(n_candidates=n3, test=t3, proxy=p3, excluded_classes=excl_ainv)
out['T4a'] = dict(n_candidates=n4, test=t4, proxy=p4)
print('T3 candidates', n3, 'T4a candidates', n4, flush=True)
# genus 2
import itertools
hs = []
for bits in range(16):
    hs.append('+'.join(['x^%d' % j for j in range(4) if bits >> j & 1]) or '0')
g2 = {}
cnt = 0
for coefs in itertools.product((-1, 0, 1), repeat=7):
    if coefs[6] == 0 and coefs[5] == 0: continue
    f = '+'.join('(%d)*x^%d' % (c, j) for j, c in enumerate(coefs) if c != 0)
    for h in hs:
        try:
            F = pari('4*(%s)+(%s)^2' % (f, h))
            if int(pari.poldegree(F)) < 5 or int(pari.poldegree(pari.gcd(F, pari.deriv(F)))) > 0: continue
            r = pari('genus2red([%s,%s])' % (f, h))
            N = int(r[0])
        except Exception:
            continue
        cnt += 1
        if N >= 1000 or N in (169, 249, 277): continue
        g2.setdefault(N, []).append([f, h])
print('genus-2 curves tried', cnt, 'conductors < 1000 found', len(g2), 'secs', round(time.time() - t0, 1), flush=True)
# one curve per L-function: key = (N, first 60 coefficients)
g2c = {}
for N in sorted(g2):
    for f, h in g2[N]:
        L = pari('lfungenus2([%s,%s])' % (f, h))
        an = pari.lfunan(L, 60)
        key = (N,) + tuple(int(x) for x in an)
        if key not in g2c: g2c[key] = [N, f, h]
print('genus-2 L-functions with N < 1000', len(g2c), 'secs', round(time.time() - t0, 1), flush=True)
g2list = []
for key, (N, f, h) in sorted(g2c.items()):
    L = pari('lfungenus2([%s,%s])' % (f, h))
    an = pari.lfunan(L, 10000)
    ps = [int(p) for p in pari.primes([2, 10000]) if N % int(p) != 0]
    m2 = float(np.mean([(float(an[p - 1]) / math.sqrt(p)) ** 2 for p in ps]))
    g2list.append([N, f, h, round(m2, 4)])
prim = [x for x in g2list if 0.8 <= x[3] <= 1.2]
t4b, p4b, n4b = pick(prim, 5, 15, 2)
out['T4b'] = dict(n_candidates=n4b, test=t4b, proxy=p4b, all_found=g2list)
print('genus-2 primitive candidates', n4b, '(non-primitive', len(g2list) - n4b, ')  secs', round(time.time() - t0, 1), flush=True)
# T5 cubic fields
fields = {}
for a in (-1, 0, 1):
    for b in range(-12, 13):
        for c in range(-12, 13):
            if c == 0: continue
            P = pari('x^3+(%d)*x^2+(%d)*x+(%d)' % (a, b, c))
            if not pari.polisirreducible(P): continue
            d = int(pari.nfdisc(P))
            if not (-500 < d < 0) or d in (-23, -31): continue
            red = str(pari.polredabs(P))
            fields[red] = d
fl5 = sorted([[d, p] for p, d in fields.items()])
t5, p5, n5 = pick(fl5, 6, 30, 2)
out['T5'] = dict(n_candidates=n5, test=t5, proxy=p5, all_found=fl5)
print('cubic fields', n5, 'secs', round(time.time() - t0, 1), flush=True)
# pairs
z = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', '..', 'data', 'study', 'r144', 'fresh_zeros.npz'))
keys = sorted(z.keys())
rng = np.random.default_rng(SEED + 7)
perm = [keys[i] for i in rng.permutation(len(keys))]
pairs_t = [[perm[2 * i], perm[2 * i + 1]] for i in range(10)]
pairs_p = [[perm[20 + 2 * i], perm[21 + 2 * i]] for i in range(2)]
out['T5pairs'] = dict(n_candidates=len(keys), test=pairs_t, proxy=pairs_p)
json.dump(out, open('testset.json', 'w'), indent=1)
print('done', round(time.time() - t0, 1))
