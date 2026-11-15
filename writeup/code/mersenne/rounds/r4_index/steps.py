"""Step rules on the derived sequences: the grammar of section 6 (grammar.py, copied unchanged), reduced to
  d1     : b_{n+1} = f(b_n), every primitive (1,187)
  inv-d1 : b_n = f(b_{n+1}) (1,187)
  d2c    : g(f(x)) with f, g both in the 'core2' set (237^2 = 56,169)
  tt     : the 16,884 two-term rules b_{n+2} = h(b_{n+1}, b_n) and their adjustments (search.py, stage tt)
All rows (9 derived sequences x (1 real + 300 null)) are scored at once. Padding is -2 (never produced by a map).
Per row: exact (max #transitions hit), run0 (consecutive hits from the first transition), runany (longest run).
OOS: transitions whose terms are all in-sample (exponent index <= 25) choose the rule (max in-sample hits);
over the tied rules, mean and max hits on the other transitions.
usage: steps.py STAGE  -> res_STAGE.pkl
"""
import sys, time, pickle, math, numpy as np
from grammar import build, apply, clean, bitlen, popcount, revbin, digit_op, ndig

STAGE = sys.argv[1]
D = np.load('data.npz')
NAMES = [str(x) for x in D['names']]
Lmax = max(D[n].shape[1] for n in NAMES)
blocks, Ms, Is = [], [], []
r0 = 0
for nm in NAMES:
    M = D[nm]; I = D[nm + '_ins']
    pad = Lmax - M.shape[1]
    Ms.append(np.pad(M, ((0, 0), (0, pad)), constant_values=-2))
    Is.append(np.pad(I, ((0, 0), (0, pad)), constant_values=False))
    blocks.append((nm, r0, r0 + M.shape[0])); r0 += M.shape[0]
M = np.vstack(Ms); INS = np.vstack(Is)
S = M.shape[0]
REALROWS = {b[1]: b[0] for b in blocks}
P = build()
names = [p[0] for p in P]
core2 = [i for i, p in enumerate(P) if 'core2' in p[2]]
ADJ = [i for i, p in enumerate(P) if p[0] in ('x-2', 'x-1', 'x+1', 'x+2', '2x+1', '2x-1', 'revbin', 'revdec',
                                                 'next≡1(mod 4)', 'next≡3(mod 4)', 'next≡1(mod 6)', 'next≡5(mod 6)',
                                                 'next binpal', 'next 2^k-1', 'next 2^k+1', 'nextsamepop', 'next decpal')]
if STAGE == 'tt':
    X = M[:, 1:-1]; X2 = M[:, :-2]; Y = M[:, 2:]; TIN = INS[:, 2:]
    N = np.broadcast_to(np.arange(2, M.shape[1]), X.shape).astype(np.int64)
elif STAGE == 'inv-d1':
    X = M[:, 1:]; Y = M[:, :-1]; TIN = INS[:, 1:]
    N = np.broadcast_to(np.arange(1, M.shape[1]), X.shape).astype(np.int64)
else:
    X = M[:, :-1]; Y = M[:, 1:]; TIN = INS[:, 1:]
    N = np.broadcast_to(np.arange(1, M.shape[1]), X.shape).astype(np.int64)
VAL = (Y >= 1) & (X >= 1) & ((X2 >= 1) if STAGE == 'tt' else True)
TOUT = VAL & ~TIN
T = X.shape[1]
KEYS = ('exact', 'run0', 'runany')
best = {k: np.zeros(S, np.int64) for k in KEYS}
oos = dict(best_in=np.full(S, -1, np.int64), n=np.zeros(S, np.int64), sum=np.zeros(S), max=np.zeros(S, np.int64))
realbest = {nm: (0, []) for nm in NAMES}      # best exact on the real row: (count, [(rule, hit cols)])
realoos = {nm: (-1, []) for nm in NAMES}      # best in-sample on the real row: (count, [(rule, oos hits)])
rr = np.array(sorted(REALROWS))
G = 0


def score(F, name_fn):
    global G
    G += 1
    E = (F == Y) & VAL
    ex = E.sum(1)
    run0 = np.argmin(np.concatenate([E, np.zeros((S, 1), bool)], 1), axis=1)
    if ex.max() >= 2:
        c = np.cumsum(E, axis=1)
        reset = np.maximum.accumulate(np.where(E, 0, c), axis=1)
        runany = (c - reset).max(1)
    else:
        runany = ex
    for k, v in (('exact', ex), ('run0', run0), ('runany', runany)):
        np.maximum(best[k], v, out=best[k])
    s_in = (E & TIN).sum(1); s_out = (E & TOUT).sum(1)
    gt = s_in > oos['best_in']; eq = s_in == oos['best_in']
    oos['max'] = np.where(gt, s_out, np.where(eq, np.maximum(oos['max'], s_out), oos['max']))
    oos['sum'] = np.where(gt, s_out, oos['sum'] + eq * s_out)
    oos['n'] = np.where(gt, 1, oos['n'] + eq)
    oos['best_in'] = np.where(gt, s_in, oos['best_in'])
    for r in rr:
        nm = REALROWS[r]
        e = int(ex[r])
        if e >= realbest[nm][0] and e > 0:
            lst = [] if e > realbest[nm][0] else realbest[nm][1]
            if len(lst) < 30:
                lst.append((name_fn(), [int(i) for i in np.nonzero(E[r])[0]]))
            realbest[nm] = (e, lst)
        si = int(s_in[r])
        if si >= realoos[nm][0]:
            lst = [] if si > realoos[nm][0] else realoos[nm][1]
            if len(lst) < 30 or (int(s_out[r]) > 0 and len(lst) < 300):
                lst.append((name_fn(), int(s_out[r])))
            realoos[nm] = (si, lst)


def ev(i, x):
    return apply(P[i][1], x, N)


t0 = time.time()
if STAGE in ('d1', 'inv-d1'):
    for i in range(len(P)):
        score(ev(i, X), lambda i=i: names[i])
elif STAGE == 'd2c':
    for i in core2:
        Z = ev(i, X)
        for j in core2:
            score(ev(j, Z), lambda i=i, j=j: f'{names[j]} o {names[i]}')
elif STAGE == 'tt':
    a, b = X, X2
    af, bf = a.astype(np.float64), b.astype(np.float64)
    H = []
    for c1 in range(1, 5):
        for c2 in range(-4, 5):
            for c3 in range(-10, 11):
                H.append((f'{c1}a{c2:+d}b{c3:+d}', lambda c1=c1, c2=c2, c3=c3: c1 * a + c2 * b + c3))
    for c in range(-5, 6):
        H.append((f'ab{c:+d}', lambda c=c: np.where(a < 2 ** 26, a * b + c, -1)))
    for c in range(-3, 4):
        H.append((f'floor(a^2/b){c:+d}', lambda c=c: np.floor(af * af / bf) + c))
        H.append((f'round(a^2/b){c:+d}', lambda c=c: np.rint(af * af / bf) + c))
    rs = {'1/4': .25, '1/3': 1 / 3, '1/2': .5, '2/3': 2 / 3, '3/4': .75, '1': 1., '4/3': 4 / 3, '3/2': 1.5, '2': 2.,
          '5/2': 2.5, '3': 3., 'e': math.e, 'pi': math.pi, 'phi': (1 + 5 ** .5) / 2, '1/phi': 2 / (1 + 5 ** .5),
          'sqrt2': 2 ** .5, 'e^gamma': math.exp(.5772156649015329), '1/e': 1 / math.e}
    for nm, r in rs.items():
        for mode, op in (('floor', np.floor), ('round', np.rint)):
            H.append((f'{mode}(a+(a-b)*{nm})', lambda r=r, op=op: op(af + (af - bf) * r)))
            H.append((f'{mode}(a+b*{nm})', lambda r=r, op=op: op(af + bf * r)))
            H.append((f'{mode}(a*{nm}+b)', lambda r=r, op=op: op(af * r + bf)))
    for mode, op in (('floor', np.floor), ('round', np.rint)):
        H.append((f'{mode}(a^(ln a/ln b))', lambda op=op: np.where(b > 1, op(np.exp(np.log(np.maximum(af, 1)) ** 2 / np.log(np.maximum(bf, 2)))), -1)))
    H += [('a xor b', lambda: a ^ b), ('a or b', lambda: a | b), ('a+(a and b)', lambda: a + (a & b)),
          ('a+(a xor b)', lambda: a + (a ^ b)),
          ('bincat(a,b)', lambda: np.where(bitlen(a) + bitlen(b) < 53, (a << bitlen(b)) | b, -1)),
          ('bincat(b,a)', lambda: np.where(bitlen(a) + bitlen(b) < 53, (b << bitlen(a)) | a, -1)),
          ('deccat(a,b)', lambda: np.where(ndig(a) + ndig(b) < 16, a * 10 ** ndig(b) + b, -1)),
          ('deccat(b,a)', lambda: np.where(ndig(a) + ndig(b) < 16, b * 10 ** ndig(a) + a, -1)),
          ('a+revdec(b)', lambda: a + digit_op(b, 'rev')), ('a+revbin(b)', lambda: a + revbin(b)),
          ('a+digitsum(b)', lambda: a + digit_op(b, 'sum')), ('a+popcount(b)', lambda: a + popcount(b)),
          ('a+bitlen(b)', lambda: a + bitlen(b)), ('a+gcd(a,b)', lambda: a + np.gcd(a, b)),
          ('a+b-gcd(a,b)', lambda: a + b - np.gcd(a, b)), ('a+(a mod b)', lambda: a + a % np.maximum(b, 1)),
          ('a+(b mod 10)', lambda: a + b % 10), ('2^bitlen(a)+b', lambda: (1 << bitlen(a)) + b)]
    for m in range(2, 31):
        H.append((f'next y>a, y≡b (mod {m})', lambda m=m: a + 1 + ((b - a - 1) % m)))
    for nm, fn in H:
        with np.errstate(all='ignore'):
            Z = clean(fn(), np.minimum(a, b))
        score(Z, lambda nm=nm: nm)
        for j in ADJ:
            score(ev(j, Z), lambda nm=nm, j=j: f'{names[j]} o [{nm}]')
pickle.dump(dict(stage=STAGE, G=G, best=best, oos=oos, blocks=blocks, realbest=realbest, realoos=realoos),
            open(f'res_{STAGE}.pkl', 'wb'))
print(STAGE, 'G =', G, 'time %.1f s' % (time.time() - t0))
for nm, lo, hi in blocks:
    print(f'  {nm:6s} exact real {best["exact"][lo]:3d} null median {np.median(best["exact"][lo + 1:hi]):4.1f} | '
          f'oos real best_in {oos["best_in"][lo]} max_oos {oos["max"][lo]} | best: {realbest[nm][1][:2]}')
