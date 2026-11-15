"""Score every map of the grammar on the real sequence (row 0) and on all null sequences at once.

usage: search.py STAGE  with STAGE in d1 d2 d3 tt  (forward)  or  inv-d1 inv-d2 inv-d3 (inverse)
  d1   : every primitive f, a_{n+1} =? f(a_n)
  d2   : g(f(x)) for f, g in the 'core2' set, plus f in all primitives with g in a small adjust set
  d3   : h(g(f(x))) for f, g, h in the 'core3' set
  tt   : two-term rules a_{n+1} =? h(a_n, a_{n-1}), and adjust(h(...))
  inv-*: the same grammar on the reversed pairs: a_n =? f(a_{n+1})
Metrics per (map, sequence), over the T transitions:
  exact  = #transitions hit exactly;   run0 = consecutive exact hits from the first transition;
  runany = longest run of consecutive exact hits;  pm1, pm2 = #hits within +-1, +-2;
  rel1 = #hits within 1% (|f - y| <= 0.01 y);
  OOS: in-sample = transitions predicting a_2..a_25 (a_3..a_25 for tt), out-of-sample = the rest.
Output: res_STAGE.pkl with per-row maxima, OOS tie statistics, G (#maps) and the real-row candidates.
"""
import sys, time, pickle, numpy as np
from grammar import build, apply, VMAX

STAGE = sys.argv[1]
INV = STAGE.startswith('inv-')
D = np.load('seqs.npz')
M = np.vstack([D['real'][None, :], D['nullA'], D['nullB']]).astype(np.int64)
S = M.shape[0]
P = build()
names = [p[0] for p in P]
core2 = [i for i, p in enumerate(P) if 'core2' in p[2]]
core3 = [i for i, p in enumerate(P) if 'core3' in p[2]]
C2SET = set(core2)
EXTRA3 = ('x-3', 'x+3', '3x', '3x+1', '3x-1', '4x+1', '4x-1', '4x+3', '2x+3', '2x-3', 'floor(x/3)', 'round(x*e)',
          'round(x*pi)', 'round(x*sqrt2)', 'round(x*4/3)', 'round(x*5/3)', 'round(x*5/2)', 'floor(x*3/2)',
          'floor(x*phi)', 'x+popcount(x)', 'x+bitlen(x)', 'x+digitsum(x)', 'x+revdec(x)', 'x+revbin(x)',
          'next 2^k+1', 'next 2^k-1', 'next≡1(mod 3)', 'next≡2(mod 3)', 'next≡1(mod 8)', 'next≡7(mod 8)',
          'gray', 'invgray', 'next pop=3', 'next square', 'prepend10', 'prepend11')
core3 = sorted(set(core3) | {i for i, p in enumerate(P) if p[0] in EXTRA3})
ADJ = [i for i, p in enumerate(P) if p[0] in ('x-2', 'x-1', 'x+1', 'x+2', '2x+1', '2x-1', 'revbin', 'revdec',
                                                 'next≡1(mod 4)', 'next≡3(mod 4)', 'next≡1(mod 6)', 'next≡5(mod 6)',
                                                 'next binpal', 'next 2^k-1', 'next 2^k+1', 'nextsamepop', 'next decpal')]

if STAGE == 'tt':
    X = M[:, 1:-1]; X2 = M[:, :-2]; Y = M[:, 2:]
    N = np.broadcast_to(np.arange(2, M.shape[1]), X.shape).astype(np.int64)
    NIN = 23  # predicts a_3..a_25
elif INV:
    X = M[:, 1:]; Y = M[:, :-1]
    N = np.broadcast_to(np.arange(1, M.shape[1]), X.shape).astype(np.int64)
    NIN = 24
else:
    X = M[:, :-1]; Y = M[:, 1:]
    N = np.broadcast_to(np.arange(1, M.shape[1]), X.shape).astype(np.int64)
    NIN = 24
T = X.shape[1]
Yf = Y.astype(np.float64)
TOL = 0.01 * Yf

KEYS = ('exact', 'run0', 'runany', 'pm1', 'pm2', 'rel1')
best = {k: np.zeros(S, np.int64) for k in KEYS}
oos = dict(best_in=np.full(S, -1, np.int64), n=np.zeros(S, np.int64), sum=np.zeros(S, np.float64),
           max=np.zeros(S, np.int64))
cands = []
G = 0
cols = np.arange(T)


def score(F, name_fn):
    global G
    G += 1
    E = F == Y
    ex = E.sum(1)
    d = np.abs(F - Y).astype(np.float64)
    d[F < 0] = np.inf
    pm1 = (d <= 1).sum(1); pm2 = (d <= 2).sum(1); rel = (d <= TOL).sum(1)
    run0 = np.where(E.all(1), T, np.argmin(E, axis=1))
    if ex.max() >= 2:
        c = np.cumsum(E, axis=1)
        reset = np.maximum.accumulate(np.where(E, 0, c), axis=1)
        runany = (c - reset).max(1)
    else:
        runany = ex
    for k, v in (('exact', ex), ('run0', run0), ('runany', runany), ('pm1', pm1), ('pm2', pm2), ('rel1', rel)):
        np.maximum(best[k], v, out=best[k])
    s_in = E[:, :NIN].sum(1); s_out = ex - s_in
    gt = s_in > oos['best_in']; eq = s_in == oos['best_in']
    oos['max'] = np.where(gt, s_out, np.where(eq, np.maximum(oos['max'], s_out), oos['max']))
    oos['sum'] = np.where(gt, s_out, oos['sum'] + eq * s_out)
    oos['n'] = np.where(gt, 1, oos['n'] + eq)
    oos['best_in'] = np.where(gt, s_in, oos['best_in'])
    r = (int(ex[0]), int(run0[0]), int(runany[0]), int(pm1[0]), int(pm2[0]), int(rel[0]), int(s_in[0]), int(s_out[0]))
    if r[0] >= 5 or r[1] >= 4 or r[2] >= 4 or r[3] >= 8 or r[4] >= 9 or r[5] >= 5 or r[6] >= 5:
        hits = [int(i) for i in np.nonzero(E[0])[0]]
        cands.append((name_fn(), r, hits))


def ev(i, x):
    return apply(P[i][1], x, N)


t0 = time.time()
if STAGE in ('d1', 'inv-d1'):
    for i in range(len(P)):
        score(ev(i, X), lambda i=i: names[i])
elif STAGE in ('d2', 'inv-d2'):
    for i in range(len(P)):
        Z = ev(i, X)
        if (Z < 0).all():
            G += len(P) if i in C2SET else len(core2)
            continue
        inner = range(len(P)) if i in C2SET else core2
        for j in inner:
            score(ev(j, Z), lambda i=i, j=j: f'{names[j]} o {names[i]}')
elif STAGE in ('d3', 'inv-d3'):
    for i in core3:
        Z = ev(i, X)
        for j in core3:
            Z2 = ev(j, Z)
            if (Z2 < 0).all():
                G += len(core3); continue
            for k in core3:
                score(ev(k, Z2), lambda i=i, j=j, k=k: f'{names[k]} o {names[j]} o {names[i]}')
elif STAGE == 'tt':
    from grammar import clean, bitlen, popcount, revbin, digit_op, ndig
    import math
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
        H.append((f'{mode}(a^(ln a/ln b))', lambda op=op: np.where(b > 1, op(np.exp(np.log(af) ** 2 / np.log(np.maximum(bf, 2)))), -1)))
    H += [('a xor b', lambda: a ^ b), ('a or b', lambda: a | b), ('a+(a and b)', lambda: a + (a & b)),
          ('a+(a xor b)', lambda: a + (a ^ b)),
          ('bincat(a,b)', lambda: np.where(bitlen(a) + bitlen(b) < 53, (a << bitlen(b)) | b, -1)),
          ('bincat(b,a)', lambda: np.where(bitlen(a) + bitlen(b) < 53, (b << bitlen(a)) | a, -1)),
          ('deccat(a,b)', lambda: np.where(ndig(a) + ndig(b) < 16, a * 10 ** ndig(b) + b, -1)),
          ('deccat(b,a)', lambda: np.where(ndig(a) + ndig(b) < 16, b * 10 ** ndig(a) + a, -1)),
          ('a+revdec(b)', lambda: a + digit_op(b, 'rev')), ('a+revbin(b)', lambda: a + revbin(b)),
          ('a+digitsum(b)', lambda: a + digit_op(b, 'sum')), ('a+popcount(b)', lambda: a + popcount(b)),
          ('a+bitlen(b)', lambda: a + bitlen(b)), ('a+gcd(a,b)', lambda: a + np.gcd(a, b)),
          ('a+b-gcd(a,b)', lambda: a + b - np.gcd(a, b)), ('a+(a mod b)', lambda: a + a % b),
          ('a+(b mod 10)', lambda: a + b % 10), ('2^bitlen(a)+b', lambda: (1 << bitlen(a)) + b)]
    for m in range(2, 31):
        H.append((f'next y>a, y≡b (mod {m})', lambda m=m: a + 1 + ((b - a - 1) % m)))
    for nm, fn in H:
        with np.errstate(all='ignore'):
            Z = clean(fn(), np.minimum(a, b))
        score(Z, lambda nm=nm: nm)
        for j in ADJ:
            score(ev(j, Z), lambda nm=nm, j=j: f'{names[j]} o [{nm}]')
else:
    sys.exit('unknown stage')

pickle.dump(dict(stage=STAGE, G=G, best=best, oos=oos, cands=cands, S=S, T=T, NIN=NIN,
                 nA=len(D['nullA']), nB=len(D['nullB'])), open(f'res_{STAGE}.pkl', 'wb'))
print(STAGE, 'G =', G, 'time %.1f s' % (time.time() - t0), 'cands', len(cands))
for k in KEYS:
    print(f'  {k:7s} real {best[k][0]:3d}   nullA median {np.median(best[k][1:251]):5.1f}   nullB median {np.median(best[k][251:]):5.1f}')
