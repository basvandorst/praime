"""Round 135, after the forty-third review (not pre-registered): the statistics of r135.py on other windows and block
lengths, and the column part from a closed formula with no sieve model:
  E col/(C-1) = tr[ Mt^+ (I - P') K (I - P) ] / (C - 1),   Mt = sum_a (diag P_a - P_a P_a'),
  K = v_q on the diagonal, -Lambda(m)/(phi(m) lam) off it (m = q/gcd(b - a, q)), P = the primes' own transition matrix.
usage: r135b.py NAME X0 [X0 ...]   (windows [X0, X0 + 1e8); block lengths 1e4, 1e5, 1e6) -> ../../data/study/r135b_NAME.json"""
import os, sys, json, math
import numpy as np
H0 = os.path.dirname(os.path.abspath(__file__)); NAME = sys.argv[1]; WIN = [int(float(a)) for a in sys.argv[2:]]
SPAN = 10 ** 8; EXTRA = 10 ** 4; BLS = (10 ** 4, 10 ** 5, 10 ** 6); GAM = 0.5772156649015329; B = 1 - GAM - math.log(2 * math.pi)
def base_primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]
BASE = base_primes(int(math.isqrt(max(WIN) + SPAN + EXTRA)) + 1)
def primes_in(lo, hi):
    seg = np.ones(hi - lo, bool)
    for p, st in zip(BASE.tolist(), ((-lo) % BASE).tolist()): seg[st::p] = False
    return lo + np.nonzero(seg)[0]
def red(q): return [a for a in range(q) if math.gcd(a, q) == 1]
def pf(q): return [p for p in (2, 3, 5, 7, 11, 13) if q % p == 0]
def mangoldt(m):
    f = pf(m)
    if len(f) != 1: return 0.0
    while m % f[0] == 0: m //= f[0]
    return math.log(f[0]) if m == 1 else 0.0
def phi(m): return len(red(m)) if m > 1 else 1
def vq(q, H, lam): return 1 - (math.log(H / q) + sum(math.log(p) / (p - 1) for p in pf(q)) - B) / lam
def s1_law(q, H, lam): return vq(q, H, lam) + (math.log(q) - sum(math.log(p) / (p - 1) for p in pf(q))) / ((phi(q) - 1) * lam)
def col_law(q, P, H, lam):
    r = red(q); C = len(r); K = np.empty((C, C))
    for i, a in enumerate(r):
        for j, b in enumerate(r):
            m = q // math.gcd((b - a) % q, q) if a != b else 1
            K[i, j] = vq(q, H, lam) if a == b else -mangoldt(m) / (phi(m) * lam)
    Mt = sum(np.diag(P[i]) - np.outer(P[i], P[i]) for i in range(C)); I = np.eye(C)
    return float(np.trace(np.linalg.pinv(Mt) @ (I - P.T) @ K @ (I - P)) / (C - 1))
S1Q = (3, 4, 5, 7, 210); S2 = ((3, 1), (3, 2), (5, 1), (7, 1), (10, 1), (10, 2))
acc = {}
def put(key, vals, law): acc.setdefault(key, {'v': [], 'law': []}); acc[key]['v'].append(np.asarray(vals)); acc[key]['law'].append(law)
for x0 in WIN:
    P = np.concatenate([primes_in(x0 + b * 10 ** 6, x0 + (b + 1) * 10 ** 6) for b in range(SPAN // 10 ** 6)] + [primes_in(x0 + SPAN, x0 + SPAN + EXTRA)])
    lam = math.log(x0 + SPAN / 2)
    for BL in BLS:
        NB = SPAN // BL; blk = (P - x0) // BL
        for q in S1Q:
            r = red(q); C = len(r); lut = np.full(q, -1); lut[r] = np.arange(C); c = lut[P % q]; ok = (blk < NB) & (c >= 0)
            n = np.zeros((NB, C)); np.add.at(n, (blk[ok], c[ok]), 1); e = n.sum(1, keepdims=True) / C
            put('S1_%d_%d' % (q, BL), ((n - e) ** 2 / e).sum(1) / (C - 1), s1_law(q, BL, lam))
        for q, k in S2:
            r = red(q); C = len(r); lut = np.full(q, -1); lut[r] = np.arange(C); c = lut[P % q]
            m = len(P) - k; a, bb, bl = c[:m], c[k:], blk[:m]; ok = (bl < NB) & (a >= 0) & (bb >= 0)
            X = np.zeros((NB, C, C)); np.add.at(X, (bl[ok], a[ok], bb[ok]), 1)
            Na = X.sum(2); Ph = X.sum(0) / X.sum(0).sum(1, keepdims=True); E = Na[:, :, None] * Ph[None]
            chi = ((X - E) ** 2 / E).sum((1, 2)); cv = (X - E).sum(1)
            M = np.einsum('ba,aij->bij', Na, np.array([np.diag(Ph[i]) - np.outer(Ph[i], Ph[i]) for i in range(C)]))
            col = np.einsum('bi,bij,bj->b', cv, np.linalg.pinv(M), cv); f = NB / (NB - 1.0)
            put('col_%d_%d_%d' % (q, k, BL), col / (C - 1) * f, col_law(q, Ph, BL, lam))
            put('win_%d_%d_%d' % (q, k, BL), (chi - col) / ((C - 1) ** 2) * f, None)
    print('window %.3g done' % x0, flush=True)
out = {'windows': WIN}
for key, d in acc.items():
    v = np.concatenate(d['v']); law = None if d['law'][0] is None else float(np.mean(d['law']))
    pw = [float(t.mean()) for t in d['v']]; m = float(np.mean(pw)); e = float(v.std() / math.sqrt(len(v)))
    out[key] = {'mean': m, 'err': e, 'law': law, 'per_window': pw}
for BL in BLS:
    print('H = %.0e' % BL)
    print('  S1   ' + '  '.join('q=%d %.4f+-%.4f (law %.4f, z %+.1f)' % (q, out['S1_%d_%d' % (q, BL)]['mean'], out['S1_%d_%d' % (q, BL)]['err'], out['S1_%d_%d' % (q, BL)]['law'], (out['S1_%d_%d' % (q, BL)]['mean'] - out['S1_%d_%d' % (q, BL)]['law']) / out['S1_%d_%d' % (q, BL)]['err']) for q in S1Q))
    print('  col  ' + '  '.join('(%d,%d) %.4f+-%.4f (formula %.4f, z %+.1f)' % (q, k, out['col_%d_%d_%d' % (q, k, BL)]['mean'], out['col_%d_%d_%d' % (q, k, BL)]['err'], out['col_%d_%d_%d' % (q, k, BL)]['law'], (out['col_%d_%d_%d' % (q, k, BL)]['mean'] - out['col_%d_%d_%d' % (q, k, BL)]['law']) / out['col_%d_%d_%d' % (q, k, BL)]['err']) for q, k in S2))
    print('  within ' + '  '.join('(%d,%d) %.4f+-%.4f' % (q, k, out['win_%d_%d_%d' % (q, k, BL)]['mean'], out['win_%d_%d_%d' % (q, k, BL)]['err']) for q, k in S2))
json.dump(out, open(os.path.join(H0, '..', '..', 'data', 'study', 'r135b_%s.json' % NAME), 'w'), indent=0)
