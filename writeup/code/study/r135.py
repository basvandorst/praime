"""Round 135 (preregister/round135.txt): primes in W1-W3 = [3e11, 6e11, 9e11] + 1e8, 100 blocks of 1e6 each.
S1: equidistribution chi2 of the block's primes over the reduced classes mod q (q = 3, 4, 10, 210).
S2: block homogeneity of the transition tables (p_n, p_{n+k}) mod q, split into a column part and a within part.
The primes are listed with a plain segmented sieve (data, not a model). usage: r135.py"""
import os, json, math
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(H, '..', '..', 'data', 'study', 'r135.json')
WIN = [3 * 10 ** 11, 6 * 10 ** 11, 9 * 10 ** 11]; SPAN = 10 ** 8; BL = 10 ** 6; NB = 100; EXTRA = 10 ** 4
def base_primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]
BASE = base_primes(int(math.isqrt(WIN[-1] + SPAN + EXTRA)) + 1)
def primes_in(lo, hi):
    seg = np.ones(hi - lo, bool)
    for p, st in zip(BASE.tolist(), ((-lo) % BASE).tolist()): seg[st::p] = False
    return lo + np.nonzero(seg)[0]
def red(q): return [a for a in range(q) if math.gcd(a, q) == 1]
S1Q = (3, 4, 10, 210); S2 = ((3, 1), (3, 2), (10, 1), (10, 2))
LAW1 = {3: 0.477, 4: 0.488, 10: 0.496, 210: 0.564}; SIEVE1 = 0.807
LAW2 = {(3, 1): (0.763, 0.564, 0.962, 0.971), (3, 2): (0.739, 0.486, 0.993, 0.917), (10, 1): (0.872, 0.578, 0.970, 0.980), (10, 2): (0.881, 0.509, 1.005, 0.968)}
SP = {3: 0.010, 10: 0.004}
s1 = {q: [] for q in S1Q}; s2 = {k: {'chi': [], 'col': [], 'win': []} for k in S2}; nprimes = []
for x0 in WIN:
    P = np.concatenate([primes_in(x0 + b * BL, x0 + (b + 1) * BL) for b in range(NB)] + [primes_in(x0 + SPAN, x0 + SPAN + EXTRA)])
    blk = (P - x0) // BL; nprimes.append(int((blk < NB).sum()))
    for q in S1Q:
        r = red(q); lut = np.full(q, -1); lut[r] = np.arange(len(r)); c = lut[P % q]
        for b in range(NB):
            n = np.bincount(c[blk == b], minlength=len(r)).astype(float); e = n.sum() / len(r)
            s1[q].append(float(((n - e) ** 2).sum() / e / (len(r) - 1)))
    for q, k in S2:
        r = red(q); C = len(r); lut = np.full(q, -1); lut[r] = np.arange(C); c = lut[P % q]
        m = len(P) - k; a, bb, bl = c[:m], c[k:], blk[:m]; ok = bl < NB
        X = np.zeros((NB, C, C)); np.add.at(X, (bl[ok], a[ok], bb[ok]), 1)
        Na = X.sum(2); Ph = X.sum(0) / X.sum(0).sum(1, keepdims=True)
        for b in range(NB):
            E = Na[b][:, None] * Ph; chi = float(((X[b] - E) ** 2 / E).sum())
            cv = (X[b] - E).sum(0); M = sum(Na[b][i] * (np.diag(Ph[i]) - np.outer(Ph[i], Ph[i])) for i in range(C))
            col = float(cv @ np.linalg.pinv(M) @ cv)
            s2[(q, k)]['chi'].append(chi / (C * (C - 1))); s2[(q, k)]['col'].append(col / (C - 1)); s2[(q, k)]['win'].append((chi - col) / ((C - 1) * (C - 1)))
    print('window %.0e: %d primes' % (x0, nprimes[-1]), flush=True)
res = {'primes': nprimes, 'S1': {}, 'S2': {}}; R = {}
print('S1: q, R1 per window, pooled +- error, law, z')
for q in S1Q:
    v = np.array(s1[q]); pw = [float(v[i * NB:(i + 1) * NB].mean()) for i in range(3)]; m = float(np.mean(pw)); e = float(v.std() / math.sqrt(len(v)))
    res['S1'][q] = {'per_window': pw, 'pooled': m, 'err': e, 'law': LAW1[q]}
    R['1_q%d' % q] = bool(abs(m - LAW1[q]) <= 3 * e)
    print('  %3d  %s  %.4f +- %.4f  law %.3f  z %+.2f' % (q, ' '.join('%.3f' % t for t in pw), m, e, LAW1[q], (m - LAW1[q]) / e))
R['1_q210_below_sieve'] = bool(res['S1'][210]['pooled'] < SIEVE1 - 5 * res['S1'][210]['err'])
print('S2: (q,k), R2 / column / within: pooled +- error, law, z; sieve model R2')
for q, k in S2:
    d = {}; f = NB / (NB - 1.0)
    for nm, law in zip(('chi', 'col', 'win'), LAW2[(q, k)][:3]):
        v = np.array(s2[(q, k)][nm]); pw = [float(v[i * NB:(i + 1) * NB].mean() * f) for i in range(3)]
        m = float(np.mean(pw)); em = float(v.std() / math.sqrt(len(v)) * f); sg = math.hypot(em, SP[q])
        d[nm] = {'per_window': pw, 'pooled': m, 'err': em, 'sigma': sg, 'law': law, 'z': (m - law) / sg}
    d['sieve'] = LAW2[(q, k)][3]; res['S2']['%d_%d' % (q, k)] = d
    R['2_law_%d_%d' % (q, k)] = bool(abs(d['chi']['z']) <= 3)
    if k == 1:
        R['2_below_sieve_%d_1' % q] = bool(d['chi']['pooled'] < d['sieve'] - 3 * d['chi']['sigma'])
        R['3_col_%d_1' % q] = bool(abs(d['col']['z']) <= 3); R['3_within_%d_1' % q] = bool(abs(d['win']['z']) <= 3)
    print('  (%d,%d)  ' % (q, k) + '  '.join('%s %.4f +- %.4f law %.3f z %+.2f' % (nm, d[nm]['pooled'], d[nm]['sigma'], d[nm]['law'], d[nm]['z']) for nm in ('chi', 'col', 'win')) + '  sieve %.3f' % d['sieve'])
w1, w2 = res['S2']['10_1'], res['S2']['10_2']
R['4_falsified'] = bool(w1['win']['z'] < -3 or w2['win']['z'] < -3 or abs(w1['chi']['pooled'] - w1['sieve']) <= 2 * w1['chi']['sigma'])
res['readings'] = R; print('readings', R)
json.dump(res, open(OUT, 'w'), indent=0)
