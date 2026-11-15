"""Theory checks on the streamed runs (runs_100_12000.pkl), plus an independent
recomputation of the late identities with plain Python integers.
Output: theory.txt
"""
import math, pickle, random
import numpy as np
from scipy.stats import kstest

EXP = {107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213}
R = pickle.load(open('runs_100_12000.pkl', 'rb'))
out = open('theory.txt', 'w')


def P(*a):
    s = ' '.join(str(x) for x in a)
    print(s); out.write(s + '\n')


ks = [r['k'] for r in R]
P('runs:', len(R), 'survivors 100 <= k < 12000; exponents among them:', sorted(EXP & set(ks)))

# 1. first k-specific iteration: i0 = first i whose unreduced value reaches M
P('\n1. First reduction i0 (before it, s_i and x_i are the same integers for every k)')
bad = 0
for r in R:
    k = r['k']; M = (1 << k) - 1
    S, i = 4, 0
    while S < M:
        S = S * S - 2; i += 1
    X, j = 3, 0
    while X < M:
        X = X * X; j += 1
    bad += (i != r['LL']['i0']) + (j != r['PRP']['i0'])
    r['_f'] = (math.ceil(math.log2(k / math.log2(2 + math.sqrt(3)))), math.ceil(math.log2(k / math.log2(3))))
P('   i0 recomputed with plain integers: %d mismatches' % bad)
P('   LL  i0 range %d..%d; formula ceil(log2(k/log2(2+sqrt3))) differs in %d cases' % (
    min(r['LL']['i0'] for r in R), max(r['LL']['i0'] for r in R), sum(r['_f'][0] != r['LL']['i0'] for r in R)))
P('   PRP i0 range %d..%d; formula ceil(log2(k/log2 3)) differs in %d cases' % (
    min(r['PRP']['i0'] for r in R), max(r['PRP']['i0'] for r in R), sum(r['_f'][1] != r['PRP']['i0'] for r in R)))

# 2. late identities (from the streamed runs)
P('\n2. Late identities (streamed runs, all survivors)')
z = {r['k'] for r in R if r['LL']['s_k2_zero']}
s3 = {r['k']: r['LL']['s_k3_sign'] for r in R if r['LL']['s_k3_sign'] != 0}
x1 = {r['k'] for r in R if r['PRP']['x_k1_is_m3']}
x9 = {r['k'] for r in R if r['PRP']['x_k_is_9']}
ex = EXP & set(ks)
P('   s_{k-2} = 0 exactly for the exponents:', z == ex)
P('   s_{k-3} = +-2^((k+1)/2) exactly for the exponents:', set(s3) == ex, ' signs:', dict(sorted(s3.items())))
P('   x_{k-1} = -3 exactly for the exponents:', x1 == ex, '| x_k = 9 exactly for the exponents (no base-3 Fermat pseudoprime):', x9 == ex)
P('   late Hamming weights of s_i for the exponents (i = k-7 .. k-2), as z = (w - k/2)/(sqrt(k)/2):')
for r in R:
    if r['k'] in ex:
        k = r['k']; L = r['LL']['late_hw']
        P('     k=%5d ' % k + ' '.join('%+.1f' % ((L[i] - k / 2) / (math.sqrt(k) / 2)) for i in sorted(L)))
nz = [((r['LL']['late_hw'][r['k'] - 3] - r['k'] / 2) / (math.sqrt(r['k']) / 2)) for r in R if r['k'] not in ex]
P('   non-exponents, z of s_{k-3}: min %.1f max %.1f (an exponent has |z| = sqrt(k) - 2/sqrt(k))' % (min(nz), max(nz)))

# 3. Jacobi symbols at every checkpoint
P('\n3. Jacobi symbols at the 8 checkpoints of every run')
jp = set(); jm = set(); js = []
for r in R:
    for f in r['LL']['snaps'].values():
        jp.add(f['chk_jac_s+2']); jm.add(f['chk_jac_s-2']); js.append(f['v_jac'])
P('   (s_i + 2 | M) values seen:', jp, ' (s_i - 2 | M) values seen:', jm, ' (s_i | M) mean %.3f over %d' % (np.mean(js), len(js)))
P('   PRP (x_i | M) values seen:', {f['v_jac'] for r in R for f in r['PRP']['snaps'].values()})

# 4. repeats and Brent / Pollard-rho signals
P('\n4. Repeats and rho (Brent) signals over the whole runs')
P('   LL exact repeats:', [r['k'] for r in R if r['LL']['repeat_at']])
P('   PRP exact repeats (k, i):', [(r['k'], r['PRP']['repeat_at']) for r in R if r['PRP']['repeat_at']])
P('   factors found by Brent gcd within the run: LL', [(r['k'], r['LL']['brent_first'], r['LL']['final_brent_gcd']) for r in R if r['LL']['final_brent_gcd'] > 1],
  ' PRP', [(r['k'], r['PRP']['brent_first'], r['PRP']['final_brent_gcd']) for r in R if r['PRP']['final_brent_gcd'] > 1])

# 5. calibration: every z-type aggregate should be ~N(0,1) over all survivors if the runs are pseudo-random
P('\n5. Calibration of the aggregate z-scores over all survivors (horizon k/2): mean, sd, KS p vs N(0,1)')
names = [n for n in R[0]['LL']['snaps']['k/2'] if n.startswith('a_') and n.endswith('_z') or n in ('a_hw_meanz', 'a_run_meanz', 'a_hw_lag1z')]
worst = 1
for seq in ('LL', 'PRP'):
    for n in names:
        v = np.array([r[seq]['snaps']['k/2'][n] for r in R])
        p = kstest(v, 'norm').pvalue
        worst = min(worst, p)
        P('   %-3s %-14s mean %+.3f sd %.3f KS p %.3f' % (seq, n, v.mean(), v.std(), p))
P('   smallest KS p over %d tests: %.3f' % (2 * len(names), worst))

# 6. independent recheck of the late identities with plain Python integers (% reduction)
P('\n6. Independent recheck (plain Python ints): s_{k-3}, s_{k-2}, x_{k-1} for the exponents and 5 random survivors, and 19937')
random.seed(1)
sample = sorted(ex) + random.sample([k for k in ks if k not in ex], 5) + [19937]
for k in sample:
    M = (1 << k) - 1
    s = 4
    for _ in range(k - 3):
        s = (s * s - 2) % M
    s3v = s; s2v = (s * s - 2) % M
    x = 3
    for _ in range(k - 1):
        x = x * x % M
    r2 = pow(2, (k + 1) // 2, M)
    P('   k=%5d s_{k-2}=0: %-5s s_{k-3}: %-2s x_{k-1}=-3: %s' % (
        k, s2v == 0, '+' if s3v == r2 else ('-' if s3v == M - r2 else 'no'), x == M - 3))
out.close()
