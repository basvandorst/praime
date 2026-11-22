"""Stress test B, part 2: count (g1, g2) for primes p_n in [1e9, 1.1e9), compare with theory.npz (written earlier)."""
import sys, math, json
sys.path.insert(0, '/home/bas/praime/writeup/code/study')
import numpy as np
import common as C
from scipy.stats import chi2 as CHI2

L0 = 10 ** 9; W = 10 ** 8
P = np.nonzero(C.primes_mask(L0, W + 2000))[0].astype(np.int64) + L0      # a little beyond, so every p_n < 1.1e9 has two successors
n = int(np.sum(P < L0 + W))
g1 = (P[1:n + 1] - P[:n]); g2 = (P[2:n + 2] - P[1:n + 1])
assert len(g2) == n
H = np.zeros((19, 19), np.int64)                 # H[g1/2, g2/2] for g <= 36
m = (g1 <= 36) & (g2 <= 36)
np.add.at(H, (g1[m] // 2, g2[m] // 2), 1)
H = H[1:, 1:].astype(float)                       # index g/2 - 1
th = np.load('/home/bas/scratch/stress-B/theory.npz')
BINS = [(2, 4), (6, 8), (10, 12), (14, 18), (20, 26), (28, 36)]
def binof(g):
    for i, (lo, hi) in enumerate(BINS):
        if lo <= g <= hi:
            return i
gs = np.arange(2, 37, 2)
def cells(A):
    Cc = np.zeros((6, 6))
    for a in gs:
        for b in gs:
            Cc[binof(a), binof(b)] += A[a // 2 - 1, b // 2 - 1]
    return Cc
def mod6(A):
    Cc = np.zeros((3, 3))
    for a in gs:
        for b in gs:
            Cc[(a % 6) // 2, (b % 6) // 2] += A[a // 2 - 1, b // 2 - 1]
    return Cc
out = {'N': n, 'mass_data': float(H.sum() / n)}
obs_s = cells(H); obs_6 = mod6(H)
print('N primes', n, 'mass g1,g2<=36 data %.6f' % (H.sum() / n))
print('DATA size cells, prob x1e4'); print(np.array2string(1e4 * obs_s / n, precision=4, max_line_width=200))
print('DATA mod6 x1e4'); print(np.array2string(1e4 * obs_6 / n, precision=4))
for name in ('T', 'Tc', 'Ti'):
    A = th[name][:18, :18]
    ex_s = n * cells(A); ex_6 = n * mod6(A)
    zs = (obs_s - ex_s) / np.sqrt(ex_s)
    c_s = float(np.sum(zs ** 2))
    k6 = ex_6 > 0
    z6 = np.zeros((3, 3)); z6[k6] = (obs_6[k6] - ex_6[k6]) / np.sqrt(ex_6[k6])
    c_6 = float(np.sum(z6 ** 2)); df6 = int(k6.sum())
    # full grid 18x18 (exploratory, not preregistered): cells with expectation > 0
    Ef = n * A; kf = Ef > 50
    c_f = float(np.sum((H[kf] - Ef[kf]) ** 2 / Ef[kf])); df_f = int(kf.sum())
    imp = float(H[~kf].sum())
    # the "shape" test: conditional on the total mass in the 36 cells (removes an overall normalisation error)
    ex_sn = ex_s * obs_s.sum() / ex_s.sum()
    c_sn = float(np.sum((obs_s - ex_sn) ** 2 / ex_sn))
    tot_z = float((obs_s.sum() - ex_s.sum()) / math.sqrt(ex_s.sum() * (1 - ex_s.sum() / n)))
    out[name] = {'chi2_size36': c_s, 'p_size36': float(CHI2.sf(c_s, 36)), 'chi2_mod6': c_6, 'df_mod6': df6, 'p_mod6': float(CHI2.sf(c_6, df6)),
                 'chi2_shape_df35': c_sn, 'p_shape': float(CHI2.sf(c_sn, 35)), 'z_total_mass': tot_z,
                 'chi2_full_grid': c_f, 'df_full_grid': df_f, 'p_full': float(CHI2.sf(c_f, df_f)), 'obs_in_cells_with_exp<=50': imp,
                 'z_size_cells': np.round(zs, 2).tolist(), 'z_mod6': np.round(z6, 2).tolist()}
    print('\n==', name, json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in out[name].items() if not k.startswith('z_') or k == 'z_total_mass'}))
    print('z size cells'); print(np.round(zs, 1))
    print('z mod6'); print(np.round(z6, 1))
# conditional (b)-type statistic
for name in ('T', 'Tc', 'Ti'):
    M6 = mod6(th[name][:18, :18])
    pr = M6[0, 0] / M6[0].sum()
    ob = obs_6[0, 0] / obs_6[0].sum(); se = math.sqrt(ob * (1 - ob) / obs_6[0].sum())
    out[name]['cond_g2_6_given_g1_6'] = {'pred': pr, 'obs': ob, 'se': se, 'z': (pr - ob) / se}
    print(name, 'P(6|g2 | 6|g1) pred %.5f obs %.5f se %.5f z %+.2f' % (pr, ob, se, (pr - ob) / se))
for a, b in ((2, 4), (4, 2), (2, 10), (6, 6), (12, 12), (30, 6), (6, 30), (4, 8)):
    o = H[a // 2 - 1, b // 2 - 1]
    print('cell (%d,%d) obs %d  law %.1f (z %+.1f) chain %.1f (z %+.1f)' % (a, b, o, n * th['T'][a // 2 - 1, b // 2 - 1], (o - n * th['T'][a // 2 - 1, b // 2 - 1]) / math.sqrt(n * th['T'][a // 2 - 1, b // 2 - 1]),
          n * th['Tc'][a // 2 - 1, b // 2 - 1], (o - n * th['Tc'][a // 2 - 1, b // 2 - 1]) / math.sqrt(n * th['Tc'][a // 2 - 1, b // 2 - 1])))
np.save('/home/bas/scratch/stress-B/H.npy', H)
json.dump(out, open('/home/bas/scratch/stress-B/result.json', 'w'), indent=1)
