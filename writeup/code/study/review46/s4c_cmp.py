import sys, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__file__))
from mylib import OUT
F = np.load(OUT + '/s4_full_100.npz'); A = np.load(OUT + '/s4_asym.npz'); S2 = np.load(OUT + '/s2_pairs.npz')
sub = F['sub']; near = F['near']; fr = S2['SW'][:, near[sub]].T
a3 = A['a3_200k'][sub]; a3b = A['a3_50k'][sub]
ref = F['f200k_du02']
def rel(x, y): return np.sqrt(((x - y) ** 2).sum() / (y ** 2).sum())
print('100 pairs (50 largest frozen + 50 random); relative rms difference:')
for k in ('f10k_du02', 'f50k_du02'):
    print('  full integral %-10s vs full 2e5: %.3f' % (k, rel(F[k], ref)))
print('  frozen (5e4, du 0.06) vs full 2e5: %.3f' % rel(fr, ref))
print('  asym L^-3 5e4 vs asym 2e5: %.4f' % rel(a3b, a3))
print('  asym 2e5 vs full 2e5: %.3f; vs frozen: %.3f' % (rel(a3, ref), rel(a3, fr)))
res = ref - a3; print('  resonant part (full 2e5 - asym) rms per w: ' + ' '.join('%.4f' % v for v in np.sqrt((res ** 2).mean(0))))
print('  asym rms per w:                           ' + ' '.join('%.4f' % v for v in np.sqrt((a3 ** 2).mean(0))))
res5 = F['f50k_du02'] - a3; print('  corr over pairs and w of the resonant parts at 5e4 and 2e5: %.3f' % np.corrcoef(res5.ravel(), res.ravel())[0, 1])
