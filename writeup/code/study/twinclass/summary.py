"""Markdown tables of the retrospective (post-fix numbers) -> summary.md"""
import json, math
import numpy as np
D = json.load(open('all_decomp.json')); L = json.load(open('levels.json')); CH = json.load(open('cross_h.json'))
o = []
o.append('| window | Q | u2X meas / pred (ratio) | u2Z meas / pred | u2R meas / pred (ratio) | beta (nominal se) |')
o.append('|---|---|---|---|---|---|')
for w in ('1e7-1e8', '1e8-1e9', '1e9-1e10', 'A', 'B', 'C'):
    for Q in (210, 2310, 30030):
        e = D[w]['h2_%d' % Q]
        o.append('| %s | %d | %.3f / %.3f (%.2f) | %.3f / %.3f | %.3f / %.3f (%.2f) | %.2f (%.2f) |' % (w, Q, e['measured_u2X'], e['pred_u2X'], e['ratio_X'],
                 e['measured_u2Z'], e['pred_u2Z'], e['measured_u2R'], e['pred_u2R'], e['ratio_R'], e['beta_hat'], e['se_beta_approx']))
o.append('')
o.append('| window | X: level 1 (210) | level 2 (2310 in 210) | level 3 (30030 in 2310) | R: level 1 | level 2 | level 3 |')
o.append('|---|---|---|---|---|---|---|')
for w in ('1e7-1e8', '1e8-1e9', '1e9-1e10', 'A', 'B', 'C'):
    x, r = L['levels'][w]['X'], L['levels'][w]['R']
    o.append('| %s | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |' % ((w,) + tuple(x) + tuple(r)))
o.append('| nominal se per window | 0.378 | 0.129 | 0.038 | 0.378 | 0.129 | 0.038 |')
# pooled beta 30030 over the 6 whole windows and over the 27 sub-windows
bw = np.array([D[w]['h2_30030']['beta_hat'] for w in ('1e7-1e8', '1e8-1e9', '1e9-1e10', 'A', 'B', 'C')])
sw = np.array([D[w]['h2_30030']['se_beta_approx'] for w in ('1e7-1e8', '1e8-1e9', '1e9-1e10', 'A', 'B', 'C')])
wt = 1 / sw ** 2
o.append('')
o.append('beta (mod 30030), six whole windows, inverse-variance: %.4f +- %.4f' % (np.sum(wt * bw) / wt.sum(), 1 / math.sqrt(wt.sum())))
sc = L['scatter']
o.append('27 sub-windows: beta %.4f (sd %.4f, nominal se %.4f); R levels %.3f +- %.3f, %.3f +- %.3f, %.3f +- %.3f; X levels %.3f, %.3f, %.3f' % (
    sc['beta30030']['mean'], sc['beta30030']['sd'], sc['beta30030']['se_nominal'], sc['R_level1']['mean'], sc['R_level1']['se_nominal'],
    sc['R_level2']['mean'], sc['R_level2']['se_nominal'], sc['R_level3']['mean'], sc['R_level3']['se_nominal'],
    sc['X_level1']['mean'], sc['X_level2']['mean'], sc['X_level3']['mean']))
o.append('sd/nominal of the 27 sub-window ratios: X %.2f %.2f %.2f; R %.2f %.2f %.2f' % tuple(
    [sc['X_level%d' % i]['sd_over_nominal'] for i in (1, 2, 3)] + [sc['R_level%d' % i]['sd_over_nominal'] for i in (1, 2, 3)]))
o.append('')
o.append('| pair gap h | cross z of R, mod 210 (3 decades; joint) | mod 2310 | mod 30030 |')
o.append('|---|---|---|---|')
for h in (2, 4, 6, 8, 10, 12):
    cells = []
    for Q in (210, 2310, 30030):
        zs = [CH['%s|h%d_%d' % (d, h, Q)]['R']['z'] for d in ('1e7-1e8', '1e8-1e9', '1e9-1e10')]
        cells.append('%s; %.2f' % (' '.join('%+.2f' % z for z in zs), sum(zs) / math.sqrt(3)))
    o.append('| %d | %s |' % (h, ' | '.join(cells)))
oth = [sum(CH['%s|h%d_%d' % (d, h, Q)]['R']['z'] for d in ('1e7-1e8', '1e8-1e9', '1e9-1e10')) / math.sqrt(3) for h in (4, 6, 8, 10, 12) for Q in (210, 2310, 30030)]
o.append('')
o.append('h = 4..12: 15 joint z, mean %.2f, sd %.2f' % (np.mean(oth), np.std(oth, ddof=1)))
for Q in (210, 2310, 30030):
    c = [CH['%s|h2_%d' % (d, Q)]['R'] for d in ('1e7-1e8', '1e8-1e9', '1e9-1e10')]
    o.append('twins Q=%d: cross meas %s, HL %s, sd %s' % (Q, ' '.join('%.3f' % x['cross_meas'] for x in c), ' '.join('%.3f' % x['cross_HL'] for x in c), ' '.join('%.3f' % x['sd'] for x in c)))
open('summary.md', 'w').write('\n'.join(o) + '\n')
print('\n'.join(o))
