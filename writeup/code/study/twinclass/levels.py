"""Level split (mod 210 | 2310 within 210 | 30030 within 2310) of the twin spread X and of the remainder R, the
within/cross split of each decade into its nine sub-windows, and the scatter of the sub-window ratios against the
nominal chi^2 noise.  Input all_decomp.json (run_all.py).  -> levels.json"""
import json, math
import numpy as np
D = json.load(open('all_decomp.json'))
K = {210: 15, 2310: 135, 30030: 1485}
def lev(e, key):          # key 'X' or 'R'; returns level ratios and their nominal se
    u = {Q: e['h2_%d' % Q]['measured_u2' + key] for Q in K}; v = {Q: e['h2_%d' % Q]['pred_u2' + key] for Q in K}
    l1 = u[210] / v[210]
    l2 = (134 * u[2310] - 14 * u[210]) / (134 * v[2310] - 14 * v[210])
    l3 = (1484 * u[30030] - 134 * u[2310]) / (1484 * v[30030] - 134 * v[2310])
    return l1, l2, l3
SE = (math.sqrt(2 / 14), math.sqrt(2 / 120), math.sqrt(2 / 1350))
out = {'levels': {}, 'cross': {}, 'scatter': {}}
print('window          | X: level1 level2 level3 | R: level1 level2 level3')
for w, e in D.items():
    lx = lev(e, 'X'); lr = lev(e, 'R')
    out['levels'][w] = {'X': lx, 'R': lr}
    if '/sub' not in w:
        print('%-15s | %.3f %.3f %.3f | %.3f %.3f %.3f' % ((w,) + lx + lr))
for dec in ('1e7-1e8', '1e8-1e9', '1e9-1e10'):
    e = D[dec]; subs = [D['%s/sub%d' % (dec, s)] for s in range(9)]
    for Q in K:
        for key in ('X', 'R'):
            tot_m = e['h2_%d' % Q]['measured_u2' + key]; tot_p = e['h2_%d' % Q]['pred_u2' + key]
            # within part in units of the decade mean: sum_s u2_s m_s / m ; m_s proportional to the pair counts
            ms = np.array([s['h2_%d' % Q]['pairs'] for s in subs], float); m = ms.sum()
            wm = float(np.sum([s['h2_%d' % Q]['measured_u2' + key] for s in subs] * ms) / m)
            wp = float(np.sum([s['h2_%d' % Q]['pred_u2' + key] for s in subs] * ms) / m)
            out['cross'].setdefault(dec, {})['%s%d' % (key, Q)] = {'total': [tot_m, tot_p], 'within': [wm, wp], 'cross': [tot_m - wm, tot_p - wp]}
print('\ncross part of each decade (measured / HL), nominal noise of the cross part ~ 1.33 v / sqrt(K-1)')
for dec, d in out['cross'].items():
    print(dec, '  '.join('%s %.3f/%.3f' % (k, v['cross'][0], v['cross'][1]) for k, v in d.items()))
print('\nscatter of the 27 sub-window ratios (sd / nominal sd) and their mean')
for key in ('X', 'R'):
    for i, L in enumerate(('level1', 'level2', 'level3')):
        r = np.array([out['levels'][w][key][i] for w in out['levels'] if '/sub' in w])
        out['scatter']['%s_%s' % (key, L)] = {'mean': float(r.mean()), 'se_nominal': SE[i] / math.sqrt(len(r)), 'sd': float(r.std(ddof=1)),
                                             'sd_over_nominal': float(r.std(ddof=1) / SE[i])}
        print(key, L, 'mean %.3f +- %.3f  sd %.3f  sd/nominal %.2f' % (r.mean(), SE[i] / math.sqrt(len(r)), r.std(ddof=1), r.std(ddof=1) / SE[i]))
    b = np.array([D[w]['h2_30030']['beta_hat'] for w in D if '/sub' in w]); sb = np.array([D[w]['h2_30030']['se_beta_approx'] for w in D if '/sub' in w])
    if key == 'X':
        print('beta (mod 30030) over 27 sub-windows: mean %.4f, sd %.4f, mean nominal se %.4f' % (b.mean(), b.std(ddof=1), sb.mean()))
        out['scatter']['beta30030'] = {'mean': float(b.mean()), 'sd': float(b.std(ddof=1)), 'se_nominal': float(sb.mean())}
json.dump(out, open('levels.json', 'w'), indent=1)
