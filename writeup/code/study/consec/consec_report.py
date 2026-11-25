"""consec_report.py: compact printout of a readings file of consec_analyze.py (and the LOS N2 numbers for its range).
usage: consec_report.py READINGS.json"""
import sys, json, math
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import numpy as np
from consec_null import tables
R = json.load(open(sys.argv[1])); x0, x1 = R['range']; nul = tables(x0, x1)
print('range [%d, %d)  primes %d  segments %d' % (x0, x1, R['nprimes'], R['segments']))
allz = []
for nm in [k for k in R if isinstance(R[k], dict) and 'corr' in R[k]]:
    E = R[nm]; print('\n=== %s ===' % nm)
    if 'z' in E['marg']: print('  marginal z vs Chebotarev:', ' '.join('%+.2f' % z for z in E['marg']['z']))
    if 'ab_b' in E:
        print('  abelian b_k (measured | LOS N2):', '  '.join('k=%s %+.5f | %+.5f' % (k, v['b'], nul[nm]['b_LOS'][k] if nm in nul and 'b_LOS' in nul[nm] else float('nan')) for k, v in E['ab_b'].items()),
              '  sd ~ %.5f' % (1 / math.sqrt(E['ab_b']['1']['n'])))
    zs = [(k, v['z'], v['r'], v['n']) for k, v in E['corr'].items()]
    print('  fibre interaction (consecutive): ' + '  '.join('%s r=%+.5f z=%+.2f' % (k, r, z) for k, z, r, n in zs))
    zsM = [(k, v['z'], v['r']) for k, v in E['corrM'].items()]
    print('  fibre interaction (next fibre prime): ' + '  '.join('%s z=%+.2f' % (k, z) for k, z, r in zsM))
    allz += [z for k, z, r, n in zs]
    if 'indep' in E: print('  contingency inside abelian fibres:', '  '.join('k=%s chi2=%.1f/%d p=%.3f' % (k, v['chi2'], v['dof'], v['p']) for k, v in E['indep'].items()))
    nb = [(k, v['z']) for k, v in E['nbr'].items()]
    if nb: print('  f vs neighbour abelian datum: ' + '  '.join('%s %+.2f' % kv for kv in nb))
    for f, g in E['gap'].items():
        s = '  gap-conditioned %s: chi2=%.1f/%d p=%.3f max|z|=%.2f' % (f, g['chi2'], g['dof'], g['p'], g['max_abs_z'])
        if 'gap_mult' in g: s += '   gaps=0 mod q: r=%+.4f z=%+.2f (n=%d); rest z=%+.2f' % (g['gap_mult']['r'], g['gap_mult']['z'], g['gap_mult']['n'], g['gap_rest']['z'])
        print(s)
    if 'triple' in E: print('  triple (1st-3rd | middle): chi2=%.1f/%d p=%.3f' % (E['triple']['chi2'], E['triple']['dof'], E['triple']['p']))
    if 'runs' in E: print('  runs vs geometric: ' + '  '.join('%s: chi2=%.1f/%d mean %.4f (geo %.4f)' % (c, v['chi2'], v['dof'], v['mean_len'], v['geo_mean']) for c, v in E['runs'].items() if c != 'total') + '  total p=%.3f' % E['runs']['total']['p'])
    w = E['w4']; print('  per-block 1e4 fibre table: within %.4f +- %.4f   column %.4f +- %.4f  (%d blocks)' % (w['within'], w['within_err'], w['col'], w['col_err'], w['blocks']))
    if 'w4a' in E: w = E['w4a']; print('  per-block 1e4 abelian sign table: within %.4f +- %.4f   column %.4f +- %.4f' % (w['within'], w['within_err'], w['col'], w['col_err']))
    for key, ent in E['disp'].items():
        H, b0, b1 = key.split('_')
        if int(b0) != x0 or int(b1) != x1: continue
        print('  dispersion %-8s ' % H + '  '.join('%s R=%.4f+-%.4f' % (f, v['R'], v['err']) for f, v in ent.items()))
print('\nall consecutive fibre-interaction z: n=%d  sum z^2 = %.1f  max|z| = %.2f' % (len(allz), sum(z * z for z in allz), max(abs(z) for z in allz)))
