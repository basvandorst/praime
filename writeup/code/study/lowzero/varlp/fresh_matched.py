"""fresh_matched.py: regression-matched frozen coefficients for the fresh sets C + Bp (all L-functions; the completeness
rule may drop a few, var_test.py re-matches on the kept ones), the theory's own numerical uncertainty matched the same way,
expected measurement se (round 142's bootstrap se scaled by sqrt(zeros_142 / zeros_fresh)), and the power of each reading.
usage: fresh_matched.py predictions_varlp.json
"""
import os, sys, json, math
sys.dont_write_bytecode = True
import numpy as np
import importlib.util
spec = importlib.util.spec_from_file_location('lprime_test', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lprime') + '/lprime_test.py')
LT = importlib.util.module_from_spec(spec); spec.loader.exec_module(LT)
PR = LT.PR
SE142 = np.array([0.00160, 0.00172, 0.00197, 0.00237, 0.00344, 0.00340]); N142 = 449710
P = json.load(open(sys.argv[1])); rows = [dict(e, set=s) for s in ('C', 'Bp') for e in P[s]]
nz = sum(r['n_zeros_expected'] for r in rows); scale = math.sqrt(N142 / nz); se = SE142 * scale
X, ns = LT.design(rows); ip = ns + 1; H = np.linalg.pinv(X.T @ X) @ X.T
m = lambda k: (H @ np.array([r[k] for r in rows]))[ip:ip + len(PR)]
M = {h: m('cum2_' + h) for h in ('3pt', 'ldm', 'lead', 'none')}
U = np.abs(H[ip:ip + len(PR)]) @ np.array([r['cum2_3pt_numunc'] for r in rows])     # conservative: |H| times per-function uncertainty
print('L-functions %d, expected zeros %.0f (round 142: %d) -> se scale %.3f' % (len(rows), nz, N142, scale))
print('groups (2,3,5,7,11,13): %s' % [sum(p in r['bad'] for r in rows) for p in PR])
print('mean L range: %.2f .. %.2f; zero-weighted <L> %.2f' % (min(math.exp(r['mean_lnL']) for r in rows), max(math.exp(r['mean_lnL']) for r in rows),
                                                          math.exp(sum(r['mean_lnL'] * r['n_zeros_expected'] for r in rows) / nz)))
print('\np:                     ' + ' '.join('%9d' % p for p in PR))
for h in ('3pt', 'ldm', 'lead', 'none'): print('matched %-5s          ' % h + ' '.join('%+9.5f' % x for x in M[h]))
print('numerical unc. of 3pt  ' + ' '.join('%9.5f' % x for x in U))
print('expected se            ' + ' '.join('%9.5f' % x for x in se))
for h in ('ldm', 'lead', 'none'): print('power |3pt - %-4s|/se  ' % h + ' '.join('%9.1f' % x for x in np.abs(M['3pt'] - M[h]) / se))
J = [0, 1, 2]; w = 1 / se[J] ** 2
for h in ('lead', 'ldm'):
    d = M['3pt'][J] - M[h][J]; sg = 1 / math.sqrt((w * d * d).sum())
    print('amplitude (0 = %s, 1 = 3pt) over p = 2, 3, 5: expected se %.4f, separation %.1f se' % (h, sg, 1 / sg))
json.dump({'primes': list(PR), 'matched': {h: list(v) for h, v in M.items()}, 'numunc_3pt': list(U), 'expected_se': list(se),
           'se_scale': scale, 'zeros_expected': nz}, open('fresh_matched.json', 'w'), indent=1)
