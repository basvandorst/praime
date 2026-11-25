"""Regression-matched predicted coefficients from predictions_lprime.json alone (no data): same design as lprime_test.py."""
import json, numpy as np, importlib.util, sys
spec = importlib.util.spec_from_file_location('lt', sys.argv[2]); lt = importlib.util.module_from_spec(spec); spec.loader.exec_module(lt)
J = json.load(open(sys.argv[1]))
for sets in (['R120'], ['Dn'], ['Fn'], ['R120', 'Dn', 'Fn'], ['B']):
    rows = [dict(e, set=s) for s in sets for e in J[s]]; X, ns = lt.design(rows); ip = ns + 1
    H = np.linalg.pinv(X.T @ X) @ X.T; pr = lambda k: np.array([r[k] for r in rows]); m = lambda v: (H @ v)[ip:ip + 6]
    print('== %s (%d L-functions; mean 1/L %.4f)' % ('+'.join(sets), len(rows), np.mean(pr('mean_invL'))))
    for key in ('cum1', 'cum2', 'cum3', 'cum4', 'logmom-0.5', 'logmom+0.5', 'logmom-1.0', 'logmom+1.0'):
        for h in ('ldm', 'lead'):
            print('  %-11s %-4s ' % (key, h) + ' '.join('%+.4f' % v for v in m(pr(key + '_' + h))))
    P1 = pr('S1_recipe')
    print('  S1 (log, vs exact recipe) none  ' + ' '.join('%+.4f' % v for v in m(np.log(pr('S1_none') / P1))) + '   lead ' + ' '.join('%+.4f' % v for v in m(np.log(pr('S1_lead') / P1))))
    Pc = pr('Sm1_conj') - 2 * lt.YC * pr('tail_C_weighted'); Pn = pr('Sm1_none') - 2 * lt.YC * pr('tail_C_weighted') * np.exp(-pr('tail_logfactor_lead'))
    print('  Sm1 (log, vs conjecture) none   ' + ' '.join('%+.4f' % v for v in m(np.log(Pn / Pc))))
    print('  groups', [int(sum(p in r['bad'] for r in rows)) for p in lt.PR])
