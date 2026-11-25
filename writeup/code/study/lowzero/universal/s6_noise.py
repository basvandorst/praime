"""measured noise factor kappa = <(G - pred)^2>/M per family (sharp, real and imaginary parts separately), and chi^2 of the slopes
against the frozen heuristic se (sqrt(0.25 M)/|pred|). Reads results.json, predictions.json."""
import json, math
import numpy as np
P = json.load(open('predictions.json')); R = json.load(open('results.json'))
fam = {}
for k, r in R['per_L'].items():
    e = P['Lfunctions'][k]
    G = complex(*r['G']); pr = complex(*e['pred_sharp']); M = e['M_expected']
    d = G - pr
    fam.setdefault(e['family'], []).append((k, d.real ** 2 / M, d.imag ** 2 / M, r.get('slope'), r.get('se_heur')))
print('family: kappa_re = mean (Re(G-pred))^2/M, kappa_im = mean (Im G-pred)^2/M  [frozen model: 0.25 each]')
allre = []; allim = []; chi2 = []; n = 0
for f, v in fam.items():
    kr = np.mean([x[1] for x in v]); ki = np.mean([x[2] for x in v])
    print('  %-32s n=%d  kappa_re %.3f  kappa_im %.3f' % (f, len(v), kr, ki))
    if not f.startswith('product'):
        allre += [x[1] for x in v]; allim += [x[2] for x in v]
        chi2 += [((x[3] - 1) / x[4]) ** 2 for x in v]
print('  all primitive: kappa_re %.3f  kappa_im %.3f ; chi^2 of slopes vs heuristic se: %.1f on %d L-functions' % (np.mean(allre), np.mean(allim), sum(chi2), len(chi2)))
