# Discovery statistics on the ranges chosen for the registration: y in (0.4, 3] (limit formula) and y in (0, 0.2] (finite-X factor)
import os, sys, json, math, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from finiteX import epsp, jack, coarse, cTT, ym13
pc = json.load(open(os.path.join(HERE, 'pred_cache.json'))); edges = np.array(pc['edges'])
HI = slice(8, 60); LO = slice(0, 4)
for sgn in (-1, 1):
    k = '%+d' % sgn; z = np.load(os.path.join(HERE, 'cubic_disc_%s.npz' % k)); num, den = z['num'], z['den']
    P = np.array(pc[k]['pred']); e = epsp(sgn, 5000.0)
    m, cov, sig, s, se = jack(num[:, HI], den[:, HI], P[HI])
    mc, covc, _, _, _ = jack(coarse(num[:, HI]), coarse(den[:, HI]), coarse(P[HI]) / 4)
    print('D%s0  y in (0.4,3]: s = %.3f +- %.3f; chi2 coarse(13) %.1f; chi2 fine diag(52) %.1f' % ('<' if sgn < 0 else '>', s, se,
          (mc - coarse(P[HI]) / 4) @ np.linalg.solve(covc, mc - coarse(P[HI]) / 4), np.sum((m - P[HI])**2 / sig**2)))
    U = np.array(pc[k]['single_line_unit']); w = 1 / sig**2; a1 = np.sum(U[HI] * m * w) / np.sum(U[HI]**2 * w)
    for name, C in (('zero', 0 * P), ('single_line_zeta3', np.array(pc[k]['single_line_zeta3'])), ('single_line_fit(%.3f)' % a1, a1 * U),
                    ('gl2_prime_level', np.array(pc[k]['gl2_prime_level'])), ('tt/(1+eps)', np.array(pc[k]['tt_secondary']) / (1 + e)),
                    ('lines_no_square_term', np.array(pc[k]['lines_no_square_term'])), ('M1', P + (1 / (1 + e) - 1) * cTT[sgn] * ym13)):
        r = mc - coarse(C[HI]) / 4
        print('     %-26s chi2 coarse %8.1f  fine diag %8.1f' % (name, r @ np.linalg.solve(covc, r), np.sum((m - C[HI])**2 / sig**2)))
    m, cov, sig, s, se = jack(num[:, LO], den[:, LO], P[LO])
    print('     y in (0,0.2]: amplitude of the limit formula %.3f +- %.3f; finite-X factor 1/(1+eps\') = %.3f' % (s, se, 1 / (1 + e)))
