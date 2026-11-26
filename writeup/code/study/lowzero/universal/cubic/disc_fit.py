# Discovery fit (5000 <= |D| < 10000): binned E[sqrt|D| a_K(p)] with jackknife over fields vs the prediction and comparators
import os, sys, json, math, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
pc = json.load(open(os.path.join(HERE, 'pred_cache.json'))); edges = np.array(pc['edges']); yc = (edges[:-1] + edges[1:]) / 2
def jack(num, den, P=None):
    n = len(num); N = num.sum(0); D = den.sum(0); m = N / D; mj = (N - num) / (D - den)
    dm = mj - mj.mean(0); cov = (n - 1) / n * dm.T @ dm; sig = np.sqrt(np.diag(cov)); out = dict(m=m, cov=cov, sig=sig, mj=mj, n=n)
    if P is not None:
        w = 1 / sig**2; s = np.sum(P * m * w) / np.sum(P * P * w); sj = (mj * (P * w)).sum(1) / np.sum(P * P * w)
        out['s'] = s; out['se'] = math.sqrt((n - 1) / n * np.sum((sj - sj.mean())**2))
    return out
def coarse(x, g=4): return x.reshape(x.shape[:-1] + (x.shape[-1] // g, g)).sum(-1)
for sgn in (-1, 1):
    k = '%+d' % sgn; z = np.load(os.path.join(HERE, 'cubic_disc_%s.npz' % k)); num, den = z['num'], z['den']
    P = np.array(pc[k]['pred']); J = jack(num, den, P)
    print('=== D %s 0: %d fields; mean |D| %.0f' % ('<' if sgn < 0 else '>', J['n'], np.abs(z['D']).mean()))
    print('amplitude s (data = s * prediction, weights 1/sig^2): %.3f +- %.3f' % (J['s'], J['se']))
    Jc = jack(coarse(num), coarse(den)); 
    models = {'prediction': P}
    for name in ('lines_no_square_term', 'single_line_zeta3', 'gl2_prime_level', 'tt_secondary'): models[name] = np.array(pc[k][name])
    models['zero'] = 0 * P
    U = np.array(pc[k]['single_line_unit']); w = 1 / J['sig']**2; a1 = np.sum(U * J['m'] * w) / np.sum(U * U * w); models['single_line_fit'] = a1 * U
    print('fitted single-line amplitude (as multiple of 2 pi K(4 pi sqrt y)): %.3f' % a1)
    for name, M in models.items():
        r = Jc['m'] - coarse(M) / 4; c2 = float(r @ np.linalg.solve(Jc['cov'], r)); d2 = float(np.sum((J['m'] - M)**2 / J['sig']**2))
        print('  %-22s chi2 coarse(15, jack cov) %7.1f   chi2 fine(60, diag) %7.1f' % (name, c2, d2))
    print('  y    data   err   pred   (fine bins)')
    for j in range(60):
        print('  %.3f %6.2f %5.2f %6.2f' % (yc[j], J['m'][j], J['sig'][j], P[j]))
