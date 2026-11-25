"""Noise model and power.  Noise model from the 190 development pairs (run_dev_data, used data): residuals of Re X after
the law decomposed into per-curve additive (u), per-curve multiplicative (v, scaling the pair's mean) and pair (e)
parts, for Re X, for Re X - measured r=1 product (X5) and for Im X; checked by reproducing the development jackknife se.
Then simulated on the fresh design (frozen predictions only, no fresh zero) under each hypothesis: expected value and
se of every registered pair reading.  One-point line readings: se scaled from development by sqrt(20/138) and checked
by a per-curve simulation.  usage: s7_power.py [NREP]"""
import sys, os, math, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ecross_test as E
nrep = int(sys.argv[1]) if len(sys.argv) > 1 else 40
rng = np.random.default_rng(20261125)
Fd = E.load_frozen(HERE, 'dev'); Ff = E.load_frozen(HERE, 'fresh')
Rd = os.path.join(HERE, 'run_dev_data'); O = np.load(os.path.join(Rd, 'onepoint.npz'))
out = open(os.path.join(HERE, 's7_power_out.txt'), 'w')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')

def decompose(I, J, n, r, p):
    A = np.zeros((len(r), n)); A[np.arange(len(r)), I] = 1; A[np.arange(len(r)), J] = 1
    B = A * p[:, None]; v = np.linalg.lstsq(np.hstack([A, B]), r, rcond=None)[0]; e = r - np.hstack([A, B]) @ v
    dof = len(r) - 2 * n
    return dict(sd=float(r.std()), sd_u=float(v[:n].std()), sd_v=float(v[n:].std()), sd_e=float(np.sqrt(np.sum(e ** 2) / dof)))

def noise(m, I, J, n, mean):
    u = rng.normal(0, m['sd_u'], n); v = rng.normal(0, m['sd_v'], n)
    return u[I] + u[J] + (v[I] + v[J]) * mean + rng.normal(0, m['sd_e'], len(I))

W1, W2 = E.WIN; T = W2 - W1
NM = {}
for s in E.SIGMAS:
    X = np.load(os.path.join(Rd, 'X_%.1f.npy' % s)); I, J = Fd['I'], Fd['J']; n = len(Fd['ids'])
    x = X[I, J]; p = Fd['pred_%.1f' % s]; p1 = Fd['pred1_%.1f' % s]
    P0 = s * math.sqrt(2 * math.pi) * np.real(O['Gw'][I, 0] * np.conj(O['Gw'][J, 0])) / T ** 2
    m_re = decompose(I, J, n, x.real - p, p); m_5 = decompose(I, J, n, x.real - P0 - (p - p1), p - p1); m_im = decompose(I, J, n, x.imag - x.imag.mean(), p)
    NM[s] = dict(re=m_re, x5=m_5, im=m_im)
    P('sigma %.1f development noise: Re X - law: sd %.4f = u %.4f, v %.3f, e %.4f | X5 residual: sd %.4f (u %.4f, v %.3f, e %.4f) | Im X: sd %.4f' % (
        s, m_re['sd'], m_re['sd_u'], m_re['sd_v'], m_re['sd_e'], m_5['sd'], m_5['sd_u'], m_5['sd_v'], m_5['sd_e'], m_im['sd']))
    # check on the development design: simulated sd of the X1 slope against the data jackknife
    fit = E.Fit(I, J, n); sl = []; sl5 = []
    for rep in range(200):
        y = p + noise(m_re, I, J, n, p); sl.append(np.linalg.lstsq(np.column_stack([np.ones(len(p)), p]), y, rcond=None)[0][1])
        y5 = (p - p1) + noise(m_5, I, J, n, p - p1); sl5.append(np.linalg.lstsq(np.column_stack([np.ones(len(p)), p - p1]), y5, rcond=None)[0][1])
    jk = fit.run(np.column_stack([np.ones(len(p)), p]), x.real)['se'][1]; jk5 = fit.run(np.column_stack([np.ones(len(p)), p - p1]), x.real - P0)['se'][1]
    P('   check on development: X1 slope sd simulated %.4f vs data jackknife %.4f; X5 %.4f vs %.4f' % (np.std(sl), jk, np.std(sl5), jk5))
json.dump({str(k): v for k, v in NM.items()}, open(os.path.join(HERE, 'noise_model.json'), 'w'), indent=1)

# fresh design
I, J = Ff['I'], Ff['J']; n = len(Ff['ids']); fit = E.Fit(I, J, n); one = np.ones(len(I))
res = {}
for s in E.SIGMAS:
    p = Ff['pred_%.1f' % s]; p1 = Ff['pred1_%.1f' % s]; pn = Ff['prednoroot_%.1f' % s]; pr = Ff['predrev_%.1f' % s]
    P('fresh sigma %.1f: %d pairs; law mean %.4f sd %.4f; H=1 part sd %.4f; H-1 part sd %.4f; corr(law, H=1) %.3f; corr(H-1 law, H-1 reversed) %.3f; corr(law, no-root) %.3f' % (
        s, len(p), p.mean(), p.std(), p1.std(), (p - p1).std(), np.corrcoef(p, p1)[0, 1], np.corrcoef(p - p1, pr - p1)[0, 1], np.corrcoef(p, pn)[0, 1]))
    hyps = {'H1 law': p, 'H0 none': 0 * p, 'H1 wrong sign': -p, 'H1 no root numbers': pn, 'H=1 only': p1, 'reversed pairing': pr, 'half': 0.5 * p}
    designs = {'X1 [icpt, law]': np.column_stack([one, p]), 'X2 [icpt, H=1, H-1]': np.column_stack([one, p1, p - p1]),
               'X3 [icpt, H=1, H-1, H-1 rev]': np.column_stack([one, p1, p - p1, pr - p1]), 'X4 [icpt, law, no-root]': np.column_stack([one, p, pn]),
               'X8 mean': one[:, None]}
    for hn, mu in hyps.items():
        for dn, A in designs.items():
            ex = np.linalg.lstsq(A, mu, rcond=None)[0]
            reps = []
            for rep in range(nrep if hn in ('H1 law', 'H0 none') else 0):
                y = mu + noise(NM[s]['re'], I, J, n, mu); reps.append(np.linalg.lstsq(A, y, rcond=None)[0])
            se = np.std(reps, axis=0) if reps else None
            res[(s, hn, dn)] = (ex, se)
            P('   %-20s %-30s expected %s%s' % (hn, dn, ', '.join('%.4f' % v for v in ex), ('   se ' + ', '.join('%.4f' % v for v in se)) if se is not None else ''))
    # X5 (measured r = 1 product removed) on the H-1 part: under H1 coefficient 1; reversed: projection; H=1 only: 0
    A5 = np.column_stack([one, p - p1])
    for hn, mu in (('H1 law', p - p1), ('H=1 only', 0 * p), ('reversed pairing', pr - p1), ('H1 no root numbers', pn - p1 * 0), ('H0 none', 0 * p)):
        ex = np.linalg.lstsq(A5, mu, rcond=None)[0]
        reps = [np.linalg.lstsq(A5, mu + noise(NM[s]['x5'], I, J, n, mu), rcond=None)[0] for _ in range(nrep)] if hn == 'H1 law' else []
        P('   %-20s %-30s expected %s%s' % (hn, 'X5 [icpt, H-1] after r=1 removal', ', '.join('%.4f' % v for v in ex), ('   se ' + ', '.join('%.4f' % v for v in np.std(reps, axis=0))) if reps else ''))
    # jackknife se on one simulated H1 data set (the test's own error estimate)
    y = p + noise(NM[s]['re'], I, J, n, p)
    P('   one simulated H1 data set: X1 slope %.4f +- %.4f (jackknife over curves)' % tuple(np.array([fit.run(designs['X1 [icpt, law]'], y)[k][1] for k in ('coef', 'se')])))
# one-point lines on fresh: per-curve simulation from the development residual rms of each line
Gl = O['Gl']; law = Fd['lineamp'] * (E.BAND[1] - E.BAND[0]) / (2 * math.pi)
rms = np.sqrt(np.mean((Gl[:, 1:, 0].real - law[:, 1:]) ** 2, axis=0))
lawf = Ff['lineamp'][:, 1:] * (E.BAND[1] - E.BAND[0]) / (2 * math.pi)
sl = []
for rep in range(400):
    y = lawf + rng.normal(0, 1, lawf.shape) * rms[None, :]; sl.append(np.sum(y * lawf) / np.sum(lawf ** 2))
P('one-point lines r != 1 (fresh, 138 x 12): pooled slope se (simulated, development residual rms per line %s): %.4f; law without root numbers: expected slope %.4f' % (
    np.round(rms, 1).tolist(), np.std(sl), np.sum(Ff['lineamp_noroot'][:, 1:] * Ff['lineamp'][:, 1:]) / np.sum(Ff['lineamp'][:, 1:] ** 2)))
out.close()
