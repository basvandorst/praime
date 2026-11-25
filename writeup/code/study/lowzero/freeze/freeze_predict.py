"""freeze_predict.py: frozen predictions for the test set (reads only the discriminant list test_set.json and the
development summary dev_summary.json, which is reference material, not a hypothesis).

Per character i (q = |D|) and prime p in {2, 3, 5, 7}:
  Lbar_i   = mean over the 143 intervals of ln(q t_j / 2 pi), t_j the interval centres (grid of freeze_test.py);
  L'_ip    = Lbar_i + 2 gamma + 2 sum_{l | q, l != p} ln l / (l - 1);
  lnR_ip   = ln(1 - 1/p) + ln(1 + 2 ln p / ((p - 1) L'_ip))   (per-prime ratio of the mean square, theorem:
             (1/T) int |L(1/2+it, chi)|^2 dt = (phi(q)/q)(ln(qt/2 pi) + 2 gamma + 2 sum_{l|q} ln l/(l-1)) + o(1)).
Hypotheses, value of 'p missing minus p present' (Delta):
  crit   frozen with the critical (second) moment's exact per-prime factor:  M: lnR/2;  lnZ_b (b > 1): b lnR;
         slope32 = lnZ_3 - lnZ_2: lnR;  lnEZ_1: lnR (theorem).
  lead   frozen with the leading Euler factor only (independent probe, GHK):  ln(1 - 1/p) in place of lnR.
  naive  not frozen, quenched = annealed at each b:  lnZ_b: -ln 2F1(b, b; 1; 1/p);  slope32 accordingly;  M: none.
  spread not frozen, Gibbs measure uniform on the interval: lnZ_b: -E_phi ln (1/theta) int_0^theta e^{2b Y_p};
         M: -E_phi max_h Y_p(h).
  zero   no dependence on the primes dividing the conductor.
Character-level predictions y_i^H(s) = sum_{p | D_i, p in 2,3,5,7} Delta^H_p(i); regression-matched coefficients =
the same OLS as freeze_test.py (columns 1, ln q - mean, I(p|D) for p = 2, 3, 5, 7, 11, 13) applied to y^H.
usage: freeze_predict.py test_set.json dev_summary.json out.json"""
import sys, json, hashlib
import numpy as np
from scipy.special import hyp2f1

EG = 0.5772156649015329
THETA = 2 * np.pi
PR = (2, 3, 5, 7)
REGP = (2, 3, 5, 7, 11, 13)
BETAS = (1.5, 2.0, 3.0)


def primes_of(n):
    n = abs(n); out = []; d = 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def spread_values():
    phis = (np.arange(4000) + 0.5) * 2 * np.pi / 4000
    h = (np.arange(4000) + 0.5) * THETA / 4000
    out = {}
    for p in PR:
        Y = -np.log(np.abs(1 - p ** -0.5 * np.exp(1j * (phis[:, None] + h[None, :] * np.log(p)))))
        d = {'M': float(-np.mean(Y.max(axis=1)))}
        for b in BETAS:
            d['lnZ_%g' % b] = float(-np.mean(np.log(np.mean(np.exp(2 * b * Y), axis=1))))
        d['slope32'] = d['lnZ_3'] - d['lnZ_2']
        out[p] = d
    return out


def design(Ds):
    lnq = np.log(np.abs(np.array(Ds, float)))
    return np.column_stack([np.ones(len(Ds)), lnq - lnq.mean()] +
                           [[1.0 if D % p == 0 else 0.0 for D in Ds] for p in REGP])


def main():
    Ds = json.load(open(sys.argv[1]))['test']
    dev = json.load(open(sys.argv[2]))
    tcen = 100 + THETA * (np.arange(143) + 0.5)
    sp = spread_values()
    stats = ['M', 'lnZ_1.5', 'lnZ_2', 'lnZ_3', 'slope32', 'lnEZ_1']
    hyps = ['crit', 'lead', 'naive', 'spread', 'zero']
    chars = []
    yH = {h: {s: np.zeros(len(Ds)) for s in stats} for h in hyps}
    for i, D in enumerate(Ds):
        q = abs(D)
        Lbar = float(np.mean(np.log(q * tcen / (2 * np.pi))))
        qp = primes_of(q)
        c = {'D': D, 'q': q, 'Lbar': Lbar, 'bad': [p for p in REGP if D % p == 0], 'lnR': {}, 'tw_note': 'see freeze_test'}
        for p in PR:
            Lp = Lbar + 2 * EG + 2 * sum(np.log(l) / (l - 1) for l in qp if l != p)
            lnR = float(np.log(1 - 1 / p) + np.log(1 + 2 * np.log(p) / ((p - 1) * Lp)))
            c['lnR'][str(p)] = lnR
            if D % p:
                continue
            lead = float(np.log(1 - 1 / p))
            vals = {
                'crit': {'M': lnR / 2, 'lnZ_1.5': 1.5 * lnR, 'lnZ_2': 2 * lnR, 'lnZ_3': 3 * lnR, 'slope32': lnR,
                         'lnEZ_1': lnR},
                'lead': {'M': lead / 2, 'lnZ_1.5': 1.5 * lead, 'lnZ_2': 2 * lead, 'lnZ_3': 3 * lead, 'slope32': lead,
                         'lnEZ_1': lead},
                'naive': {'M': float('nan'), 'lnZ_1.5': -np.log(hyp2f1(1.5, 1.5, 1, 1 / p)),
                          'lnZ_2': -np.log(hyp2f1(2, 2, 1, 1 / p)), 'lnZ_3': -np.log(hyp2f1(3, 3, 1, 1 / p)),
                          'slope32': -np.log(hyp2f1(3, 3, 1, 1 / p)) + np.log(hyp2f1(2, 2, 1, 1 / p)),
                          'lnEZ_1': lead},
                'spread': dict(sp[p], lnEZ_1=lead),
                'zero': {s: 0.0 for s in stats}}
            for h in hyps:
                for s in stats:
                    yH[h][s][i] += vals[h][s]
        chars.append(c)
    X = design(Ds)
    coef = {}
    for h in hyps:
        coef[h] = {}
        for s in stats:
            y = yH[h][s]
            if np.any(np.isnan(y)):
                coef[h][s] = None
                continue
            b, *_ = np.linalg.lstsq(X, y, rcond=None)
            coef[h][s] = {str(p): float(b[2 + k]) for k, p in enumerate(REGP)}
    out = {'what': 'frozen predictions for the per-prime readings (p missing minus p present); see freeze_predict.py',
           'test_set': Ds, 'n': len(Ds), 'regression_columns': ['1', 'lnq-mean'] + ['I(%d|D)' % p for p in REGP],
           'chars': chars,
           'per_char': {h: {s: [float(x) for x in yH[h][s]] for s in stats} for h in hyps},
           'coef': coef,
           'closed_forms': {'spread': {str(p): sp[p] for p in PR},
                            'naive': {str(p): {'lnZ_%g' % b: float(-np.log(hyp2f1(b, b, 1, 1 / p))) for b in BETAS}
                                      for p in PR},
                            'lead': {str(p): float(np.log(1 - 1 / p)) for p in PR}},
           'relational': {'H_W': 'c_p = r_p (removal mean on the test set) for M exactly; for lnZ_b up to the warp '
                                 'Jacobian (|.| < 0.01, theory uncertainty)',
                          'H_G': 'c_p = a_p (addition mean on the test set)'},
           'dev_reference': {k: {'rem': dev[k]['rem'], 'add': dev[k]['add'], 'coef': dev[k]['coef'],
                                 'se': dev[k]['se'], 'resid_sd': dev[k]['resid_sd']} for k in dev}}
    s = json.dumps(out, indent=1, sort_keys=True)
    open(sys.argv[3], 'w').write(s)
    print('sha256', hashlib.sha256(s.encode()).hexdigest())
    print('regression-matched predicted coefficients (p = 2, 3, 5, 7):')
    for h in hyps:
        for st in stats:
            if coef[h][st] is None:
                continue
            print('  %-6s %-8s %s' % (h, st, ' '.join('%+.4f' % coef[h][st][str(p)] for p in PR)))
    print('mean Lbar %.3f (range %.3f..%.3f)' % (np.mean([c['Lbar'] for c in chars]), min(c['Lbar'] for c in chars),
                                                 max(c['Lbar'] for c in chars)))


if __name__ == '__main__':
    main()
