"""dev_closure2.py: development characters (32): per-character statistics and surgeries (chartools.py), the indicator
regression (intercept, ln q, I(p|D) for p = 2, 3, 5, 7), removal and addition means per prime, and the theorem
control for the annealed mean square: the removal ratio of the mean square must be
  1 - 1/p + 2 ln p / (p L'),  L' = ln(q t/2 pi) + 2 gamma + 2 sum_{l | q} ln l/(l - 1)       (twisted second moment)
and the addition ratio 1/(1 - 1/p) exactly (chi(p) = 0 kills every twist by p).  Writes dev_chars.json."""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import obslib, chartools as ct

HERE = os.path.dirname(os.path.abspath(__file__))
DEV8 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'r149', 'dev_values')
DEVX = os.path.join(HERE, 'dev_extra_values')
EG = 0.5772156649015329


def load(D):
    for d in (DEV8, DEVX):
        f = os.path.join(d, 'absL_%d.npy' % D)
        if os.path.exists(f):
            return np.load(f)


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


def twisted_ratio(D, p):
    """theorem value of mean(|L|^2 |E_p|^2)/mean(|L|^2) over [100, 1000] (main terms)."""
    q = abs(D)
    t = obslib.grid()
    Lp = np.log(q * t / (2 * np.pi)) + 2 * EG + 2 * sum(np.log(l) / (l - 1) for l in primes_of(q))
    num = np.trapezoid((1 + 1 / p) * Lp - (2 / p) * (Lp - np.log(p)), t)
    return float(num / np.trapezoid(Lp, t))


def main():
    dev8 = [1001, -1003, 1005, -1048, 1068, -1195, 1253, -1267]
    Ds = dev8 + json.load(open(os.path.join(HERE, 'dev_extra.json')))['dev_extra']
    rows = []
    for D in Ds:
        v = load(D)
        s = ct.stats(v)
        sg = ct.surgery(D, v, s)
        rows.append({'D': D, 'stats': s, 'surg': {str(p): sg[p] for p in ct.PR},
                     'bad': [p for p in ct.PR if D % p == 0], 'lnq': float(np.log(abs(D))),
                     'tw': {str(p): twisted_ratio(D, p) for p in ct.PR if D % p}})
    json.dump(rows, open(os.path.join(HERE, 'dev_chars.json'), 'w'))
    n = len(rows)
    X = np.column_stack([np.ones(n), [r['lnq'] - np.mean([r['lnq'] for r in rows]) for r in rows]] +
                        [[1.0 if p in r['bad'] else 0.0 for r in rows] for p in ct.PR])
    print('n = %d, groups %s' % (n, {p: sum(1 for r in rows if p in r['bad']) for p in ct.PR}))
    print('%-8s | %-36s | %-36s | %-36s | resid sd' % ('stat', 'regression c_p (p missing - present)',
                                                       'removal mean', 'addition mean'))
    out = {}
    for k in ct.KEYS:
        y = np.array([r['stats'][k] for r in rows])
        b, *_ = np.linalg.lstsq(X, y, rcond=None)
        res = y - X @ b; s2 = res @ res / (n - X.shape[1])
        se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
        rem = {p: [r['surg'][str(p)][k] for r in rows if r['surg'][str(p)]['kind'] == 'rem'] for p in ct.PR}
        add = {p: [r['surg'][str(p)][k] for r in rows if r['surg'][str(p)]['kind'] == 'add'] for p in ct.PR}
        out[k] = {'coef': b[2:].tolist(), 'se': se[2:].tolist(), 'resid_sd': float(np.sqrt(s2)),
                  'rem': [float(np.mean(rem[p])) for p in ct.PR],
                  'rem_se': [float(np.std(rem[p], ddof=1) / np.sqrt(len(rem[p]))) for p in ct.PR],
                  'add': [float(np.mean(add[p])) for p in ct.PR],
                  'add_se': [float(np.std(add[p], ddof=1) / np.sqrt(len(add[p]))) for p in ct.PR]}
        o = out[k]
        print('%-8s | %s | %s | %s | %.4f' % (
            k, ' '.join('%+.3f(%.3f)' % (o['coef'][i], o['se'][i]) for i in range(4)),
            ' '.join('%+.3f' % x for x in o['rem']), ' '.join('%+.3f' % x for x in o['add']), o['resid_sd']))
    print('\nse of removal / addition means:')
    for k in ct.KEYS:
        print('%-8s rem %s   add %s' % (k, ' '.join('%.4f' % x for x in out[k]['rem_se']),
                                       ' '.join('%.4f' % x for x in out[k]['add_se'])))
    print('\ntheorem control (annealed mean square), per prime: removal ratio measured / twisted-moment value; '
          'addition ratio x (1 - 1/p)')
    for p in ct.PR:
        rr = [np.exp(r['surg'][str(p)]['lnEZ_1']) / r['tw'][str(p)] for r in rows if r['surg'][str(p)]['kind'] == 'rem']
        aa = [np.exp(-r['surg'][str(p)]['lnEZ_1']) * (1 - 1 / p) for r in rows if r['surg'][str(p)]['kind'] == 'add']
        print('p=%d removal %.4f +- %.4f (n=%d)   addition %.4f +- %.4f (n=%d)   [leading-only removal value '
              '(1-1/p)/twisted = %.4f]' % (p, np.mean(rr), np.std(rr) / np.sqrt(len(rr)), len(rr), np.mean(aa),
                                          np.std(aa) / np.sqrt(len(aa)), len(aa),
                                          np.mean([(1 - 1 / p) / r['tw'][str(p)] for r in rows if str(p) in r['tw']])))
    json.dump(out, open(os.path.join(HERE, 'dev_summary.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
