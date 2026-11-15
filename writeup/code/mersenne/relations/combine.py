"""Empirical p-values and family-wise correction (Westfall-Young single-step min-p,
computed with leave-one-out p-values of the null samples themselves), per null family,
per angle and over all statistics of all angles."""
import json, os
import numpy as np


def p_real(r, v, side):
    n = v.size
    hi = (1 + np.sum(v >= r)) / (n + 1)
    lo = (1 + np.sum(v <= r)) / (n + 1)
    return {'hi': hi, 'lo': lo, 'two': min(1.0, 2 * min(hi, lo))}[side]


def p_loo(v, side):
    n = v.size
    s = np.sort(v)
    ge = n - np.searchsorted(s, v, 'left') - 1
    le = np.searchsorted(s, v, 'right') - 1
    hi = (1 + ge) / n
    lo = (1 + le) / n
    if side == 'hi':
        return hi
    if side == 'lo':
        return lo
    return np.minimum(1.0, 2 * np.minimum(hi, lo))


def main():
    stats = []
    for ang, fn in ((1, 'angle1.json'), (2, 'angle2.json'), (2, 'angle2b.json'), (3, 'angle3.json'),
                    (4, 'angle4.json'), (5, 'angle5.json')):
        if not os.path.exists(fn):
            continue
        d = json.load(open(fn))
        for name, rec in d.items():
            stats.append((ang, name, rec))
    lines = []
    for fam in ('a', 'b'):
        names = [(a, n) for a, n, _ in stats]
        pr = np.array([p_real(float(rec['real']), np.array(rec[fam], float), rec['side']) for _, _, rec in stats])
        PL = np.array([p_loo(np.array(rec[fam], float), rec['side']) for _, _, rec in stats]).T  # (nnull, nstat)
        nn = PL.shape[0]
        lines.append('=== null (%s): %d null sequences, %d statistics ===' % (fam, nn, len(stats)))
        minp_null = PL.min(1)
        adj_all = (1 + np.sum(minp_null[None, :] <= pr[:, None], axis=1)) / (nn + 1)
        lines.append('global min-p of real: %.4f (%s), family-wise adjusted p = %.3f; Bonferroni = %.3f'
                     % (pr.min(), names[int(pr.argmin())], adj_all[int(pr.argmin())], min(1, pr.min() * len(pr))))
        for ang in range(1, 6):
            idx = [i for i, (a, _) in enumerate(names) if a == ang]
            if not idx:
                continue
            mp = PL[:, idx].min(1)
            i0 = idx[int(np.argmin(pr[idx]))]
            adj_ang = (1 + np.sum(mp <= pr[i0])) / (nn + 1)
            lines.append('angle %d: %3d stats; smallest raw p %.4f (%s, real %s, null median %s); '
                         'adjusted within angle %.3f, over all %.3f'
                         % (ang, len(idx), pr[i0], names[i0][1], fmt(stats[i0][2]['real']),
                            fmt(np.median(np.array(stats[i0][2][fam], float))), adj_ang, adj_all[i0]))
        order = np.argsort(pr)
        lines.append('ten smallest raw p:')
        for i in order[:10]:
            v = np.array(stats[i][2][fam], float)
            lines.append('  A%d %-14s side=%-3s real=%-10s null median=%-10s 5%%-95%%=[%s, %s]  p=%.4f  adj=%.3f'
                         % (names[i][0], names[i][1], stats[i][2]['side'], fmt(stats[i][2]['real']),
                            fmt(np.median(v)), fmt(np.percentile(v, 5)), fmt(np.percentile(v, 95)), pr[i], adj_all[i]))
        lines.append('all raw p (angle, stat, real, median null, p):')
        for i, (a, nme) in enumerate(names):
            v = np.array(stats[i][2][fam], float)
            lines.append('  A%d %-14s real=%-10s med=%-10s p=%.4f' % (a, nme, fmt(stats[i][2]['real']),
                                                                    fmt(np.median(v)), pr[i]))
    open('combined.txt', 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(l for l in lines if not l.startswith('  A') or 'adj=' in l))


def fmt(x):
    x = float(x)
    return ('%d' % x) if x == int(x) and abs(x) < 1e9 else ('%.4g' % x)


if __name__ == '__main__':
    main()
