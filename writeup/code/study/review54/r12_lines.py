# Lines r = 2, 1/2, 3, 1/3: own sums, block-jackknife errors (K = 4..30 averaged), pooled slopes.
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+"/")
import json, math, numpy as np
from r1_sums import theta, own_data, D
pred = json.load(open(D + 'predictions.json'))['Lfunctions']
L = own_data()
for key in ('X5', 'X5d', 'X7', 'X7d'):
    L[key] = dict(Q=float(pred[key]['Q']), mu=pred[key]['mu'], eps=complex(*pred[key]['eps']), fam='GL1')
acc = {}
for key, v in L.items():
    P = pred[key]; T2 = float(P['T2'])
    g = np.load(D + 'zeros/zeros_%s.npz' % key)['g']; g = np.sort(g[(g > 0) & (g <= T2)])
    e = np.exp(2j * theta(g, v['Q'], v['mu']))
    for rs, lv in P['lines'].items():
        r = eval(rs.replace('/', '/1.0/')) if '/' in rs else float(rs)
        pr = complex(*lv["sharp"])
        if abs(pr) < 1e-6: print("   zero line", key, rs); continue
        w = e * r ** (-1j * g)
        G = w.sum()
        ses = []
        for K in (4, 6, 8, 12, 20, 30):
            edges = np.linspace(0, T2, K + 1)
            Gk = np.array([w[(g > edges[j]) & (g <= edges[j + 1])].sum() for j in range(K)])
            ses.append(((Gk * np.conj(pr) / abs(pr)).real.std(ddof=1) * math.sqrt(K)))
        var = np.mean(np.square(ses))
        a = acc.setdefault(rs, [0.0, 0.0, 0.0])
        a[0] += (G * np.conj(pr)).real / var; a[1] += abs(pr) ** 2 / var
        b = acc.setdefault('all', [0.0, 0.0, 0.0]); b[0] += (G * np.conj(pr)).real / var; b[1] += abs(pr) ** 2 / var
        fam = acc.setdefault('fam ' + v['fam'], [0.0, 0.0, 0.0]); fam[0] += (G * np.conj(pr)).real / var; fam[1] += abs(pr) ** 2 / var
for k, (num, den, _) in acc.items():
    print('%-28s inverse-variance slope %.3f +- %.3f' % (k, num / den, 1 / math.sqrt(den)))
