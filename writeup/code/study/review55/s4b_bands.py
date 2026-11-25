"""Slope by height band (fraction of the window) for the Artin factors and for the primitive families; jackknife of sharp - Hann."""
import json, math, os
import numpy as np
from scipy.special import loggamma
B = os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/universal/blind'; D = os.path.dirname(os.path.abspath(__file__))
P = json.load(open(B + '/predictions_blind.json'))
def theta(t, Q, mu):
    t = np.asarray(t, float); s = t / 2 * math.log(Q / math.pi ** len(mu))
    for m in mu: s = s + np.imag(loggamma((0.5 + m + 1j * t) / 2))
    return s
zz = np.load(D + '/zeta_zeros_310.npy')
edges = np.linspace(0, 1, 7)
def bands(objs):
    # objs: list of (zeros, Q, mu, T2, pred_complex)  -> per band pooled slope with jackknife over objects
    nb = len(edges) - 1; num = np.zeros((len(objs), nb)); den = np.zeros((len(objs), nb)); numh = np.zeros(len(objs)); denh = np.zeros(len(objs))
    nums = np.zeros(len(objs)); dens = np.zeros(len(objs))
    for i, (g, Q, mu, T2, pr) in enumerate(objs):
        g = g[(g > 0) & (g <= T2)]; ph = np.exp(2j * theta(g, Q, mu))
        for b in range(nb):
            m = (g > edges[b] * T2) & (g <= edges[b + 1] * T2)
            p = pr * (edges[b + 1] - edges[b])
            num[i, b] = (complex(ph[m].sum()) * p.conjugate()).real; den[i, b] = abs(p) ** 2
        w = np.sin(math.pi * g / T2) ** 2
        numh[i] = (complex((w * ph).sum()) * (pr / 2).conjugate()).real; denh[i] = abs(pr / 2) ** 2
        nums[i] = (complex(ph.sum()) * pr.conjugate()).real; dens[i] = abs(pr) ** 2
    out = []
    n = len(objs)
    for b in range(nb):
        s = num[:, b].sum() / den[:, b].sum()
        loo = np.array([(num[:, b].sum() - num[i, b]) / (den[:, b].sum() - den[i, b]) for i in range(n)])
        out.append((s, math.sqrt((n - 1) / n * np.sum((loo - loo.mean()) ** 2))))
    sh = nums.sum() / dens.sum(); hh = numh.sum() / denh.sum()
    loo = np.array([(nums.sum() - nums[i]) / (dens.sum() - dens[i]) - (numh.sum() - numh[i]) / (denh.sum() - denh[i]) for i in range(n)])
    return out, sh, hh, math.sqrt((n - 1) / n * np.sum((loo - loo.mean()) ** 2))
fam = {}
for key, e in P['L'].items():
    z = np.load(B + '/zeros_blind/' + key + '.npz')
    fam.setdefault(e['cls'], []).append((z['g'], e['Q'], e['mu'], e['T2'], complex(*e['pred'])))
art = []
for key, e in P['products'].items():
    if e['cls'] != 'T5' or e['rho']['c'] == 0: continue
    g = np.load(B + '/zeros_blind/' + key + '.npz')['g']
    isz = np.array([np.min(np.abs(zz - x)) < 1e-6 for x in g])
    art.append((g[~isz], e['rho']['Q'], e['rho']['mu'], e['T2'], complex(e['readings']['R_artin']['pred'], 0)))
fam['Artin'] = art
print('band (fraction of window): ' + '  '.join('%.2f-%.2f' % (edges[b], edges[b + 1]) for b in range(len(edges) - 1)))
for name in ['T1', 'T2', 'T3', 'T4a', 'T4b', 'Artin']:
    out, sh, hh, sd = bands(fam[name])
    print('%-6s ' % name + '  '.join('%.3f±%.3f' % o for o in out) + '   sharp %.4f Hann %.4f  sharp-Hann %+.4f ± %.4f (jack)' % (sh, hh, sh - hh, sd))
