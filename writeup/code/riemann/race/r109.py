"""Round 109 (preregister/round109.txt): the race variance R_chi(X, d) of round 102 from the zeros of L(s, chi) alone,
at the same 200 X in [1e6, 1e9], averaged; compared with the measured mean and the long-run law. usage: r109.py"""
import os, json, math
import numpy as np, cypari2
pari = cypari2.Pari(); pari.allocatemem(800 * 10 ** 6)
H = os.path.dirname(os.path.abspath(__file__))
r102 = json.load(open(os.path.join(H, 'r102.json')))['rows']
Xs = np.exp(np.linspace(math.log(1e6), math.log(1e9), 200))
T = 1500.0; out = {}
for name, D, q in (('chi4', -4, 4), ('chi3', -3, 3)):
    g = np.array([float(t) for t in pari.lfunzeros(pari.lfuncreate(D), T)]); g = g[g > 0]
    G = np.concatenate([g, -g]); RHO = 0.5 + 1j * G
    for d in (0.1, 0.01):
        c = ((1 + d) ** RHO - 1) / RHO
        Rz = np.zeros(len(Xs))
        for i0 in range(0, len(G), 400):
            dg = G[i0:i0 + 400, None] - G[None, :]
            cc = c[i0:i0 + 400, None] * np.conj(c[None, :]) / (2 + 1j * dg)
            for k, X in enumerate(Xs):
                Rz[k] += float(np.real(np.sum(cc * (np.exp((2 + 1j * dg) * math.log(X)) - 1)))) / (X * X)
        Rz = 2 * Rz / d
        tail = (2 + d) * (math.log(q * T / (2 * math.pi)) + 1) / (math.pi * T * d)
        Rz = Rz + tail
        meas = np.array(r102['%s_%s' % (name, d)]['R'])
        row = {'zeros_mean': float(Rz.mean()), 'measured_mean': float(meas.mean()), 'law_D': r102['%s_%s' % (name, d)]['D'],
               'diff_zeros_minus_measured': float(Rz.mean() - meas.mean()), 'max_abs_pointwise': float(np.max(np.abs(Rz - meas))),
               'corr_pointwise': float(np.corrcoef(Rz, meas)[0, 1]), 'nzeros': int(len(g)), 'tail': tail}
        out['%s_%s' % (name, d)] = row
        print(name, d, {k: round(v, 4) if isinstance(v, float) else v for k, v in row.items()}, flush=True)
json.dump(out, open(os.path.join(H, 'r109.json'), 'w'), indent=1)
