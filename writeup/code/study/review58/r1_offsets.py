"""Review 152, step 1: own theta, own sums, offsets per class / root number, jackknife; offset as a function of T."""
import json, math, os, sys
import numpy as np
from scipy.special import loggamma
B4 = '/home/bas/praime/writeup/code/study/lowzero/universal/blind4/'
B1 = '/home/bas/praime/writeup/code/study/lowzero/universal/blind/'

def theta(t, Q, mu):
    t = np.asarray(t, float)
    out = t / 2 * np.log(Q / np.pi ** len(mu))
    for m in mu:
        out = out + np.imag(loggamma((0.5 + m + 1j * t) / 2))
    return out

def load():
    rows = []
    P = json.load(open(B4 + 'predictions_blind4.json'))['L']
    for k, e in P.items():
        z = np.load(B4 + 'zeros_blind4/' + k + '.npz'); info = json.loads(str(z['info']))
        rows.append(dict(key=k, rnd=152, cls=e['cls'], Q=e['Q'], mu=e['mu'], eps=e['eps'][0], c=e['c'], T2=e['T2'], g=z['g'], info=info, spec=e['spec']))
    P1 = json.load(open(B1 + 'predictions_blind.json'))['L']
    for k, e in P1.items():
        if e['cls'] not in ('T4a', 'T4b'): continue
        z = np.load(B1 + 'zeros_blind/' + k + '.npz'); info = json.loads(str(z['info']))
        rows.append(dict(key=k, rnd=151, cls={'T4a': 'sym3', 'T4b': 'g2'}[e['cls']], Q=e['Q'], mu=e['mu'], eps=e['eps'][0], c=e['c'], T2=e['T2'], g=z['g'], info=info, spec=e['spec']))
    return rows

def jk_mean(x):
    x = np.asarray(x, float); n = len(x); loo = (x.sum() - x) / (n - 1)
    return x.mean(), math.sqrt((n - 1) / n * np.sum((loo - loo.mean()) ** 2))

if __name__ == '__main__':
    rows = load()
    for r in rows:
        g = r['g']; g = g[(g > 0) & (g <= r['T2'])]
        r['ph'] = np.exp(2j * theta(g, r['Q'], r['mu'])); r['gg'] = g
        r['G'] = r['ph'].sum(); r['pred'] = -r['eps'] * r['c'] * r['T2'] / (2 * math.pi)
        r['res'] = r['G'] - r['pred']
    # compare with the project's results
    RB = json.load(open(B4 + 'results_blind4.json'))['per_L']
    md = max(abs(complex(*RB[r['key']]['G']) - r['G']) for r in rows if r['rnd'] == 152)
    print('max |G_own - G_project| round 152:', md)
    for rnd in (152, 151):
        for cls in ('g2', 'sym3', 'rs'):
            for sg in (None, 1, -1):
                S = [r for r in rows if r['rnd'] == rnd and r['cls'] == cls and (sg is None or r['eps'] == sg)]
                if len(S) < 3: continue
                m, s = jk_mean([r['res'].real for r in S]); mi, si = jk_mean([r['res'].imag for r in S])
                print('round %d %-5s eps %-4s n=%2d  Re(G-pred) %+6.2f +- %5.2f   Im %+6.2f +- %5.2f   mean c %.3f  mean T2 %.0f' % (
                    rnd, cls, sg, len(S), m, s, mi, si, np.mean([r['c'] for r in S]), np.mean([r['T2'] for r in S])))
    # slopes per class (registered estimator) for reference
    def slope(S):
        sx = sum((r['G'] * r['pred']).real for r in S); sy = sum(r['pred'] ** 2 for r in S); return sx / sy
    G2 = [r for r in rows if r['cls'] == 'g2']
    print('genus-2 slope round 152 all %.3f, eps+1 %.3f, eps-1 %.3f' % (slope([r for r in G2 if r['rnd'] == 152]),
          slope([r for r in G2 if r['rnd'] == 152 and r['eps'] > 0]), slope([r for r in G2 if r['rnd'] == 152 and r['eps'] < 0])))
    # offset as a function of T, genus 2, 55 curves
    Ts = np.array([5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150], float)
    D = np.zeros((len(G2), len(Ts)), complex)
    for i, r in enumerate(G2):
        for j, T in enumerate(Ts):
            D[i, j] = r['ph'][r['gg'] <= T].sum() + r['eps'] * r['c'] * T / (2 * math.pi)
    print('\nGenus 2 (55 curves): mean residual Re/Im vs T (jackknife), and by eps')
    for j, T in enumerate(Ts):
        m, s = jk_mean(D[:, j].real); mi, si = jk_mean(D[:, j].imag)
        ep = np.array([r['eps'] for r in G2])
        mp, sp = jk_mean(D[ep > 0, j].real); mm, sm = jk_mean(D[ep < 0, j].real)
        print('T=%5.0f  Re %+6.2f +- %4.2f  Im %+6.2f +- %4.2f   eps+1 %+6.2f +- %4.2f  eps-1 %+6.2f +- %4.2f   Re/(T/2pi) %+.3f' % (T, m, s, mi, si, mp, sp, mm, sm, m / (T / 2 / math.pi)))
    np.savez('r1_D.npz', D=D, Ts=Ts,
             keys=np.array([r['key'] for r in G2]), eps=np.array([r['eps'] for r in G2]), c=np.array([r['c'] for r in G2]), Q=np.array([r['Q'] for r in G2]))
    # fits of the mean residual (Re) on T: models a*T/2pi, a*sqrt(T), a (constant), a + b T
    m = D.real.mean(0)
    # covariance across curves (curves independent)
    C = np.cov(D.real.T) / len(G2)
    Ci = np.linalg.pinv(C)
    for name, X in (('T/2pi', np.c_[Ts / 2 / math.pi]), ('sqrtT', np.c_[np.sqrt(Ts)]), ('const', np.c_[np.ones_like(Ts)]),
                    ('logT', np.c_[np.log(Ts)]), ('const+T/2pi', np.c_[np.ones_like(Ts), Ts / 2 / math.pi]),
                    ('sqrtT+T/2pi', np.c_[np.sqrt(Ts), Ts / 2 / math.pi])):
        A = X.T @ Ci @ X; b = np.linalg.solve(A, X.T @ Ci @ m); res = m - X @ b; chi2 = res @ Ci @ res
        print('fit %-12s coef %s  se %s  chi2 %.1f (dof %d)' % (name, np.round(b, 3), np.round(np.sqrt(np.diag(np.linalg.inv(A))), 3), chi2, len(Ts) - X.shape[1]))
