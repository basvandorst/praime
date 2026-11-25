"""Own loader for the 100 genus-2 L-functions (rounds 151, 152 and the review-58 replication) and the 75 Sym^3 / 20 RS
of rounds 151-152 (for comparison). Own theta (Gamma_R convention, PARI mu = Vga + (k-1)/2), own sums."""
import json, math, os, glob
import numpy as np
from scipy.special import loggamma
B4 = '/home/bas/praime/writeup/code/study/lowzero/universal/blind4/'
B1 = '/home/bas/praime/writeup/code/study/lowzero/universal/blind/'
RV = '/home/bas/praime/writeup/code/study/review58/rep/'

def theta(t, Q, mu):
    t = np.asarray(t, float)
    out = t / 2 * np.log(Q / np.pi ** len(mu))
    for m in mu:
        out = out + np.imag(loggamma((0.5 + m + 1j * t) / 2))
    return out

def load(classes=('g2',)):
    rows = []
    P = json.load(open(B4 + 'predictions_blind4.json'))['L']
    for k, e in P.items():
        if e['cls'] not in classes: continue
        z = np.load(B4 + 'zeros_blind4/' + k + '.npz'); info = json.loads(str(z['info']))
        rows.append(dict(key=k, rnd='152', cls=e['cls'], Q=e['Q'], mu=e['mu'], eps=e['eps'][0], c=e['c'], T2=e['T2'],
                         g=np.array(z['g']), m0=info.get('m0'), spec=e['spec'], status=info['status']))
    P1 = json.load(open(B1 + 'predictions_blind.json'))['L']
    for k, e in P1.items():
        cl = {'T4a': 'sym3', 'T4b': 'g2'}.get(e['cls'])
        if cl not in classes: continue
        z = np.load(B1 + 'zeros_blind/' + k + '.npz'); info = json.loads(str(z['info']))
        rows.append(dict(key=k, rnd='151', cls=cl, Q=e['Q'], mu=e['mu'], eps=e['eps'][0], c=e['c'], T2=e['T2'],
                         g=np.array(z['g']), m0=info.get('m0'), spec=e['spec'], status=info.get('status')))
    if 'g2' in classes:
        for fn in sorted(glob.glob(RV + 'R_N*.npz')):
            z = np.load(fn); info = json.loads(str(z['info']))
            rows.append(dict(key=info['key'], rnd='rep', cls='g2', Q=info['Q'], mu=info['mu'], eps=info['eps'][0], c=info['c'],
                             T2=150, g=np.array(z['g']), m0=info.get('m0'), spec=info['spec'], status=info['status']))
    for r in rows:
        g = r['g']; g = g[(g > 0) & (g <= r['T2'])]; r['gg'] = g
        r['ph'] = np.exp(2j * theta(g, r['Q'], r['mu']))
        r['G'] = complex(r['ph'].sum())
    return rows

def slope(G, pred):
    """registered estimator b = sum Re(G conj p)/sum |p|^2 and its leave-one-out jackknife error."""
    G = np.asarray(G, complex); pred = np.asarray(pred, complex)
    x = (G * pred.conj()).real; y = np.abs(pred) ** 2; sx, sy = x.sum(), y.sum(); n = len(x)
    loo = (sx - x) / (sy - y); se = math.sqrt((n - 1) / n * np.sum((loo - loo.mean()) ** 2))
    return sx / sy, se

if __name__ == '__main__':
    rows = load()
    print('genus-2 rows', len(rows), {k: sum(r['rnd'] == k for r in rows) for k in ('151', '152', 'rep')},
          'status', set(r['status'] for r in rows), 'T2', set(r['T2'] for r in rows))
    print('mu sets', set(tuple(r['mu']) for r in rows))
    RB = json.load(open(B4 + 'results_blind4.json'))['per_L']
    print('max |G_own - G_project| (round 152):', max(abs(complex(*RB[r['key']]['G']) - r['G']) for r in rows if r['rnd'] == '152'))
    for sub, S in (('all', rows), ('151+152', [r for r in rows if r['rnd'] != 'rep']), ('rep', [r for r in rows if r['rnd'] == 'rep']),
                   ('eps+1', [r for r in rows if r['eps'] > 0]), ('eps-1', [r for r in rows if r['eps'] < 0])):
        b, se = slope([r['G'] for r in S], [-r['eps'] * r['c'] * r['T2'] / (2 * math.pi) for r in S])
        print('%-8s n=%3d slope %.3f +- %.3f' % (sub, len(S), b, se))
    # weight concentration
    w = np.array([r['c'] ** 2 for r in rows]); ws = np.sort(w)[::-1]
    print('effective n (sum w)^2/sum w^2 = %.1f ; top-10 weight share %.2f' % (w.sum() ** 2 / (w ** 2).sum(), ws[:10].sum() / w.sum()))
    print('c range %.3f .. %.3f ; eps+1 mean c %.3f, eps-1 mean c %.3f' % (min(r['c'] for r in rows), max(r['c'] for r in rows),
          np.mean([r['c'] for r in rows if r['eps'] > 0]), np.mean([r['c'] for r in rows if r['eps'] < 0])))
