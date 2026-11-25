import json, math, numpy as np, sys
sys.path.insert(0, '.')
from r1_offsets import load, theta, jk_mean
rows = [r for r in load() if r['cls'] == 'g2']
out = []
for r in rows:
    g = r['g']; g = g[(g > 0) & (g <= r['T2'])]
    ph = np.exp(2j * theta(g, r['Q'], r['mu'])); G = ph.sum(); pred = -r['eps'] * r['c'] * r['T2'] / (2 * math.pi)
    i = r['info']
    out.append(dict(key=r['key'], rnd=r['rnd'], N=int(r['Q']), eps=int(r['eps']), m0=i['m0'], c=r['c'], pred=pred, G=[G.real, G.imag], res=[(G - pred).real, (G - pred).imag],
                    n=len(g), nfirst=i['n_first'], added=i['added'], extra=i['added_extra'], dup=i['duplicates'], drift=i['drift'], rmean=i['r_mean'], cexp=i['c_exp'],
                    g1=float(g[0]), spec=r['spec']))
    o = out[-1]
    print('%-14s r%d N=%5d eps %+d m0 %d c %+.3f pred %+7.2f G %+7.2f %+6.2f res %+6.2f n %3d first %3d add %d/%d dup %d drift %+.3f rmean %+.3f (cexp %+.2f) g1 %.3f' % (
        o['key'], o['rnd'], o['N'], o['eps'], o['m0'], o['c'], o['pred'], G.real, G.imag, o['res'][0], o['n'], o['nfirst'], o['added'], o['extra'], o['dup'], o['drift'], o['rmean'], o['cexp'], o['g1']))
json.dump(out, open('r2_table.json', 'w'), indent=1)
m0 = np.array([o['m0'] for o in out]); res = np.array([o['res'][0] for o in out]); eps = np.array([o['eps'] for o in out])
for v in sorted(set(m0)):
    s = res[m0 == v]; print('m0 = %d: n = %d, mean res %+.2f +- %.2f' % (v, len(s), s.mean(), s.std(ddof=1) / math.sqrt(len(s)) if len(s) > 1 else float('nan')))
