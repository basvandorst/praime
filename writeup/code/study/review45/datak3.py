# Own kappa3 from the zeros: u = (g/2pi) ln(|D| g/(2 pi e)); windows [x, x+l), x = arange(u0, uend - l, 0.01)
import numpy as np, json, sys, math
from scipy.special import loggamma
SCR = __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/'
LZ = __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/data/riemann/lzeros/'
ls = np.arange(0.5, 8.01, 0.5)

def k3stats(g, q, step=0.01, trim=0.0):
    u = g/(2*math.pi)*np.log(q*g/(2*math.pi*math.e))
    out = np.zeros((4, len(ls)))
    for j, l in enumerate(ls):
        x = np.arange(u[0] + trim, u[-1] - trim - l, step)
        c = (np.searchsorted(u, x + l, 'left') - np.searchsorted(u, x, 'left')).astype(float)
        m1 = c.mean(); d = c - m1; m2 = (d*d).mean(); m3 = (d**3).mean()
        out[:, j] = (m3, m2, m1, ((c - l)**3).mean())
    return out

def k3fix(g, q, lo, hi, step=0.01):
    """windows over the fixed band [u(lo), u(hi) - l): m3 with sample mean, M3 with exact mean l"""
    u = g/(2*math.pi)*np.log(q*g/(2*math.pi*math.e))
    ulo = lo/(2*math.pi)*math.log(q*lo/(2*math.pi*math.e)); uhi = hi/(2*math.pi)*math.log(q*hi/(2*math.pi*math.e))
    out = np.zeros((4, len(ls)))
    for j, l in enumerate(ls):
        x = np.arange(ulo, uhi - l, step)
        c = (np.searchsorted(u, x + l, 'left') - np.searchsorted(u, x, 'left')).astype(float)
        m1 = c.mean(); d = c - m1
        out[:, j] = ((d**3).mean(), ((c - l)**3).mean(), m1 - l, ((c - l)**2).mean())
    return out

def theta(t, q, a):
    return t/2*np.log(q/np.pi) + np.imag(loggamma((0.5 + a + 1j*t)/2))

def completeness(gall, q, a, lo, hi):
    g = np.sort(gall)
    ts = np.concatenate([[lo], np.linspace(lo, hi, 200), [hi]])
    N = np.searchsorted(g, ts, 'right')
    dev = N - theta(ts, q, a)/np.pi
    gb = g[(g >= lo) & (g <= hi)]
    u = gb/(2*math.pi)*np.log(q*gb/(2*math.pi*math.e))
    return float(np.max(np.abs(dev))), float(np.max(np.diff(u))), dev

if __name__ == '__main__':
    name = sys.argv[1]; files = sys.argv[2].split(','); lo, hi = float(sys.argv[3]), float(sys.argv[4])
    Z = {}
    for f in files:
        z = np.load(LZ + f); Z.update({k: z[k] for k in z.files})
    res = {}
    for k in Z:
        D = int(k[1:]); q = abs(D); a = 0 if D > 0 else 1
        g = np.sort(Z[k]); gb = g[(g >= lo) & (g <= hi)]
        st = k3stats(gb, q)
        st2 = k3stats(gb, q, trim=3.0)
        dmax, gmax, dev = completeness(g, q, a, lo, hi)
        fx = k3fix(gb, q, lo, hi)
        res[D] = dict(k3=st[0].tolist(), m2=st[1].tolist(), m1=st[2].tolist(), k3l=st[3].tolist(), k3trim=st2[0].tolist(),
                      nz=int(len(gb)), m3fix=fx[0].tolist(), M3fix=fx[1].tolist(), dfix=fx[2].tolist(), M2fix=fx[3].tolist(), devmax=dmax, devmean=float(np.mean(dev)), gapmax=gmax, g0=float(g[0]))
    json.dump(res, open(SCR + f'data_{name}.json', 'w'))
    bad = [(D, r['devmax'], r['gapmax']) for D, r in res.items() if r['devmax'] > 1.5 or r['gapmax'] > 3.0]
    print(name, len(res), 'zeros in band', sum(r['nz'] for r in res.values()), 'flagged', bad)
    print('max dev over sets', max(r['devmax'] for r in res.values()), 'mean of devmean', np.mean([r['devmean'] for r in res.values()]))
