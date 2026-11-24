"""Own Sigma^2 from zeros, plus completeness check against smooth N(t)."""
import numpy as np, json, sys
from scipy.special import loggamma
OUT = __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/'
LZ = __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/data/riemann/lzeros/'
d = json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/bk/predictions_per_lfunction.json'))
lgrid = np.array(d['l'])

def unfold(g, k):
    return (g / (2 * np.pi)) * np.log(k * g / (2 * np.pi * np.e))

def theta_count(T, k, a):
    # (1/pi) * Im ln[(k/pi)^{s/2} Gamma((s+a)/2)] at s=1/2+iT  -> smooth count of zeros in (0, T]
    s = 0.5 + 1j * T
    ph = (T / 2) * np.log(k / np.pi) + np.imag(loggamma((s + a) / 2))
    return ph / np.pi

def sigma2(u, ls, step=0.01, u0=None, u1=None, ddof=0):
    if u0 is None: u0 = u[0]
    if u1 is None: u1 = u[-1]
    out = []; nw = []
    for l in ls:
        x = np.arange(u0, u1 - l, step)
        c = np.searchsorted(u, x + l) - np.searchsorted(u, x)
        out.append(np.var(c, ddof=ddof)); nw.append(len(x))
    return np.array(out), np.array(nw)

if __name__ == '__main__':
    which = sys.argv[1]
    if which == 'B':
        z = np.load(LZ + 'r134_B.npz')
        Ds = d['B']['D']
        S = []; info = []
        for D in Ds:
            g = z['D%d' % D]
            g = np.sort(g)
            a = 0 if D > 0 else 1
            # completeness: count in (0, T] vs theta/pi; drift
            Ts = np.linspace(20, 600, 59)
            cnt = np.searchsorted(g, Ts)
            sm = theta_count(Ts, abs(D), a)
            dev = cnt - sm
            gg = g[(g >= 100) & (g <= 600)]
            u = unfold(gg, abs(D))
            u0 = unfold(100.0, abs(D)); u1 = unfold(600.0, abs(D))
            s, nw = sigma2(u, lgrid, u0=u0, u1=u1)
            S.append(s)
            # min gap (unfolded) as a check for close pairs
            info.append((D, len(gg), float(dev.mean()), float(dev.min()), float(dev.max()), float(np.diff(u).min()), float(np.diff(u).max()), len(g)))
        np.savez(OUT + 'S_B.npz', S=np.array(S), D=np.array(Ds))
        json.dump(info, open(OUT + 'info_B.json', 'w'))
        dv = np.array([i[2] for i in info])
        print('mean dev', dv.mean(), 'range of mean dev', dv.min(), dv.max())
        print('max |dev|', max(max(abs(i[3]), abs(i[4])) for i in info))
        print('min gap', min(i[5] for i in info), 'max gap', max(i[6] for i in info))
        print('zeros in band: mean', np.mean([i[1] for i in info]))
    else:
        # test A: complex characters, both sides
        S = []; Ns = []; info = []
        keys = []
        for part in range(4):
            z = np.load(LZ + 'r121_part%d.npz' % part)
            for kk in z.keys():
                keys.append((kk, z[kk]))
        kd = dict(keys)
        names = list(kd.keys())
        # order to match prediction N list: name 'N101_m2'
        for N in d['A']['N']:
            cand = [n for n in names if n.startswith('N%d_' % N)]
            assert len(cand) == 1, (N, cand)
            g = np.sort(kd[cand[0]])
            ss = []; ws = []
            for side in (1, -1):
                h = np.sort(side * g[(side * g >= 100) & (side * g <= 500)])
                u = unfold(h, N)
                s, nw = sigma2(u, lgrid, u0=unfold(100.0, N), u1=unfold(500.0, N))
                ss.append(s); ws.append(nw)
                Ts = np.linspace(20, 500, 49)
                hh = np.sort(side * g[side * g > 0])
                cnt = np.searchsorted(hh, Ts)
                info.append((cand[0], side, float((cnt - Ts / (2 * np.pi) * np.log(N * Ts / (2 * np.pi * np.e))).mean()), len(h)))
            ss = np.array(ss); ws = np.array(ws)
            S.append((ss * ws).sum(0) / ws.sum(0))
            Ns.append(N)
        np.savez(OUT + 'S_A.npz', S=np.array(S), N=np.array(Ns))
        json.dump(info, open(OUT + 'info_A.json', 'w'))
        dv = np.array([i[2] for i in info])
        print('A mean dev', dv.mean(), dv.min(), dv.max())
