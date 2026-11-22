import sys, json, numpy as np
sys.path.insert(0, '/home/bas/praime/writeup/code/study')
import r47
d = json.load(open('/home/bas/praime/writeup/data/study/r81.json'))
L0 = d['lambda']; DMAX = d['DMAX']; NZ = 48; KS=(1,2,3,4); QS=(3,5,7,11,13)
Z = np.exp(2j*np.pi*np.arange(NZ)/NZ); Ds = np.arange(2, DMAX+1, 2)
shifts = np.round(np.arange(-1.0, 1.01, 0.1), 2)
sts = [r47.setup(D) for D in Ds]
for v in ('noR', 'flat_p'):
    best = None
    for sh in shifts:
        lam = L0+sh; rho0 = 1/(lam*r47.RHO_DEN)
        dist = {k: np.zeros(len(Ds)) for k in KS}
        for t, st in enumerate(sts):
            M = st['M']; n = M.shape[1]
            s1 = st['s1'] if v == 'noR' else np.ones(n)
            p = rho0*s1
            logG = M @ np.log1p(p[:, None]*(Z[None, :]-1)) if n else np.zeros((M.shape[0], NZ))
            coef = np.fft.fft(np.exp(logG), axis=1).real/NZ
            for k in KS: dist[k][t] = rho0*st['S2']*coef[:, k-1].sum()
        zs = np.array([[(dist[k][Ds % q == 0].sum()/dist[k].sum()-d['data'][str(k)][str(q)]['p'])/d['data'][str(k)][str(q)]['se'] for q in QS] for k in KS])
        s = (zs**2).sum()
        if best is None or s < best[0]: best = (s, sh, zs)
    print(v, 'best shift', best[1], 'sumz2', round(best[0],1), 'max|z|', round(abs(best[2]).max(),2)); print(np.round(best[2],1))
