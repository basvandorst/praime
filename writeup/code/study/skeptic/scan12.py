import sys, json, numpy as np
sys.path.insert(0, '/home/bas/praime/writeup/code/study')
import r47
d = json.load(open(sys.argv[1]))
L0 = d['lambda']; DMAX = d['DMAX']; NZ = 48
KS = tuple(int(k) for k in d['data']); QS = tuple(int(q) for q in d['data']['1'])
Z = np.exp(2j*np.pi*np.arange(NZ)/NZ); Ds = np.arange(2, DMAX+1, 2)
sts = [r47.setup(D) for D in Ds]
def run(v, sh):
    lam = L0+sh; rho0 = 1/(lam*r47.RHO_DEN)
    dist = {k: np.zeros(len(Ds)) for k in KS}
    for t, st in enumerate(sts):
        M = st['M']; n = M.shape[1]
        p = rho0*(st['s1'] if v in ('noR','full') else np.ones(n))
        logG = M @ np.log1p(p[:, None]*(Z[None, :]-1)) if n else np.zeros((M.shape[0], NZ), complex)
        if v == 'full' and n > 1:
            W = p[:, None]*(Z[None, :]-1)/(1+p[:, None]*(Z[None, :]-1))
            for m in range(NZ):
                X = M*W[:, m][None, :]
                V = X.real @ st['R'] + 1j*(X.imag @ st['R'])
                logG[:, m] += 0.5*np.sum(V*X, axis=1)
        coef = np.fft.fft(np.exp(logG), axis=1).real/NZ
        for k in KS: dist[k][t] = rho0*st['S2']*coef[:, k-1].sum()
    return np.array([[(dist[k][Ds % q == 0].sum()/dist[k].sum()-d['data'][str(k)][str(q)]['p'])/d['data'][str(k)][str(q)]['se'] for q in QS] for k in KS])
def rep(name, zs):
    sm = zs[:, :5]; print(name, 'q<=13: sumz2 %.1f max %.2f | all: sumz2 %.1f max %.2f' % ((sm**2).sum(), abs(sm).max(), (zs**2).sum(), abs(zs).max()), flush=True)
rep('full sh0', run('full', 0.0))
for v in ('noR', 'flat_p'):
    rep(v+' sh0', run(v, 0.0))
    best = min(((run(v, s)[:, :5]**2).sum(), s) for s in np.round(np.arange(-1.0, 1.01, 0.1), 2))
    zs = run(v, best[1]); rep(v+' best sh %.1f' % best[1], zs); print(np.round(zs, 1))
