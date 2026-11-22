import sys, json, math, numpy as np
sys.path.insert(0, '/home/bas/praime/writeup/code/study')
import r47
PP = r47.PP
d = json.load(open('/home/bas/praime/writeup/data/study/r81.json'))
LAM = d['lambda']; DMAX = d['DMAX']; NZ = 48
KS = (1, 2, 3, 4); QS = (3, 5, 7, 11, 13)
Z = np.exp(2j*np.pi*np.arange(NZ)/NZ)
VAR = ['full', 'noR', 'flat_p', 'cramer_rough']
dist = {v: {k: np.zeros(DMAX//2) for k in KS} for v in VAR}
def lam_for(v): return LAM
for D in range(2, DMAX+1, 2):
    st = r47.setup(D); M = st['M']; n = M.shape[1]
    for v in VAR:
        rho0 = 1/(LAM*r47.RHO_DEN)
        s1 = st['s1'] if v in ('full', 'noR') else np.ones(n)
        S2 = st['S2'] if v != 'cramer_rough' else 1.0
        p = rho0*s1
        logG = np.zeros((M.shape[0], NZ), complex)
        if n:
            logG += M @ np.log1p(p[:, None]*(Z[None, :]-1))
        if n > 1 and v == 'full':
            W = p[:, None]*(Z[None, :]-1)/(1+p[:, None]*(Z[None, :]-1))
            for m in range(NZ):
                X = M*W[:, m][None, :]
                V = X.real @ st['R'] + 1j*(X.imag @ st['R'])
                logG[:, m] += 0.5*np.sum(V*X, axis=1)
        coef = np.fft.fft(np.exp(logG), axis=1).real/NZ
        for k in KS:
            dist[v][k][D//2-1] = rho0*S2*coef[:, k-1].sum()
Ds = np.arange(2, DMAX+1, 2)
out = {}
for v in VAR:
    out[v] = {}
    print('==', v)
    for k in KS:
        dk = dist[v][k]
        row = []
        for q in QS:
            th = dk[Ds % q == 0].sum()/dk.sum()
            dd = d['data'][str(k)][str(q)]
            z = (th-dd['p'])/dd['se']; row.append(z)
            out[v][f'{k},{q}'] = [th, dd['p'], z]
        print(k, ' '.join('%+7.1f' % z for z in row), ' sumz2 %.1f' % sum(z*z for z in row))
    print('total sum z^2 k>=2:', sum(out[v][f'{k},{q}'][2]**2 for k in (2,3,4) for q in QS), ' all:', sum(x[2]**2 for x in out[v].values()))
json.dump(out, open('/home/bas/scratch/skeptic-A/variants.json', 'w'), indent=1)
