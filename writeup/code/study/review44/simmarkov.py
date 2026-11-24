"""Null check of the statistic: a first-order Markov chain on the classes (P = the primes' pooled matrix) cut into blocks
with Poisson(N) transitions, statistic as in r135b (pooled P-hat per 'window', NB/(NB-1)).  Expect 1 + O(1/N)."""
import numpy as np, math, sys
rng = np.random.default_rng(12345)
SRC = open('mywin.py').read(); exec(SRC[SRC.index('def red'):SRC.index('base = small_primes')])
WIN = [6*10**9 + i*10**8 for i in range(24)]
for q in (3, 5, 7):
    r = red(q); C = len(r)
    T = sum(np.load('out/w_%d.npz' % x)['T_%d_1' % q] for x in WIN).astype(float); P = T/T.sum(1, keepdims=True)
    w, v = np.linalg.eig(P.T); pi = np.real(v[:, np.argmin(abs(w-1))]); pi /= pi.sum()
    cum = np.cumsum(P, 1)
    for Nmean, NB, reps in ((441, 10000, 6), (4407, 1000, 6)):
        vals = []
        for rep in range(reps):
            n = rng.poisson(Nmean, NB); L = n.max() + 1
            s = np.empty((NB, L), np.int64); s[:, 0] = rng.choice(C, NB, p=pi)
            for t in range(1, L):
                u = rng.random(NB); s[:, t] = (u[:, None] > cum[s[:, t-1]]).sum(1)
            # flatten into one chain per 'window' with block labels; pairs (s_t, s_t+1) for t < n_b
            idx = np.arange(L)[None, :] < n[:, None]
            a = s[:, :-1][idx[:, 1:]]; b = s[:, 1:][idx[:, 1:]]; bl = np.repeat(np.arange(NB), np.maximum(n - 1, 0))
            X = np.bincount((bl*C + a)*C + b, minlength=NB*C*C).reshape(NB, C, C).astype(float)
            Na = X.sum(2); Tt = X.sum(0); Ph = Tt/Tt.sum(1, keepdims=True); E = Na[:, :, None]*Ph[None]
            chi = ((X-E)**2/E).sum((1, 2)); cv = (X-E).sum(1)
            M = np.einsum('ba,aij->bij', Na, np.stack([np.diag(Ph[i]) - np.outer(Ph[i], Ph[i]) for i in range(C)]))
            sol = np.linalg.solve(M[:, :-1, :-1], cv[:, :-1, None])[:, :, 0]; col = (cv[:, :-1]*sol).sum(1)
            vals.append(((chi-col)/(C-1)**2).mean()*NB/(NB-1))
        vals = np.array(vals)
        print('q=%d N~%d: Markov-null within/dof = %.5f +- %.5f (%d reps of %d blocks)' % (q, Nmean, vals.mean(), vals.std(ddof=1)/math.sqrt(reps), reps, NB), flush=True)
