import sys, numpy as np
from loop15 import c_exact
from loop9 import cth
rng = np.random.default_rng(21); N = int(sys.argv[1]); M = int(sys.argv[2]); cuts = [float(v) for v in sys.argv[3].split(',')]
kall = np.arange(1, int(max(cuts)*N) + 1)
call = np.array([float(c_exact(N, int(k))) if k <= N else float(cth(np.array([k/N]))[0]) for k in kall])
err = {c: [] for c in cuts}
for it in range(M):
    Z = (rng.standard_normal((N, N)) + 1j*rng.standard_normal((N, N)))/np.sqrt(2); Q, R = np.linalg.qr(Z)
    x = np.sort(np.angle(np.linalg.eigvals(Q*(np.diag(R)/np.abs(np.diag(R)))[None, :])) % (2*np.pi))*N/(2*np.pi)
    rho = np.exp(-2j*np.pi*np.outer(kall, x)/N).sum(axis=1); n = np.arange(N)
    for c in cuts:
        K = int(c*N); ks = kall[:K]; coef = call[:K]*rho[:K]/(2j*np.pi*ks)
        P = lambda y: 2*np.real(np.exp(2j*np.pi*np.outer(y, ks)/N) @ coef)
        const = np.mean(x + P(x) - n); y = n + const
        for _ in range(60): y = 0.5*y + 0.5*(n + const - P(y))
        err[c].append(np.mean((y - x)**2))
for c in cuts: print('matrices N %d, modes to u = %.1f: rms error of a rebuilt eigenvalue %.4f spacings' % (N, c, np.sqrt(np.mean(err[c]))), flush=True)
