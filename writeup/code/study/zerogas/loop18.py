import sys, numpy as np
from loop15 import c_exact
rng = np.random.default_rng(18); N = int(sys.argv[1]); M = int(sys.argv[2]); K = int(0.9*N)
ck = np.array([float(c_exact(N, k)) for k in range(1, K + 1)]); ks = np.arange(1, K + 1); kh = np.arange(1, 2*N + 1)
num = np.zeros(len(kh), complex); d1 = np.zeros(len(kh)); d2 = np.zeros(len(kh)); res = []
for it in range(M):
    Z = (rng.standard_normal((N, N)) + 1j*rng.standard_normal((N, N)))/np.sqrt(2); Q, R = np.linalg.qr(Z)
    x = np.sort(np.angle(np.linalg.eigvals(Q*(np.diag(R)/np.abs(np.diag(R)))[None, :])) % (2*np.pi))*N/(2*np.pi)
    rho = np.exp(-2j*np.pi*np.outer(ks, x)/N).sum(axis=1); coef = ck*rho/(2j*np.pi*ks)      # the low modes, with the response coefficients
    P = lambda y: 2*np.real(np.exp(2j*np.pi*np.outer(y, ks)/N) @ coef)
    # S(x_n) = n - x_n + const: fix the constant from the true mean (one number), then solve y + P(y) = n + const
    n = np.arange(N); const = np.mean(x + P(x) - n)
    y = n + const
    for _ in range(60): y = n + const - P(y)*0.7 + 0.3*(y - (n + const - P(y)))*0 if False else 0.5*y + 0.5*(n + const - P(y))
    res.append(np.sqrt(np.mean((y - x)**2)))
    a = np.exp(-2j*np.pi*np.outer(kh, x)/N).sum(axis=1); b = np.exp(-2j*np.pi*np.outer(kh, y)/N).sum(axis=1)
    num += a*np.conj(b); d1 += np.abs(a)**2; d2 += np.abs(b)**2
corr = np.real(num)/np.sqrt(d1*d2); u = kh/N
print('N %d, %d matrices, modes used: k <= %d (u <= 0.9). rms error of a rebuilt position: %.4f spacings' % (N, M, K, np.mean(res)))
for lo in np.arange(0.0, 2.0, 0.1):
    sel = (u > lo) & (u <= lo + 0.1)
    print('  u in (%.1f, %.1f]: correlation of rebuilt and true modes %.3f | power of rebuilt / true %.3f' % (lo, lo + 0.1, np.real(num[sel]).sum()/np.sqrt(d1[sel].sum()*d2[sel].sum()), d2[sel].sum()/d1[sel].sum()))
