# Is the O(1/N^2) CUE correction to every n-point function a "random rescaling" (mixture) correction?
# f(rho) = det[ sin(pi rho x_ij)/(pi x_ij) ]; mixture: delta R_n = (sig2/2) f''(1); CUE(N): K_N = S + D/(6N^2), D = -pi x sin(pi x)*(-1)?? check
import numpy as np
rng = np.random.default_rng(1)
def Kr(x, rho):
    with np.errstate(invalid='ignore', divide='ignore'):
        v = np.sin(np.pi*rho*x)/(np.pi*x)
    v[x == 0] = rho
    return v
def f(xs, rho):
    X = xs[:, None]-xs[None, :]
    return np.linalg.det(Kr(X, rho))
def KN(x, N):
    with np.errstate(invalid='ignore', divide='ignore'):
        v = np.sin(np.pi*x)/(N*np.sin(np.pi*x/N))
    v[x == 0] = 1.0
    return v
for n in (2, 3, 4, 5):
    for trial in range(3):
        xs = np.sort(rng.uniform(0, 3, n))
        h = 1e-3
        fpp = (f(xs, 1+h) - 2*f(xs, 1) + f(xs, 1-h))/h**2
        # CUE: R_n(N) = det K_N ; extract 1/N^2 coefficient by Richardson
        X = xs[:, None]-xs[None, :]
        Ns = np.array([200., 400.])
        vals = [np.linalg.det(KN(X, N)) for N in Ns]
        base = f(xs, 1.0)
        c = [(v-base)*N**2 for v, N in zip(vals, Ns)]
        c_extrap = (4*c[1]-c[0])/3
        # mixture prediction with sig2 = -1/(6N^2): coefficient of 1/N^2 is -(1/12) f''
        print(n, trial, 'CUE 1/N^2 coef %.6f   mixture -(1/12)f\'\' %.6f' % (c_extrap, -fpp/12))
