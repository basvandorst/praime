# Stage 2 float optimisation: mu = sum_j mu_j sigma_{j/K}, nu = sum_m w_m sigma_{s_m} (M free spheres); f = g.
import numpy as np, sys, time, json
from scipy.optimize import minimize


def T(rho, rhop, A, B):
    """rho^3 rho'^3 Phi_ss(a/rho, b/rho'): sphere rho at 0 (scale a), sphere rho' at e1 (scale b)"""
    X1 = A*rhop + B*rho - A*B
    X2 = np.abs(A*rhop - B*rho) - A*B
    return A*A*B*B/(4*rho*rhop)*(np.maximum(X1, 0)**2 - np.maximum(X2, 0)**2)


def make_grid(N1=150, N2=60, X=30.0):
    x1 = np.linspace(0, 3, N1 + 1)[1:]
    x2 = 3*np.exp(np.linspace(0, np.log(X/3), N2 + 1)[1:])
    xs = np.concatenate([x1, x2])
    e = np.concatenate([[0.0], xs])
    w = (e[1:] - e[:-1])/2
    w[:-1] += (e[2:] - e[1:-1])/2
    return xs, w, X


xs, wq, X = make_grid()
A = xs[:, None]; B = xs[None, :]
TRIL = np.tril(np.ones((len(xs), len(xs)), dtype=bool))


def value(radii, mu, nu, return_g=False):
    """radii: atoms; mu, nu: weights on atoms (sum 1 each)"""
    n = len(radii)
    c = np.outer(mu, nu) + np.outer(nu, mu) - np.outer(nu, nu)
    Psi = np.zeros((len(xs), len(xs)))
    for j in range(n):
        for k in range(n):
            if abs(c[j, k]) > 0:
                Psi += c[j, k]*T(radii[j], radii[k], A, B)
    d = np.maximum(np.diag(Psi), 0)/2
    H = np.where(TRIL, Psi - d[None, :], -np.inf)
    g = H.max(axis=1)
    bb = np.linspace(0, 1, 20001)
    bV = (mu[None, :]*np.minimum(1.0, bb[:, None]/radii[None, :])).sum(1)
    kap = (bb**3*(1 - bV)).max()
    J = 4*np.pi*((wq*g/xs**5).sum() + kap/X)
    D = 0.5*(mu[:, None]*mu[None, :]/np.maximum(radii[:, None], radii[None, :])).sum()
    cval = 1.5*(2*J*D*D)**(1/3)
    return (cval, J, D, kap, g) if return_g else cval


if __name__ == '__main__':
    K = int(sys.argv[1]); M = int(sys.argv[2])
    grid_r = np.arange(1, K + 1)/K

    def unpack(p):
        em = np.exp(p[:K] - p[:K].max()); mu = em/em.sum()
        s = 1/(1 + np.exp(-p[K:K + M]))
        if M > 1:
            ew = np.exp(p[K + M:] - p[K + M:].max()); wn = ew/ew.sum()
        else:
            wn = np.ones(1)
        radii = np.concatenate([grid_r, s])
        return radii, np.concatenate([mu, np.zeros(M)]), np.concatenate([np.zeros(K), wn])

    TSTART = time.time()
    best = [1e9, None]

    def obj(p):
        v = value(*unpack(p))
        if v < best[0]:
            best[0], best[1] = v, p.copy()
        if time.time() - TSTART > float(sys.argv[3]) if len(sys.argv) > 3 else 500:
            raise StopIteration
        return v
    # start: ball-like mu, nu sphere(s) near 0.69
    mu0 = (np.arange(1, K + 1)**3 - np.arange(0, K)**3)/K**3
    s0 = np.array([0.69, 0.95, 0.4][:M])
    p0 = np.concatenate([np.log(mu0), np.log(s0/(1 - s0)), np.zeros(M if M > 1 else 0)])
    if M > 1:
        p0[K + M:] = np.log(np.array([0.9, 0.05, 0.05][:M]))
    t = time.time()
    print('start', obj(p0), time.time() - t, flush=True)
    if len(sys.argv) > 4:
        dj = json.load(open(sys.argv[4]))
        if dj['K'] == K and dj['M'] == M:
            p0 = np.array(dj['p'])
        else:   # interpolate the cumulative mass of mu onto the new grid, keep the nu radius
            rj = np.array(dj['radii'][:dj['K']]); mj = np.array(dj['mu'][:dj['K']])
            cdf = np.concatenate([[0.0], np.cumsum(mj)]); rr0 = np.concatenate([[0.0], rj])
            Fn = np.interp(grid_r, rr0, cdf); mun = np.diff(np.concatenate([[0.0], Fn])); mun = np.maximum(mun, 1e-9)
            s_old = np.array(dj['radii'][dj['K']:dj['K'] + M])
            p0 = np.concatenate([np.log(mun/mun.sum()), np.log(s_old/(1 - s_old))])
        print('restart from', sys.argv[4], obj(p0), flush=True)
    try:
        res = minimize(obj, p0, method='L-BFGS-B', options={'maxiter': 400, 'maxfun': 40000})
        print('L-BFGS-B', res.fun, res.nit, res.message, time.time() - t, flush=True)
        res2 = minimize(obj, res.x, method='Powell', options={'maxiter': 20000, 'xtol': 1e-6, 'ftol': 1e-10})
        print('Powell', res2.fun, time.time() - t, flush=True)
    except StopIteration:
        print('time limit', flush=True)
    print('best', best[0], time.time() - t, flush=True)
    radii, mu, nu = unpack(best[1])
    np.set_printoptions(precision=5, suppress=True, linewidth=150)
    print('radii', radii); print('mu', mu); print('nu', nu)
    cval, J, D, kap, g = value(radii, mu, nu, True)
    print('c', cval, 'J', J, 'D', D, 'kappa', kap)
    json.dump({'p': best[1].tolist(), 'K': K, 'M': M, 'radii': radii.tolist(), 'mu': mu.tolist(), 'nu': nu.tolist(), 'c_float_grid': cval},
              open(f's2_opt_K{K}_M{M}.json', 'w'))
