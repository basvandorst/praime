"""De Bruijn-Newman flow of zeta zeros, in heights gamma:
   d gamma_k/dt = (1/2) sum_{j != k} 1/(gamma_k - gamma_j), sum over all zeros and mirrors -gamma_j.
First M zeros move; zeros above Y(t) = Y0 - (pi/8) t are a smooth density (1/2pi) ln(y/2pi), both signs.
Force by a block multipole method (near blocks direct, far blocks by expansion); see check_fmm.py."""
import numpy as np
from scipy.special import spence, comb

Z1 = '/home/bas/praime/writeup/data/riemann/odlyzko_zeros1.txt'


def load(M):
    g = np.loadtxt(Z1, skiprows=2, max_rows=M + 1)
    return g[:M].copy(), 0.5 * (g[M - 1] + g[M])


def tail_force(x, Y):
    """(1/2) int_Y^inf rho(y) [1/(x-y) + 1/(x+y)] dy, rho = ln(y/2pi)/(2pi), 0 < x < Y (closed form)."""
    w = Y / x
    A = np.log(x) * (np.log(w - 1) - np.log(w + 1)) - spence(w) - np.log(w) * np.log1p(w) - spence(1 + w)
    I_ln = np.pi ** 2 / 3 - A                      # int_Y^inf ln y [1/(y-x) - 1/(y+x)] dy
    I_c = np.log((Y + x) / (Y - x))                # int_Y^inf [1/(y-x) - 1/(y+x)] dy
    return -(I_ln - np.log(2 * np.pi) * I_c) / (4 * np.pi)


class Force:
    def __init__(self, M, Y0, B=60, P=24, eta=0.4):
        assert M % B == 0
        self.M, self.Y0, self.B, self.P, self.eta = M, Y0, B, P, eta
        self.nb = M // B
        p = np.arange(P)
        self.sgn = (-1.0) ** p
        # C[p, q] = binom(p+q, q)
        self.C = comb(p[:, None] + p[None, :], p[None, :])

    def pair_sum(self, x):
        """sum_{j != k} [1/(x_k - x_j)] + sum_j 1/(x_k + x_j) for all k (mirror sources included)."""
        B, P, nb = self.B, self.P, self.nb
        Xb = x.reshape(nb, B)
        lo, hi = Xb[:, 0], Xb[:, -1]
        c = 0.5 * (lo + hi); r = 0.5 * (hi - lo)
        U = Xb - c[:, None]
        Mt = np.empty((nb, P)); pw = np.ones_like(U)
        for q in range(P):
            Mt[:, q] = pw.sum(1); pw = pw * U
        cs = np.concatenate([c, -c]); rs = np.concatenate([r, r])
        Ms = np.concatenate([Mt, Mt * self.sgn[None, :]])
        D = c[:, None] - cs[None, :]                        # (nb, 2nb)
        far = (r[:, None] + rs[None, :]) <= self.eta * np.abs(D)
        Dinv = np.where(far, 1.0 / np.where(far, D, 1.0), 0.0)
        # G[t, p, n] = sum_s Ms[s, p] Dinv[t, s]^(n+1)
        L = np.zeros((nb, P))
        Dn = Dinv.copy()
        for n in range(P):
            Gn = Dn @ Ms[:, :n + 1]                         # (nb, n+1): p = 0..n, q = n - p
            q = n - np.arange(n + 1)
            L[:, q] += Gn * (self.C[np.arange(n + 1), q] * self.sgn[q])[None, :]
            Dn = Dn * Dinv
        # evaluate local expansions
        out = np.zeros_like(U); pw = np.ones_like(U)
        for q in range(P):
            out += L[:, q:q + 1] * pw; pw = pw * U
        # near pairs, direct
        ti, si = np.nonzero(~far)
        Sb = np.concatenate([Xb, -Xb])
        diff = Xb[ti][:, :, None] - Sb[si][:, None, :]
        selfp = np.nonzero(ti == si)[0]
        ar = np.arange(B)
        diff[selfp[:, None], ar[None, :], ar[None, :]] = np.inf
        contrib = (1.0 / diff).sum(2)
        np.add.at(out, ti, contrib)
        return out.ravel()

    def __call__(self, x, t):
        Y = self.Y0 - np.pi * t / 8
        return 0.5 * self.pair_sum(x) + tail_force(x, Y)


def integrate(x0, Y0, tsnap, c, dtmax, log=None):
    F = Force(len(x0), Y0)
    x = x0.copy(); t = 0.0; snaps = [x.copy()]; nsteps = 0; neval = 0
    for tn in tsnap[1:]:
        while t < tn - 1e-15:
            gmin = np.min(np.diff(x))
            dt = min(dtmax, c * gmin ** 2, tn - t)
            k1 = F(x, t); k2 = F(x + 0.5 * dt * k1, t + 0.5 * dt)
            k3 = F(x + 0.5 * dt * k2, t + 0.5 * dt); k4 = F(x + dt * k3, t + dt)
            x = x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            t += dt; nsteps += 1; neval += 4
        snaps.append(x.copy())
        if log: log('t=%.3f steps=%d gmin=%.4f' % (t, nsteps, np.min(np.diff(x))))
    return np.array(snaps), nsteps
