"""Referee's own implementation of Psi_{mu nu} for the Lewin-Lieb-Seiringer reduction.

Definition (paper eq. 19): Phi_{mu nu}(a,b) = a^3 b^3 (1 - W), W = iint dmu(u) dnu(v) / |e1 + u/a - v/b|,
i.e. W = Coulomb interaction of two unit charges: mu dilated to scale R1 = 1/a centred at 0 and
nu dilated to scale R2 = 1/b centred at -e1 (equivalently +e1 by symmetry).
Psi = Phi_{mu nu} + Phi_{nu mu} - Phi_{nu nu}.

Method (written here, not taken from the repo):
 * Newton: potential of a unit sphere of radius R at distance x is 1/max(x,R);
   potential of a uniform unit ball of radius R is (3R^2-x^2)/(2R^3) for x<=R, 1/x beyond.
 * Archimedes: for y uniform on the sphere of radius rho about e1, |y|^2 is uniform on
   [(1-rho)^2,(1+rho)^2] (width 4 rho).  So the shell average of V(|y|) is
   (P(thi) - P(tlo))/(4 rho), P = antiderivative of t -> V(sqrt t).
 * Sphere-sphere: closed form from that.  Ball-ball: shells of the second ball, Gauss-Legendre on the
   pieces between breakpoints (integrand polynomial on each piece, so the rule is exact up to rounding).
"""
import numpy as np

# ---------------- sphere - sphere ----------------
def W_ss(R1, R2):
    """Interaction of unit-mass uniform spheres radius R1 (centre 0) and R2 (centre e1)."""
    R1 = np.asarray(R1, float); R2 = np.asarray(R2, float)
    tlo = (1.0 - R2) ** 2
    thi = (1.0 + R2) ** 2
    def P(t):
        # antiderivative of 1/max(sqrt t, R1): t/R1 below R1^2, 2 sqrt t - R1 above
        return np.where(t <= R1 * R1, t / R1, 2.0 * np.sqrt(t) - R1)
    return (P(thi) - P(tlo)) / (4.0 * R2)


def Phi_sphere(a, b, r, q):
    """Phi_{sigma_r sigma_q}(a,b) for spheres of radii r (first, scale a) and q (second, scale b)."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    with np.errstate(divide='ignore', invalid='ignore'):
        w = W_ss(r / a, q / b)
        out = a ** 3 * b ** 3 * (1.0 - w)
    return np.where((a > 0) & (b > 0), out, 0.0)

# ---------------- ball - ball ----------------
_GX, _GW = np.polynomial.legendre.leggauss(4)

def _P_ball(t, R1):
    # antiderivative in t of V_ball(sqrt t): (3R1^2 t - t^2/2)/(2R1^3) for t<=R1^2 ; 2 sqrt t - 3R1/4 beyond
    return np.where(t <= R1 * R1, (3.0 * R1 * R1 * t - 0.5 * t * t) / (2.0 * R1 ** 3),
                    2.0 * np.sqrt(t) - 0.75 * R1)


def W_BB(R1, R2):
    """Interaction of unit-mass uniform balls radius R1 (centre 0) and R2 (centre e1), shells of ball 2."""
    R1, R2 = np.broadcast_arrays(np.asarray(R1, float), np.asarray(R2, float))
    shp = R1.shape
    R1 = R1.ravel(); R2 = R2.ravel()
    # breakpoints in y = rho/R2 in [0,1]: rho = |1-R1|, 1+R1, R1-1, 1
    bps = np.stack([np.zeros_like(R1), np.abs(1 - R1) / R2, (1 + R1) / R2, (R1 - 1) / R2,
                    1.0 / R2, np.ones_like(R1)], axis=1)
    bps = np.clip(bps, 0.0, 1.0)
    bps.sort(axis=1)
    total = np.zeros_like(R1)
    for k in range(bps.shape[1] - 1):
        lo = bps[:, k]; hi = bps[:, k + 1]
        half = 0.5 * (hi - lo); mid = 0.5 * (hi + lo)
        for x, w in zip(_GX, _GW):
            y = mid + half * x
            rho = R2 * y
            with np.errstate(divide='ignore', invalid='ignore'):
                avg = (_P_ball((1 + rho) ** 2, R1) - _P_ball((1 - rho) ** 2, R1)) / (4.0 * rho)
            avg = np.where(rho > 0, avg, 0.0)
            total += w * half * 3.0 * y * y * avg
    return total.reshape(shp)


def Psi_ball(a, b):
    a, b = np.broadcast_arrays(np.asarray(a, float), np.asarray(b, float))
    out = np.zeros(a.shape)
    m = (a > 0) & (b > 0)
    out[m] = a[m] ** 3 * b[m] ** 3 * (1.0 - W_BB(1.0 / a[m], 1.0 / b[m]))
    return out

# ---------------- mixtures of spheres ----------------
class Mixture:
    """mu = sum_i mw[i] sigma_{mr[i]},  nu = sum_j nw[j] sigma_{nr[j]}."""
    def __init__(self, mr, mw, nr, nw):
        self.mr = np.asarray(mr, float); self.mw = np.asarray(mw, float)
        self.nr = np.asarray(nr, float); self.nw = np.asarray(nw, float)
        # list of (first radius r, second radius q, coefficient): Psi = sum c * Phi_{s_r s_q}
        pairs = {}
        def add(r, q, c):
            pairs[(r, q)] = pairs.get((r, q), 0.0) + c
        for r, wr in zip(self.mr, self.mw):
            for q, wq in zip(self.nr, self.nw):
                add(r, q, wr * wq)          # Phi_{mu nu}
                add(q, r, wq * wr)          # Phi_{nu mu}
        for r, wr in zip(self.nr, self.nw):
            for q, wq in zip(self.nr, self.nw):
                add(r, q, -wr * wq)         # -Phi_{nu nu}
        self.pairs = [(r, q, c) for (r, q), c in pairs.items() if c != 0.0]

    def Psi(self, a, b):
        a, b = np.broadcast_arrays(np.asarray(a, float), np.asarray(b, float))
        out = np.zeros(a.shape)
        for r, q, c in self.pairs:
            out += c * Phi_sphere(a, b, r, q)
        return out

    def D_mumu(self):
        # D(sigma_r, sigma_q) = 1/(2 max(r,q))
        R = np.maximum.outer(self.mr, self.mr)
        return 0.5 * float(self.mw @ (1.0 / R) @ self.mw)
