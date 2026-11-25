"""Own implementation of the swapped-term covariance in equal windows (independent of swaplib.py).
R_swap(t1, t2) = (c/4pi^2) { kappa e^{2i(theta_2(t2) - theta_1(t1))} H(t2 - t1) + c.c. }   (registered c = 1)
C_swap(w) = mean_t int int_[0,w]^2 R_swap(t + y1, t + y2) dy1 dy2
          = (c/2pi^2) mean_t Re[ kappa e^{i(Th(tc) - w dL/2)} I(Lm(tc)) ],  I = (1/2) int_{-w}^{w} e^{i u Lm} H(u) V(|u|) du,
V(a) = int_a^{2w-a} e^{i v dL/2} dv  (theta linearised at the window centre tc = t + w/2; dL = ln(q2/q1), Lm = ln(sqrt(q1 q2) tc/2pi)).
Two evaluations of I: 'full' (trapezoid in u, H from an Euler product on a u-grid) and 'asym' (the endpoint/kink
expansion to order Lm^-3, which uses H, H' at u = 0, +-w only and so is blind to the spectral content of H near the
frequency Lm, i.e. to the Dirichlet coefficients of H at n ~ sqrt(q1 q2) t/2pi)."""
import math
import numpy as np
from mylib import theta, T1, T2, H as HSTEP, grid


def hp_factors(psi, P):
    """returns arrays for the Euler factors of H restricted to p not | k1 k2 (psi != 0)."""
    keep = psi != 0
    return P[keep].astype(float), psi[keep].astype(float)


def H_and_dH(u, Pf, s):
    """H(u)/kappa (Euler product over Pf with psi values s) and H'(u)/kappa, u array. Handles the zero factor at p = 3, psi = -1, u = 0."""
    u = np.atleast_1d(np.asarray(u, float)); lp = np.log(Pf)
    e = np.exp(1j * np.outer(u, lp))                                 # p^{iu}
    num = 1 - 2 / Pf + s * e / Pf; den = (1 - s * np.conj(e) / Pf) * (1 - 1 / Pf) ** 2
    dnum = 1j * lp * s * e / Pf; dden_over = -(1j * lp * s * np.conj(e) / Pf) / (1 - s * np.conj(e) / Pf)   # d/du log(1 - s p^{-1-iu})^{-1}
    # factor-by-factor: h = num/den ; dlog h = dnum/num + dden_over  (num may vanish only for p = 3, s = -1, u = 0)
    is3 = (Pf == 3)
    h = num / den
    Hrest = np.prod(np.where(is3[None, :], 1.0, h), axis=1)
    dl_rest = np.sum(np.where(is3[None, :], 0.0, dnum / np.where(is3[None, :], 1.0, num) + dden_over), axis=1)
    if is3.any():
        j = np.nonzero(is3)[0][0]
        h3 = h[:, j]; dh3 = (dnum[:, j] / den[:, j]) + h3 * dden_over[:, j]
    else:
        h3 = np.ones(len(u), complex); dh3 = np.zeros(len(u), complex)
    Hv = h3 * Hrest; dH = dh3 * Hrest + h3 * Hrest * dl_rest
    return Hv, dH


def pair_setup(D1, D2, W):
    q1, q2 = abs(D1), abs(D2); a1 = 0 if D1 > 0 else 1; a2 = 0 if D2 > 0 else 1
    dL = math.log(q2 / q1)
    out = []
    for w in W:
        tg = grid(w)[::5]                     # t-step 0.1 for the average (phases vary slowly)
        tc = tg + w / 2
        Th = 2 * theta(tc, q2, a2) - 2 * theta(tc, q1, a1)
        Lm = np.log(math.sqrt(q1 * q2) * tc / (2 * math.pi))
        out.append((tc, Th, Lm))
    return dL, out


def cswap_asym(kappa, dL, setup, W, H0, dH0, Hw, dHw, Hmw, dHmw, order=3):
    """order-2 (1/L^2) or order-3 (adds 1/L^3) endpoint expansion; H values per width arrays (16,)"""
    res = np.zeros(len(W))
    for k, w in enumerate(W):
        tc, Th, Lm = setup[k]
        V0 = (np.exp(1j * dL * w) - 1) / (1j * dL / 2) if abs(dL) > 1e-12 else 2 * w
        Vp0 = -np.exp(1j * dL * w) - 1; Vpw = -2 * np.exp(1j * dL * w / 2)
        eL = np.exp(1j * w * Lm)
        S = (Hw[k] * Vpw * eL - 2 * H0 * Vp0 + Hmw[k] * Vpw / eL) / Lm ** 2
        if order >= 3:
            S = S + 1j * (2 * dHw[k] * Vpw * eL - 4 * dH0 * Vp0 + 2 * dHmw[k] * Vpw / eL) / Lm ** 3
        I = 0.5 * S
        ph = np.exp(1j * (Th - w * dL / 2))
        res[k] = np.mean((kappa * ph * I).real) / (2 * math.pi ** 2)
    return res


def cswap_full(kappa, dL, setup, W, Pf, s, du=0.02, nL=120, Hconst=None):
    """full trapezoid integral in u with H from the Euler product (or a constant Hconst)"""
    res = np.zeros(len(W)); wmax = max(W)
    ug = np.arange(0, wmax + du / 2, du)
    if Hconst is None:
        Hu, _ = H_and_dH(ug, Pf, s)
    else:
        Hu = np.full(len(ug), Hconst, complex)
    for k, w in enumerate(W):
        tc, Th, Lm = setup[k]
        sel = ug <= w + 1e-9; u = ug[sel]; Hs = Hu[sel]
        a = u; b = 2 * w - u
        V = (np.exp(1j * dL * b / 2) - np.exp(1j * dL * a / 2)) / (1j * dL / 2) if abs(dL) > 1e-12 else (b - a).astype(complex)
        wq = np.full(len(u), du); wq[0] = du / 2; wq[-1] = du / 2
        Lg = np.linspace(Lm.min(), Lm.max(), nL)
        E = np.exp(1j * np.outer(Lg, u))
        # I(L) = (1/2)[ int_0^w e^{iuL} H(u) V du + int_0^w e^{-iuL} conj(H(u)) V du ]  (H(-u) = conj H(u), real psi)
        Jg = 0.5 * (E @ (Hs * V * wq) + np.conj(E) @ (np.conj(Hs) * V * wq))
        I = np.interp(Lm, Lg, Jg.real) + 1j * np.interp(Lm, Lg, Jg.imag)
        ph = np.exp(1j * (Th - w * dL / 2))
        res[k] = np.mean((kappa * ph * I).real) / (2 * math.pi ** 2)
    return res
