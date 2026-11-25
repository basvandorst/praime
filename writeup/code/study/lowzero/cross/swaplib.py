"""The swapped (Gram-lattice) term of the cross two-point function (derivation.md section 3.3-3.4, 4.5) for two real
primitive characters chi_D1, chi_D2, turned into the equal-window covariance on the registered grid:
  R_swap(t1, t2) = (1/4pi^2) { kappa e^{2i(theta_2(t2) - theta_1(t1))} H(t2 - t1) + c.c. },  kappa = (k1/phi(k1))(k2/phi(k2)),
  H(u) = L(1-iu, psi) L(1+iu, psi) A12(-iu) = prod_p (1 - c1/p - c2/p + psi p^(-1+iu)) / ((1 - psi p^(-1-iu))(1 - c1/p)(1 - c2/p))
  (c_j = 1[p not | k_j]; Euler product truncated at PMAX; real psi so L(1+iu) = conj L(1-iu)).
  C_swap(w) = mean over the grid t of (1/2pi^2) Re[ kappa int int_{[0,w]^2} e^{i(2 theta_2(t+y2) - 2 theta_1(t+y1))} H(y2-y1) dy ].
The double integral is done exactly in (u, v) = (y2 - y1, y1 + y2) with theta_j linearised around the window centre
(theta_j(t+y) = theta_j(t+w/2) + (y - w/2) L_j/2; the curvature term (y-w/2)^2/(4t) <= 0.03 rad is dropped)."""
import math
import numpy as np
from crosslib import theta, primes_upto, kron_table

def phi_ratio(k):
    r = 1.0; n = k; p = 2
    while p * p <= n:
        if n % p == 0:
            r *= p / (p - 1)
            while n % p == 0: n //= p
        p += 1
    if n > 1: r *= n / (n - 1)
    return r

def H_of_u(D1, D2, u, P, x1, x2):
    """H(u) on an array u >= 0 (H(-u) = conj H(u)); x1, x2 = chi_D(p) rows for primes P. Only primes not dividing
    k1 k2 enter the product (for p | k1 k2 the factor is 1 after the conductor factors k_j/phi(k_j))."""
    ps = (x1 * x2).astype(float); keep = ps != 0
    Pf = P[keep].astype(float); s = ps[keep]; lp = np.log(Pf)
    e = np.exp(1j * np.outer(np.asarray(u, float), lp))                   # (nu, np)
    num = 1 - 2 / Pf + s * e / Pf; den = (1 - s * np.conj(e) / Pf) * (1 - 1 / Pf) ** 2
    out = np.exp(np.log(num / den).sum(1))
    return phi_ratio(abs(D1)) * phi_ratio(abs(D2)) * out

def cswap(D1, D2, W, T1, T2, h, P, x1, x2, du=0.02, tstep=0.24, nL=96):
    """predicted swapped-term covariance at the window widths W (array). The window integral depends on t only through
    Lm(t) = ln(sqrt(q1 q2) t/2pi) (dL = ln(q2/q1) is constant): it is tabulated on nL values of Lm and interpolated."""
    q1, q2 = abs(D1), abs(D2); a1 = 0 if D1 > 0 else 1; a2 = 0 if D2 > 0 else 1
    umax = float(np.max(W)); ug = np.arange(0, umax + du / 2, du); Hu = H_of_u(D1, D2, ug, P, x1, x2)
    dL = math.log(q2 / q1); out = np.zeros(len(W))
    Lg = np.linspace(math.log(math.sqrt(q1 * q2) * T1 / (2 * math.pi)), math.log(math.sqrt(q1 * q2) * T2 / (2 * math.pi)), nL)
    for k, w in enumerate(W):
        M = int(math.floor((T2 - w - T1) / h + 1e-9)) + 1; tg = T1 + h * np.arange(0, M, max(1, int(round(tstep / h))))
        tc = tg + w / 2
        Th = 2 * theta(tc, q2, a2) - 2 * theta(tc, q1, a1)
        Lm = np.log(math.sqrt(q1 * q2) * tc / (2 * math.pi))
        sel = ug <= w + 1e-12; u = ug[sel]; Hs = Hu[sel]
        a = u; b = 2 * w - u
        V = (b - a) if abs(dL) < 1e-12 else (np.exp(1j * dL * b / 2) - np.exp(1j * dL * a / 2)) / (1j * dL / 2)
        wq = np.full(len(u), du); wq[0] = du / 2
        if abs(u[-1] - w) < 1e-9: wq[-1] = du / 2
        E = np.exp(1j * np.outer(Lg, u))
        Jg = 0.5 * ((E * (Hs * V * wq)[None, :]).sum(1) + (np.conj(E) * (np.conj(Hs) * V * wq)[None, :]).sum(1))
        I = np.interp(Lm, Lg, Jg.real) + 1j * np.interp(Lm, Lg, Jg.imag)
        ph = np.exp(1j * (Th - (w / 2) * dL))
        out[k] = float(np.mean((ph * I).real)) / (2 * math.pi ** 2)
    return out
