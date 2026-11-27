"""Reviewer SDP for t*(r): exact moment relaxation (no angle grid), own barrier (path-following Newton) solver.

Unknown per prime p <= X: trigonometric moments c_{p,l} = E cos(l theta), l = 1..L_p, of ONE probability
measure on [0, pi]; the moment set is exactly {Toeplitz [c_|i-j|]_{i,j=0..L_p} PSD} (Caratheodory-Toeplitz,
symmetrised). Forms k in ks: M_k(c) - k t I PSD (conductor e^{-kt}); split into even/odd Legendre blocks.
maximise t.  The barrier optimum t_mu satisfies t_mu <= t* <= t_mu + mu * (total block size).
"""
import sys; sys.path.insert(0, "/home/bas/praime/writeup/code/study/bootstrap/maass/review")
import numpy as np
import weilmat as W


def build(geo, r, eps, ks=(1, 2, 3, 4), parity_split=True, central=False):
    """central=True: odd L(f) has root number -1, so (GRH) a zero at gamma = 0 of odd order; use Q - h(0) >= 0,
    h(0) = (int phi)^2 = 2a c_0^2 in the orthonormal Legendre basis."""
    X = geo.X
    pps = W.prime_powers(X)
    primes = sorted({p for p, m in pps})
    Lp = {p: max(k * m for (q, m) in pps if q == p for k in ks) for p in primes}
    var = []  # (p, l)
    for p in primes:
        for l in range(1, Lp[p] + 1):
            var.append((p, l))
    var.append(("t", 0))
    vidx = {v: i for i, v in enumerate(var)}
    nv = len(var)
    blocks = []  # list of (C, F) with F shape (nv, n, n)
    N = geo.N
    idx_sets = [np.arange(0, N, 2), np.arange(1, N, 2)] if parity_split else [np.arange(N)]
    for k in ks:
        C = W.arch(geo, W.gammas(k, r, eps))
        if central and eps == 1 and k == 1:
            C = C.copy(); C[0, 0] -= 2 * geo.a
        F = np.zeros((nv, N, N))
        for p, m in pps:
            w = 2 * np.log(p) * p ** (-m / 2)
            G = geo.G(m * np.log(p))
            for l, mult in W.sym_chars(k, m).items():
                if l == 0:
                    C = C - w * mult * G
                else:
                    F[vidx[(p, l)]] -= w * mult * G
        F[vidx[("t", 0)]] = -k * np.eye(N)
        for ix in idx_sets:
            blocks.append((C[np.ix_(ix, ix)], F[:, ix][:, :, ix].copy(), ("form", k)))
    for p in primes:
        n = Lp[p] + 1
        C = np.eye(n)
        F = np.zeros((nv, n, n))
        for l in range(1, n):
            E = np.eye(n, k=l) + np.eye(n, k=-l)
            F[vidx[(p, l)]] = E
        blocks.append((C, F, ("toep", p)))
    return var, blocks


def Fval(blk, x):
    C, F, _ = blk
    return C + np.tensordot(x, F, axes=(0, 0))


def solve(var, blocks, mu_end=1e-11, verbose=False):
    nv = len(var)
    obj = np.zeros(nv); obj[-1] = 1.0
    x = np.zeros(nv)
    # starting t: strictly feasible
    tmin = np.inf
    for blk in blocks:
        if blk[2][0] == "form":
            k = blk[2][1]
            tmin = min(tmin, np.linalg.eigvalsh(blk[0])[0] / k)
    x[-1] = tmin - 1.0
    ntot = sum(b[0].shape[0] for b in blocks)
    mu = 1.0
    newton_total = 0

    def phi(x, mu):
        s = obj @ x / mu
        for blk in blocks:
            try:
                L = np.linalg.cholesky(Fval(blk, x))
            except np.linalg.LinAlgError:
                return -np.inf
            s += 2 * np.sum(np.log(np.diag(L)))
        return s

    while True:
        for it in range(200):
            g = obj / mu
            H = np.zeros((nv, nv))
            for blk in blocks:
                Fm = Fval(blk, x)
                L = np.linalg.cholesky(Fm)
                Li = np.linalg.inv(L)
                S = np.einsum("ab,ibc,dc->iad", Li, blk[1], Li, optimize=True)
                g += np.einsum("iaa->i", S)
                Sf = S.reshape(nv, -1)
                H -= Sf @ Sf.T
            # Newton direction (H negative definite)
            try:
                d = np.linalg.solve(-H, g)
            except np.linalg.LinAlgError:
                d = np.linalg.lstsq(-H, g, rcond=None)[0]
            lam2 = g @ d
            newton_total += 1
            if lam2 < 1e-12:
                break
            s = 1.0
            f0 = phi(x, mu)
            while s > 1e-12:
                xn = x + s * d
                fn = phi(xn, mu)
                if fn > f0 + 0.25 * s * lam2:
                    break
                s *= 0.5
            x = xn
        if verbose:
            print("mu %.1e t %.12f newton %d" % (mu, x[-1], newton_total))
        if mu * ntot < mu_end:
            break
        mu /= 6.0
    gap = mu * ntot
    return x, gap


def tstar(geo, r, eps, ks=(1, 2, 3, 4), verbose=False, central=False):
    var, blocks = build(geo, r, eps, ks, central=central)
    x, gap = solve(var, blocks, verbose=verbose)
    mom = {}
    for (p, l), v in zip(var, x):
        if p != "t":
            mom.setdefault(p, {})[l] = v
    ap = {p: 2 * mom[p][1] for p in mom}
    # independent check of the returned point: lambda_min of each form at the returned moments
    momarr = {p: np.array([1.0] + [mom[p][l] for l in range(1, len(mom[p]) + 1)]) for p in mom}
    lam = {}
    for k in ks:
        lam[k] = np.linalg.eigvalsh(W.form_matrix(geo, k, r, eps, momarr))[0] / k
    toepmin = min(np.linalg.eigvalsh(np.array([[momarr[p][abs(i - j)] for j in range(len(momarr[p]))]
                                                  for i in range(len(momarr[p]))]))[0] for p in momarr)
    return dict(t=x[-1], gap=gap, ap=ap, lam=lam, toepmin=toepmin, mom=momarr)
