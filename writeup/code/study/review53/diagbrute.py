"""Brute-force 'diagonal-only' models for a triple at height t:
  prodAFE(t)  = 2 sum_{N <= sqrt(C)} a(N)^2 / N,  a = chi1 * chi2 * chi3 (full Dirichlet convolution), C = X1 X2 X3
                (the registered H_diag density, sharp cutoff; checks its residue formula)
  perfactor(t) = 8 sum_N b(N)^2 / N, b(N) = sum_{m1 m2 m3 = N, m_j <= sqrt(X_j)} chi1(m1) chi2(m2) chi3(m3)
                (diagonal of the product of the three per-factor approximate functional equations, sharp cutoffs)
usage: diagbrute.py D1 D2 D3 t1 [t2 ...]   prints both densities and the recipe P(t) and H_diag(t) of my own code."""
import sys, os, math, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import myrecipe as M

pari = M.pari


def chi_vec(D, N):
    v = pari('vector(%d, n, kronecker(%d, n))' % (N, D))
    return np.concatenate([[0.0], np.array([int(x) for x in v], dtype=float)])   # index n


def prod_afe(Ds, t):
    X = [abs(D) * t / (2 * math.pi) for D in Ds]
    Nmax = int(math.floor(math.sqrt(X[0] * X[1] * X[2])))
    ch = [chi_vec(D, Nmax) for D in Ds]
    # a = chi1 * chi2 * chi3 up to Nmax by slicing
    def conv(f, g):
        out = np.zeros(Nmax + 1)
        for m in range(1, Nmax + 1):
            if f[m] != 0:
                out[m::m] += f[m] * g[1:Nmax // m + 1]
        return out
    a = conv(conv(ch[0], ch[1]), ch[2])
    n = np.arange(1, Nmax + 1)
    return 2 * np.sum(a[1:] ** 2 / n)


def per_factor(Ds, t, B=4_000_000):
    X = [abs(D) * t / (2 * math.pi) for D in Ds]
    x = [int(math.floor(math.sqrt(XX))) for XX in X]
    ch = [chi_vec(D, xx) for D, xx in zip(Ds, x)]
    # c(M) = sum_{m1 m2 = M, m1 <= x1, m2 <= x2} chi1 chi2
    m1 = np.arange(1, x[0] + 1); m2 = np.arange(1, x[1] + 1)
    Mm = (m1[:, None] * m2[None, :]).ravel(); W = (ch[0][1:][:, None] * ch[1][1:][None, :]).ravel()
    Mmax = x[0] * x[1]
    c = np.bincount(Mm, weights=W, minlength=Mmax + 1)
    Ntot = Mmax * x[2]
    tot = 0.0
    for N0 in range(1, Ntot + 1, B):
        N1 = min(N0 + B, Ntot + 1)
        b = np.zeros(N1 - N0)
        for m3 in range(1, x[2] + 1):
            if ch[2][m3] == 0:
                continue
            lo = -(-N0 // m3); hi = min(Mmax, (N1 - 1) // m3)
            if lo > hi:
                continue
            b[m3 * lo - N0: m3 * hi - N0 + 1: m3] += ch[2][m3] * c[lo:hi + 1]
        n = np.arange(N0, N1)
        tot += np.sum(b * b / n)
    return 8 * tot


if __name__ == '__main__':
    Ds = tuple(int(a) for a in sys.argv[1:4]); ts = [float(a) for a in sys.argv[4:]]
    hm, _ = M.multilinear(Ds); g = M.diag_coeffs(Ds)
    for t in ts:
        P = M.P_at(Ds, hm, t); Dd = M.diag_density(Ds, g, t)
        lead = hm[()] * np.prod([math.log(abs(D) * t / (2 * math.pi)) for D in Ds])
        pa = prod_afe(Ds, t) if math.sqrt(np.prod([abs(D) * t / (2 * math.pi) for D in Ds])) < 3e7 else float('nan')
        pf = per_factor(Ds, t)
        print(json.dumps({'Ds': Ds, 't': t, 'P': P, 'Hdiag_formula': Dd, 'prodAFE_brute': pa, 'perfactor_brute': pf, 'lead': lead,
                          'prodAFE/P': pa / P, 'Hdiag/P': Dd / P, 'perfactor/P': pf / P, 'lead/P': lead / P}), flush=True)
