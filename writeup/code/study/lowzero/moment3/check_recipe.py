"""check_recipe.py: independent check of the triple recipe density P(t) (tmlib, multilinear jets + constant term)
against the literal sum over the 2^3 swap terms at small finite shifts, with every factor evaluated directly
(no jets): zeta(1 +- 2a) by mpmath, Lam_ij(1 + z_i + z_j) by PARI lfun at real arguments, A(z) by a direct Euler
product of the residue/series R_p in plain floats, c_j by its definition.  Shifts alpha_j = beta_j = a_j = eps*(1, 1.31, 1.73)
for eps = 1e-3, 5e-4, 2.5e-4, 1.25e-4; three-level Richardson extrapolation eps -> 0.
Also: the direct per-prime A_p from the residue formula (distinct x) vs the series (tmlib's method)."""
import sys, os, math, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, mpmath as mp
import tmlib as T

pari = T.pari
PM = 10**5


def h_direct(Ds, z):
    k = len(Ds)
    val = 1.0
    for j in range(k):
        for p in T.prime_factors(Ds[j]):
            val *= (1 - p ** (-1.0 - 2 * z[j]))
    for i in range(k):
        for j in range(i + 1, k):
            Dsx = int(pari.coredisc(Ds[i] * Ds[j])); L = pari.lfuncreate(Dsx)
            s = 1 + z[i] + z[j]
            v = float(pari.lfun(L, pari(repr(s))))
            for p in T.prime_factors(Ds[i] * Ds[j]):
                c = int(pari.kronecker(Dsx, p))
                if c: v *= (1 - c * p ** (-s))
            val *= v ** 2
    # A(z): direct product over p <= PM of R_p * prod_{i,j} (1 - x_i x_j), R_p by plain series; tail by -sum_{i<j} p^{-2-2(zi+zj)}
    P = T.primes_upto(PM).astype(float)
    chis = [T.kron_primes(D, PM).astype(float) for D in Ds]
    x = [chis[j] * P ** (-0.5 - z[j]) for j in range(k)]
    e1 = x[0] + x[1] + x[2]; e2 = x[0] * x[1] + x[0] * x[2] + x[1] * x[2]; e3 = x[0] * x[1] * x[2]
    h = [np.ones_like(P), e1.copy(), e1 * e1 - e2]
    R = 1 + h[1] ** 2 + h[2] ** 2
    for n in range(3, 160):
        hn = e1 * h[-1] - e2 * h[-2] + e3 * h[-3]; h.append(hn); h.pop(0); R += hn ** 2
    lA = np.log(R)
    for i in range(k):
        for j in range(k):
            lA += np.log(1 - x[i] * x[j])
    tail = 0.0
    for i in range(k):
        for j in range(i + 1, k):
            tail -= float(mp.quad(lambda u: u ** (-2 - 2 * (z[i] + z[j])) / mp.log(u), [PM, 10 * PM, mp.inf]))
    return val * math.exp(lA.sum() + tail)


def swap_sum(Ds, a, t):
    k = len(Ds); X = [abs(D) * t / (2 * math.pi) for D in Ds]
    tot = 0.0
    for S in itertools.product((0, 1), repeat=k):
        z = [(-a[j] if S[j] else a[j]) for j in range(k)]
        term = h_direct(Ds, z)
        for j in range(k):
            if S[j]:
                term *= X[j] ** (-2 * a[j]) * float(mp.zeta(1 - 2 * a[j]))
            else:
                term *= float(mp.zeta(1 + 2 * a[j]))
        tot += term
    return tot


def residue_Ap(x):
    """A_p via the residue formula sum_k x_k^2 / (prod_i (1 - x_i x_k) prod_{j != k} (x_k - x_j)) (distinct x)."""
    R = 0.0
    for kk in range(3):
        den = 1.0
        for i in range(3):
            den *= (1 - x[i] * x[kk])
        for j in range(3):
            if j != kk:
                den *= (x[kk] - x[j])
        R += x[kk] ** 2 / den
    pr = 1.0
    for i in range(3):
        for j in range(3):
            pr *= (1 - x[i] * x[j])
    return R * pr


def main():
    rng = np.random.default_rng(5)
    print('A_p: residue formula vs series (random distinct x):')
    for _ in range(4):
        x = rng.uniform(-0.6, 0.6, 3)
        e1, e2, e3 = x.sum(), x[0] * x[1] + x[0] * x[2] + x[1] * x[2], x.prod()
        h = [1.0, e1, e1 * e1 - e2]; R = 1 + h[1] ** 2 + h[2] ** 2
        for n in range(3, 400):
            hn = e1 * h[-1] - e2 * h[-2] + e3 * h[-3]; h.append(hn); h.pop(0); R += hn ** 2
        ser = R * np.prod([1 - x[i] * x[j] for i in range(3) for j in range(3)])
        print('   x = %s  residue %.15f  series %.15f  quartic approx exp(-e2^2) %.6f' % (np.round(x, 3), residue_Ap(x), ser, math.exp(-e2 ** 2)))
    for Ds in [(1001, -1003, 1005), (-1003, 1253, -1267), (1005, 1068, -1267)]:
        tr = T.Triple(Ds, pmax=PM)
        for t in [100.0, 1000.0]:
            vals = []
            for eps in [1e-3, 5e-4, 2.5e-4, 1.25e-4]:
                a = [eps, 1.31 * eps, 1.73 * eps]
                vals.append(swap_sum(Ds, a, t))
            tab = [vals]
            for lev in range(1, 4):
                prev = tab[-1]; f = 2.0 ** lev
                tab.append([(f * prev[i + 1] - prev[i]) / (f - 1) for i in range(len(prev) - 1)])
            rr = tab[-1][0]
            print('%s t = %6.0f  swap sums eps=1e-3..1.25e-4: %s  Richardson %.6f   jets P(t) = %.6f  rel diff %.1e' % (
                Ds, t, np.round(vals, 3), rr, float(tr.P(t)), (rr - float(tr.P(t))) / float(tr.P(t))))


if __name__ == '__main__':
    main()
