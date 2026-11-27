"""Loop 2: the limit law K_tau(u) (Mayer series, r <= 4 clusters) against Haar matrices of size N under the flow."""
import sys, numpy as np
from itertools import combinations
from scipy.stats import unitary_group

def conn_all(x):
    """x: array (..., r, r) of edge weights; returns the sum over connected graphs of prod of edge weights (vectorised, r <= 4)."""
    r = x.shape[-1]
    full = {}
    for mask in range(1, 1 << r):
        idx = [i for i in range(r) if mask >> i & 1]
        w = 1.0
        for a, b in combinations(idx, 2): w = w * (1.0 + x[..., a, b])
        full[mask] = w
    full[0] = 1.0
    C = {}
    for mask in range(1, 1 << r):
        low = mask & -mask; rest = mask ^ low; s = full[mask]; sub = rest
        while True:
            if sub != rest: s = s - C[low | sub] * full[mask ^ (low | sub)]
            if sub == 0: break
            sub = (sub - 1) & rest
        C[mask] = s
    return C[(1 << r) - 1]

def K_limit(u, tau, rmax=4, n=200000, seed=0):
    """u^2 e^{-2 tau u} sum_r (1/r!) int_{sum a = u} prod e^{2 tau a_i^2}/a_i Phi_r^2; terms r = 1..rmax.
    r >= 2 by Monte Carlo on the simplex (uniform Dirichlet), integrand vanishes like a_i at the faces."""
    rng = np.random.default_rng(seed)
    terms = [u * np.exp(-2 * tau * u * (1 - u))]
    fact = 1.0
    for r in range(2, rmax + 1):
        fact *= r
        a = rng.dirichlet(np.ones(r), size=n) * u          # uniform on the simplex, volume u^{r-1}/(r-1)!
        x = np.expm1(2 * tau * a[:, :, None] * a[:, None, :])
        phi = conn_all(x)
        f = np.exp(2 * tau * np.sum(a * a, axis=1)) / np.prod(a, axis=1) * phi ** 2
        vol = u ** (r - 1) / np.prod(np.arange(1, r))
        terms.append(u * u * np.exp(-2 * tau * u) * vol * f.mean() / fact)
    return terms

def simulate(N, t, n_mat, ks, seed=1):
    rng = np.random.default_rng(seed)
    j = np.arange(N + 1); w = np.exp(-t * j * (N - j))
    vals = np.zeros((n_mat, len(ks)))
    for m, U in enumerate(unitary_group.rvs(N, size=n_mat, random_state=rng)):
        z = np.linalg.eigvals(U)
        r = np.roots(np.poly(-z) * w)
        ang = np.angle(-r)
        vals[m] = np.abs(np.exp(1j * np.outer(ks, ang)).sum(axis=1)) ** 2
    return vals.mean(axis=0) / N, vals.std(axis=0) / np.sqrt(n_mat) / N

if __name__ == "__main__":
    N = int(sys.argv[1]); n_mat = int(sys.argv[2])
    ks = np.array([int(round(u * N)) for u in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.25, 1.5, 2.0)])
    for tau in [float(a) for a in sys.argv[3:]]:
        mean, err = simulate(N, tau / N ** 2, n_mat, ks)
        print(f"N = {N}, tau = {tau}, {n_mat} matrices")
        print("   u    simulated K      series r<=4   terms r = 1, 2, 3, 4")
        for k, m, e in zip(ks, mean, err):
            u = k / N
            if u < 1:
                T = K_limit(u, tau)
                print(f" {u:5.2f} {m:8.4f}+-{e:6.4f}  {sum(T):8.4f}    " + " ".join(f"{x:8.5f}" for x in T))
            else:
                print(f" {u:5.2f} {m:8.4f}+-{e:6.4f}")
