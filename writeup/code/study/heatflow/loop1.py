"""Loop 1: Haar unitary matrices under the backward heat flow of the characteristic polynomial.
E_t(x) = sum_j exp(-t j (N-j)) e_j x^j ; p_k(t) = power sums of the flowed roots (normalised to the circle).
Compare E|p_k(t)|^2 with (a) the exact partition (Mayer) formula, (b) its one-part term k exp(-2 t k (N-k))."""
import sys, numpy as np
from itertools import combinations
from math import factorial
from scipy.stats import unitary_group

def partitions(k, m=None):
    m = m or k
    if k == 0:
        yield (); return
    for i in range(min(k, m), 0, -1):
        for r in partitions(k - i, i):
            yield (i,) + r

def conn_weight(parts, t):
    """sum over connected graphs on the parts of prod_edges (exp(2 t i j) - 1), by the recursion
    C(S) = W(S) - sum_{T subset S containing the first vertex, T != S} C(T) W(S \\ T), W(S) = prod_{pairs in S} exp(2 t i j)."""
    r = len(parts)
    W = {}; C = {}
    for mask in range(1, 1 << r):
        idx = [i for i in range(r) if mask >> i & 1]
        W[mask] = np.exp(2 * t * sum(parts[a] * parts[b] for a, b in combinations(idx, 2)))
    W[0] = 1.0
    for mask in range(1, 1 << r):
        low = mask & -mask
        rest = mask ^ low
        s = W[mask]
        sub = rest
        while True:                       # proper subsets T = low | sub2 with sub2 strict subset of rest
            if sub != rest:
                T = low | sub
                s -= C[T] * W[mask ^ T]
            if sub == 0: break
            sub = (sub - 1) & rest
        C[mask] = s
    return C[(1 << r) - 1]

def exact(k, N, t):
    """E|p_k(t)|^2 for k <= N: p_k(t)/k = e^{-tNk} sum_lambda (sign) prod_s e^{t i_s^2} Conn(lambda) prod_i (p_i/i)^{m_i}/m_i!"""
    tot = 0.0; one = 0.0
    for lam in partitions(k):
        coef = k * np.exp(t * sum(i * i for i in lam)) * conn_weight(lam, t)
        mult = {}
        for i in lam: mult[i] = mult.get(i, 0) + 1
        den = 1.0
        for i, m in mult.items(): den *= i ** m * factorial(m)
        c = coef / den                    # coefficient of p_lambda (up to sign)
        z = den                           # z_lambda = prod i^m m!
        term = c * c * z * np.exp(-2 * t * N * k)
        tot += term
        if len(lam) == 1: one = term
    return tot, one

def simulate(N, t, n_mat, kmax, seed=1):
    rng = np.random.default_rng(seed)
    j = np.arange(N + 1)
    w = np.exp(-t * j * (N - j))
    acc = np.zeros(kmax + 1); acc2 = np.zeros(kmax + 1); offc = 0.0
    for U in unitary_group.rvs(N, size=n_mat, random_state=rng):
        z = np.linalg.eigvals(U)
        e = np.poly(-z)[::1]              # prod (x + z_i): coefficients high to low = e_0..e_N
        c = e * w                         # sum_j w_j e_j x^{N-j}
        r = np.roots(c)                   # roots = -z_i(t)
        offc = max(offc, np.max(np.abs(np.abs(r) - 1)))
        zt = -r / np.abs(r)
        for k in range(1, kmax + 1):
            v = abs(np.sum(zt ** k)) ** 2
            acc[k] += v; acc2[k] += v * v
    mean = acc / n_mat
    err = np.sqrt((acc2 / n_mat - mean ** 2) / n_mat)
    return mean, err, offc

if __name__ == "__main__":
    N = int(sys.argv[1]); n_mat = int(sys.argv[2]); kmax = int(sys.argv[3])
    for tau in [float(a) for a in sys.argv[4:]]:
        t = tau / N ** 2
        mean, err, offc = simulate(N, t, n_mat, kmax)
        print(f"N = {N}, tau = {tau}, {n_mat} matrices; largest | |root| - 1 | = {offc:.1e}")
        print("  k   u     simulated        exact     z     one-part term   share of multi-part terms")
        for k in range(1, kmax + 1):
            ex, one = exact(k, N, t)
            print(f" {k:2d} {k/N:5.3f} {mean[k]:9.4f}+-{err[k]:6.4f} {ex:9.4f} {(mean[k]-ex)/err[k]:+5.1f} {one:9.4f}   {1-one/ex:6.3f}")
