"""Loop 2b: the limit law against the exact finite-N value, computed without root finding.
For k <= N the traces p_i = sqrt(i) g_i (independent complex Gaussians) have the exact joint moments (Diaconis-Shahshahani),
so E|p_k(t)|^2 is the Gaussian expectation of |k [x^k] log sum_j e^{t j^2} e_j x^j|^2 e^{-2tNk}; power series by Newton's identities."""
import sys, numpy as np
from loop2 import K_limit

def gauss_K(N, tau, kmax, n_s, seed=3):
    t = tau / N ** 2
    rng = np.random.default_rng(seed)
    g = (rng.standard_normal((n_s, kmax + 1)) + 1j * rng.standard_normal((n_s, kmax + 1))) / np.sqrt(2)
    i = np.arange(kmax + 1)
    p = np.sqrt(i) * g                                   # p_0 unused
    sgn = (-1.0) ** (i - 1)
    e = np.zeros((n_s, kmax + 1), complex); e[:, 0] = 1
    for j in range(1, kmax + 1):                         # j e_j = sum_{i=1..j} (-1)^{i-1} p_i e_{j-i}
        e[:, j] = np.sum(sgn[1:j + 1] * p[:, 1:j + 1] * e[:, j - 1::-1][:, :j], axis=1) / j
    a = e * np.exp(t * i * i)
    q = np.zeros((n_s, kmax + 1), complex)               # q_k = (-1)^{k-1} p_k(t) e^{tNk}
    for k in range(1, kmax + 1):                         # k a_k = sum_{i=1..k} q_i a_{k-i}
        q[:, k] = k * a[:, k] - np.sum(q[:, 1:k] * a[:, k - 1:0:-1], axis=1)
    v = np.abs(q) ** 2 * np.exp(-2 * t * N * i)
    return v.mean(axis=0) / N, v.std(axis=0) / np.sqrt(n_s) / N

if __name__ == "__main__":
    N = int(sys.argv[1]); n_s = int(sys.argv[2])
    for tau in [float(x) for x in sys.argv[3:]]:
        kmax = int(0.9 * N)
        m, s = gauss_K(N, tau, kmax, n_s)
        print(f"N = {N}, tau = {tau}, {n_s} Gaussian samples (exact moments for k <= N)")
        print("   u    finite N         series r<=4    z      terms r = 1, 2, 3, 4")
        for u in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9):
            k = int(round(u * N)); T = K_limit(u, tau)
            print(f" {u:5.2f} {m[k]:8.4f}+-{s[k]:6.4f}  {sum(T):8.4f}  {(m[k]-sum(T))/s[k]:+5.1f}   " + " ".join(f"{x:8.5f}" for x in T))
