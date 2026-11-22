"""Independent re-implementation (reviewer) of the conditioned gap law, y = 13.
ln S_{>13}(set) for a set of even offsets: exact local factors for 13 < q <= g/2 (only there can two even offsets
in [0, g] collide mod an odd q), generic factor (1-k/q)/(1-1/q)^k for q > g/2 summed with primes up to 2e7."""
import math, sys, json, time
import numpy as np

def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.nonzero(s)[0]

PR = primes_upto(20_000_000).astype(np.float64)
Y = 13
SMALL = [2, 3, 5, 7, 11, 13]
PP = 30030
PHI = 5760
RHO_DEN = np.prod([1 - 1 / q for q in SMALL])


class Gap:
    def __init__(self, g):
        self.g = g
        self.pts = np.arange(2, g, 2)
        self.n = n = len(self.pts)
        mid = [int(q) for q in PR if 13 < q <= g // 2]
        big = PR[PR > g // 2]
        self.tail = {k: float(np.sum(np.log1p(-k / big) - k * np.log1p(-1 / big))) for k in (2, 3, 4, 5)}
        self.mid = mid
        # residues: r0 = 0, rg = g mod q, ri = pts mod q
        self.res = [(q, g % q, self.pts % q) for q in mid]
        # ln S2
        lS2 = self.tail[2]
        for q, rg, _ in self.res:
            nu = 1 if rg == 0 else 2
            lS2 += math.log((1 - nu / q) / (1 - 1 / q) ** 2)
        self.lS2 = lS2
        # ln S3(i) for {0,g,i}
        lS3 = np.full(n, self.tail[3])
        for q, rg, ri in self.res:
            nu = (1 if rg == 0 else 2) + ((ri != 0) & (ri != rg)).astype(int)
            lS3 += np.log((1 - nu / q) / (1 - 1 / q) ** 3)
        self.lS3 = lS3
        self.s1 = np.exp(lS3 - lS2)
        # pairs
        lS4 = np.full((n, n), self.tail[4])
        for q, rg, ri in self.res:
            base = (1 if rg == 0 else 2)
            newi = ((ri != 0) & (ri != rg)).astype(int)
            newj_given_i = ((ri[None, :] != 0) & (ri[None, :] != rg) & (ri[None, :] != ri[:, None])).astype(int)
            nu = base + newi[:, None] + newj_given_i
            lS4 += np.log((1 - nu / q) / (1 - 1 / q) ** 4)
        self.r = np.exp(lS4 - lS3[:, None] - lS3[None, :] + lS2) - 1
        np.fill_diagonal(self.r, 0.0)
        a = np.array([x for x in range(PP) if math.gcd(x, PP) == 1])
        a = a[np.gcd(a + g, PP) == 1]
        self.M = (np.gcd(a[:, None] + self.pts[None, :], PP) == 1).astype(np.float64)

    def kappa_i(self, i):
        """kappa_ijk for j,k > i (n-i-1 square, zero diag)."""
        n = self.n; idx = np.arange(i + 1, n)
        m = len(idx)
        lS5 = np.full((m, m), self.tail[5])
        for q, rg, ri in self.res:
            base = (1 if rg == 0 else 2)
            x = ri[i]
            ni = int(x != 0 and x != rg)
            rj = ri[idx]
            nj = ((rj != 0) & (rj != rg) & (rj != x)).astype(int)
            nk = ((rj[None, :] != 0) & (rj[None, :] != rg) & (rj[None, :] != x) & (rj[None, :] != rj[:, None])).astype(int)
            nu = base + ni + nj[:, None] + nk
            lS5 += np.log(np.maximum(1 - nu / q, 1e-300) / (1 - 1 / q) ** 5)
        s3 = np.exp(lS5 - self.lS2)
        s1 = self.s1
        K = s3 / (s1[i] * s1[idx][:, None] * s1[idx][None, :]) - 1 - self.r[i, idx][:, None] - self.r[i, idx][None, :] - self.r[np.ix_(idx, idx)]
        np.fill_diagonal(K, 0.0)
        return K

    def lnshare2(self, lam):
        rho0 = 1 / (lam * RHO_DEN)
        p = rho0 * self.s1; w = p / (1 - p)
        C = w[:, None] * w[None, :] * self.r
        L = self.M @ np.log1p(-p) + 0.5 * np.einsum('ai,ai->a', self.M @ C, self.M)
        return math.log(rho0) + self.lS2 + math.log(np.sum(np.exp(L - L.max())) / PHI) + L.max()

    def delta3(self, lams):
        """ln(share order3+overlap / share order2) at each lam."""
        M = self.M; n = self.n
        out = []
        pre = []
        for lam in lams:
            rho0 = 1 / (lam * RHO_DEN)
            p = rho0 * self.s1; w = p / (1 - p)
            C = w[:, None] * w[None, :] * self.r
            L2 = M @ np.log1p(-p) + 0.5 * np.einsum('ai,ai->a', M @ C, M)
            # overlap: -1/2 sum_{i<j} C^2 - sum_{paths i-j-k, i<k} C_ij C_jk  (within R_a)
            v = M @ C; v2 = M @ (C * C)
            ov = -0.25 * np.einsum('aj,aj->a', v2, M) - 0.5 * np.einsum('aj,aj->a', v * v - v2, M)
            pre.append((L2, ov, M * w[None, :]))
        o3 = [np.zeros(M.shape[0]) for _ in lams]
        for i in range(n - 2):
            K = self.kappa_i(i)
            for t, (L2, ov, W) in enumerate(pre):
                Wi = W[:, i + 1:]
                o3[t] += W[:, i] * 0.5 * np.einsum('aj,aj->a', Wi @ K, Wi)
        for t, (L2, ov, W) in enumerate(pre):
            m = L2.max()
            out.append(math.log(np.sum(np.exp(L2 - o3[t] + ov - m)) / np.sum(np.exp(L2 - m))))
        return out
