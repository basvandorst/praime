"""The within part of block transition tables from the Hardy-Littlewood conjecture, finite-lambda evaluation.
Own code (no repo imports).  Two pieces:
 (a) the fine-class chain: the next-prime law per rough class u mod Q (order-1 conditioned HL expansion, primes <= y
     exact by CRT, primes > y through the HL pair factor S_{>y}(t) with the predecessor); within-projected marks of the
     mod-q transitions; their lag covariances tau_j computed exactly by powers of the kernel (no simulation);
 (b) the long-range tail: the Hardy-Littlewood pair correlations of the primes > y between non-adjacent primes, acting
     on the insertion responses r_u of the within-projected table (linear response):
       W_tail(H) = (2/((C-1)^2 lam)) (Q/phi(Q)) sum_{h<H} (1 - h/H) w2(h) A_W(h mod Q) (S_{>y}(h) - 1),
     A_W(h) = mean_u <Pi_W r_u, Pi_W r_{u+h}>,  w2(h) = 1 - exp(-h/lam) (share of pairs at distance h that are not
     adjacent primes; heuristic).
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, math
from scipy import sparse
from hlcore import primorial, rough_classes, S_tail_table, primes_upto

def proj_within(Pm, pi):
    C = len(pi); R = np.zeros((C * C, C * C)); V = np.zeros((C * C, C))
    for a in range(C):
        sp = np.sqrt(Pm[a]); R[a * C:(a + 1) * C, a * C:(a + 1) * C] = np.eye(C) - np.outer(sp, sp)
        for b in range(C): V[a * C + b, b] = math.sqrt(pi[a] * Pm[a, b])
    return R - R @ V @ np.linalg.pinv(V.T @ R @ V) @ V.T @ R

def increment_vectors(Pm, pi):
    """Yv[a*C+b] = y(a,b) = (e_ab - sum_b' P_ab' e_ab')/sqrt(pi_a P_ab .) in normalized coordinates"""
    C = len(pi); Yv = np.zeros((C * C, C * C))
    for a in range(C):
        for b in range(C):
            for b2 in range(C): Yv[a * C + b, a * C + b2] -= Pm[a, b2] / math.sqrt(pi[a] * Pm[a, b2])
            Yv[a * C + b, a * C + b] += 1 / math.sqrt(pi[a] * Pm[a, b])
    return Yv

class Model:
    def __init__(self, lam, y=13, T=800, Smax=1_100_000):
        self.lam = lam; self.y = y
        Q = primorial(y); U = rough_classes(Q); n = len(U); rho0 = (Q / n) / lam
        self.Q, self.U, self.n, self.rho0 = Q, U, n, rho0
        St = S_tail_table(T, y)[1:]; t = np.arange(1, T + 1)
        isr = np.zeros(Q, bool); isr[U] = True
        index = -np.ones(Q, np.int64); index[U] = np.arange(n); self.index = index
        rows, cols, vals = [], [], []
        WL = np.zeros((n, T)); WR = np.zeros((n, T))
        for side, W in ((1, WR), (-1, WL)):
            pos = (U[:, None] + side * t[None, :]) % Q
            pt = np.where(isr[pos], rho0 * St[None, :], 0.0)
            surv = np.cumprod(1 - pt, axis=1); surv = np.concatenate([np.ones((n, 1)), surv[:, :-1]], axis=1)
            W[:] = pt * surv
            if side == 1:
                nz = W > 1e-15
                ii, jj = np.nonzero(nz)
                rows = ii; cols = index[pos[ii, jj]]; vals = W[ii, jj]
        self.K = sparse.csr_matrix((vals, (rows, cols)), shape=(n, n))
        self.WL, self.WR, self.t = WL, WR, t
        self.Stab = S_tail_table(Smax, y)

    def classes(self, q):
        red = [a for a in range(q) if math.gcd(a, q) == 1]; lut = -np.ones(q, np.int64); lut[red] = np.arange(len(red))
        return red, lut

    def neighbour_laws(self, q):
        red, lut = self.classes(q); C = len(red); n = self.n; T = len(self.t)
        PR = np.zeros((n, C)); PL = np.zeros((n, C))
        for side, W, P in ((1, self.WR, PR), (-1, self.WL, PL)):
            cls = lut[(self.U[:, None] + side * self.t[None, :]) % q]
            for b in range(C): P[:, b] = (W * (cls == b)).sum(1)
        cu = lut[self.U % q]
        return red, C, cu, PR, PL

    def setup_q(self, q):
        red, C, cu, PR, PL = self.neighbour_laws(q)
        Pm = np.array([PR[cu == a].mean(0) for a in range(C)]); Pm /= Pm.sum(1, keepdims=True)
        pi = np.full(C, 1.0 / C)
        PiW = proj_within(Pm, pi); w, Ub = np.linalg.eigh(PiW); Ub = Ub[:, w > 0.5]
        Yv = increment_vectors(Pm, pi)
        Fv = (Yv @ Ub).reshape(C, C, -1)           # within coordinates of y(a,b)
        return dict(red=red, C=C, cu=cu, PR=PR, PL=PL, Pm=Pm, Fv=Fv, D=Fv.shape[2])

    def chain_tau(self, q, J=60):
        """tau_j / (C-1)^2, j = 0..J, of the fine-class chain (k = 1 tables), exactly by kernel powers."""
        s = self.setup_q(q); C, cu, PR, Fv, D = s['C'], s['cu'], s['PR'], s['Fv'], s['D']
        K = self.K
        G = np.einsum('ub,ubd->ud', PR, Fv[cu])     # E[F(u, u') | u]
        tau = [float(np.einsum('ub,ubd->', PR, Fv[cu] ** 2) / self.n)]
        Phi = G.copy()
        # A_j(u) = sum_{u'} K(u,u') F(c_u, c_u') . Phi(u')  = sum_b Fv[c_u, b] . (K (1[c=b] Phi))(u)
        for j in range(1, J + 1):
            acc = 0.0
            for b in range(C):
                Mb = K @ (Phi * (cu == b)[:, None])
                acc += np.einsum('ud,ud->', Fv[cu, b], Mb)
            tau.append(float(acc / self.n))
            Phi = K @ Phi
        return np.array(tau) / (C - 1) ** 2, s

    def responses(self, q, s=None):
        """within coordinates of the insertion response r_u = E[y(c_L, c_u) + y(c_u, c_R) - y(c_L, c_R) | x = u]"""
        if s is None: s = self.setup_q(q)
        C, cu, PR, PL, Fv = s['C'], s['cu'], s['PR'], s['PL'], s['Fv']
        r = np.einsum('ua,uad->ud', PL, Fv[:, cu].transpose(1, 0, 2)) + np.einsum('ub,ubd->ud', PR, Fv[cu]) \
            - np.einsum('ua,ub,abd->ud', PL, PR, Fv)
        return r

    def tail(self, q, Hs, s=None, lagbins=None):
        r = self.responses(q, s); C = len([a for a in range(q) if math.gcd(a, q) == 1])
        Q, n = self.Q, self.n
        full = np.zeros((Q, r.shape[1])); full[self.U] = r
        F = np.fft.rfft(full, axis=0)
        A = np.fft.irfft((np.conj(F) * F).sum(1), n=Q) / n       # A_W(h) = mean_u <w_u, w_{u+h}>, h mod Q
        V = float(A[0])
        lam = self.lam; pref = (Q / n) / ((C - 1) ** 2 * lam)
        out = {'V': V, 'V_per_dir': V / (C - 1) ** 2}
        Hmax = max(Hs); h = np.arange(1, Hmax + 1)
        base = A[h % Q] * (self.Stab[1:Hmax + 1] - 1) * (1 - np.exp(-h / lam))
        out['W_tail'] = {H: float(2 * pref * np.sum((1 - h[:H] / H) * base[:H])) for H in Hs}
        if lagbins:
            # index-lag view: pairs at distance h have lag j with Poisson(h/lam) weights (j - 1 primes between)
            res = {}
            hh = np.arange(1, 200 * int(lam) * 1 + 1)
            hh = np.arange(1, min(Hmax, 4000 * int(lam)) + 1)
            b0 = A[hh % Q] * (self.Stab[hh] - 1)
            mu = hh / lam
            from scipy.stats import poisson
            for lo, hi in lagbins:
                wgt = poisson.cdf(hi - 2, mu) - poisson.cdf(lo - 2, mu)     # P(lo <= lag < hi), lag - 1 ~ Poisson(mu)
                res[(lo, hi)] = float(pref * np.sum(wgt * b0))
            out['bins'] = res
        return out

# ---------------- additions: sieve-only option, k = 2 chain covariances, k = 2 responses ----------------
class Model2(Model):
    """Model with an optional sieve-only kernel (S_{>y} == 1: the sieve chain M_Q) and a left kernel."""
    def __init__(self, lam, y=13, T=800, Smax=1_100_000, sieve_only=False):
        if not sieve_only:
            super().__init__(lam, y, T, Smax)
        else:
            import hlcore
            orig = hlcore.S_tail_table
            import chainlaw as _cl
            _cl.S_tail_table = lambda hmax, yy: np.where(np.arange(hmax + 1) == 0, np.nan, 1.0)
            try:
                super().__init__(lam, y, T, Smax=10)
            finally:
                _cl.S_tail_table = orig
        n = self.n; Q = self.Q; t = self.t
        pos = (self.U[:, None] - t[None, :]) % Q
        nz = self.WL > 1e-15; ii, jj = np.nonzero(nz)
        self.KL = sparse.csr_matrix((self.WL[ii, jj], (ii, self.index[pos[ii, jj]])), shape=(n, n))

    def chain_tau_k2(self, q, J=60):
        """tau_j/(C-1)^2 for the k = 2 table (pairs c_n, c_{n+2}) in the fine-class chain, exactly."""
        red, C, cu, PR, PL = self.neighbour_laws(q)
        K = self.K; n = self.n
        P2u = K @ PR                                   # P(c_{n+2} = b | u_n)
        Pm = np.array([P2u[cu == a].mean(0) for a in range(C)]); Pm /= Pm.sum(1, keepdims=True); pi = np.full(C, 1 / C)
        PiW = proj_within(Pm, pi); w, Ub = np.linalg.eigh(PiW); Ub = Ub[:, w > 0.5]
        Fv = (increment_vectors(Pm, pi) @ Ub).reshape(C, C, -1); D = Fv.shape[2]
        G2 = np.einsum('ub,ubd->ud', P2u, Fv[cu])      # E[F2(u_n, u_{n+2}) | u_n]
        pi_u = np.full(n, 1.0 / n)
        tau = [float(np.einsum('ub,ubd->', P2u, Fv[cu] ** 2) / n)]
        # lag 1: E[F2(c0,c2) F2(c1,c3)] = sum_{a0,a1} sum_{u2} w_{a0 a1}(u2) <Fv[a0, c(u2)], Hm(a1, u2)>,
        #   Hm(a1, u2) = sum_b PR(u2, b) Fv[a1, b],  w_{a0 a1} = K'(1[a1] * K'(pi 1[a0]))
        KT = K.T.tocsr()
        Hm = np.einsum('ub,abd->aud', PR, Fv)          # (C, n, D)
        acc = 0.0
        for a0 in range(C):
            v0 = KT @ (pi_u * (cu == a0))
            for a1 in range(C):
                w01 = KT @ (v0 * (cu == a1))
                acc += np.einsum('u,ud,ud->', w01, Fv[a0, cu], Hm[a1])
        tau.append(acc)
        # lag j >= 2: E[F2(u0,u2) . (K^{j-2} G2)(u2)]
        Phi = G2.copy()
        K2 = K @ K
        for j in range(2, J + 1):
            acc = 0.0
            for b in range(C):
                Mb = K2 @ (Phi * (cu == b)[:, None])
                acc += np.einsum('ud,ud->', Fv[cu, b], Mb)
            tau.append(float(acc / n))
            Phi = K @ Phi
        return np.array(tau) / (C - 1) ** 2, dict(red=red, C=C, cu=cu, PR=PR, PL=PL, P2u=P2u, Pm=Pm, Fv=Fv)

    def responses_k2(self, q, s):
        C, cu, PR, PL, Fv = s['C'], s['cu'], s['PR'], s['PL'], s['Fv']
        PR2 = self.K @ PR; PL2 = self.KL @ PL
        r = (np.einsum('ua,uad->ud', PL2, Fv[:, cu].transpose(1, 0, 2))
             + np.einsum('ua,ub,abd->ud', PL, PR, Fv)
             + np.einsum('ub,ubd->ud', PR2, Fv[cu])
             - np.einsum('ua,ub,abd->ud', PL2, PR, Fv)
             - np.einsum('ua,ub,abd->ud', PL, PR2, Fv))
        return r
