# Truncated Weil forms for L-functions with complex Gamma shifts (float64, own code; independent of ../weil.py).
# Test functions f on [-a, a], basis exp(i pi k u / a), |k| <= N; X = e^(2a) is the support.
# For a Gamma factor prod_i Gamma_R(s + kappa_i) (the multiset of kappa closed under conjugation), conductor 1,
# pole of order `pole` at s = 1 and Lambda-coefficients c(n) (c(p^m) = weight * log p):
#   Q(f) = pole * 2 int_{-2a}^{2a} g cosh(w/2) dw + sum_i A_kappa_i(g) - 2 sum_n c(n) n^(-1/2) g(log n),
#   A_kappa(g) = -(gamma + log pi) g(0) + int_0^inf [2 e^(-2u) g(0) - 2 Re e^(-(1/2+kappa)u) g(u)] / (1 - e^(-2u)) du,
# with g the autocorrelation of f.  For a true L-function Q(f) = sum_gamma |fhat(gamma)|^2 (>= 0 under its GRH).
import numpy as np
EULER = 0.5772156649015329

def S_stack(a, N, w):
    """S_jk(w) = g_jk(w) + g_jk(-w) for an array of w in [0, 2a]; shape (len(w), n, n)."""
    w = np.atleast_1d(np.asarray(w, float)); ks = np.arange(-N, N + 1)
    sn = np.sin(np.pi*np.outer(w, ks)/a); cs = np.cos(np.pi*np.outer(w, ks)/a)
    d = ks[:, None] - ks[None, :]
    sign = np.where(d % 2 == 0, 1.0, -1.0)
    with np.errstate(divide='ignore', invalid='ignore'):
        fac = np.where(d != 0, sign*2*a/(np.pi*np.where(d == 0, 1, d)), 0.0)
    S = fac[None, :, :]*(sn[:, None, :] - sn[:, :, None])          # (sin(pi k w/a) - sin(pi j w/a)), j row, k column
    idx = np.arange(len(ks))
    S[:, idx, idx] = 2*(2*a - w)[:, None]*cs
    return S

def gl(lo, hi, panels, deg=8):
    x, wt = np.polynomial.legendre.leggauss(deg)
    edges = np.linspace(lo, hi, panels + 1); h = (hi - lo)/panels
    nodes = (edges[:-1, None] + (x[None, :] + 1)*h/2).ravel()
    return nodes, np.tile(wt*h/2, panels)

class Basis:
    def __init__(self, X, N, maxfreq=0.0):
        self.X = X; self.a = a = np.log(X)/2; self.N = N; self.n = 2*N + 1
        panels = int(40 + 2*a*(2*np.pi*N/a + maxfreq)/np.pi*1.5)
        self.u, self.wt = gl(0.0, 2*a, panels)
        self._S = None
    def S(self, w):
        return S_stack(self.a, self.N, w)
    def integrate(self, fun_scalar, chunk=400):
        """int_0^2a fun(u) * S(u) du as a matrix; fun_scalar gives an array over nodes."""
        vals = fun_scalar(self.u)*self.wt; out = np.zeros((self.n, self.n))
        for i in range(0, len(self.u), chunk):
            out += np.einsum('i,ijk->jk', vals[i:i+chunk], self.S(self.u[i:i+chunk]))
        return out
    def arch(self, kappas, pole=0):
        """Matrix of pole + archimedean part, conductor 1."""
        a = self.a; n = self.n; ng = len(kappas); I = np.eye(n)
        K = lambda u: sum(np.exp(-(0.5 + k.real)*u)*np.cos(k.imag*u) for k in map(complex, kappas))
        # int_0^2a [2 ng e^(-2u) g(0) - 2 K(u) g(u)]/(1-e^(-2u)) du, with g -> S/2 and g(0) -> 2a I
        A = -self.integrate(lambda u: K(u)/(1 - np.exp(-2*u)))
        # the g(0) part, with the K-part's singularity at u = 0 cancelling: add and subtract by quadrature of scalars
        A += I*2*a*2*ng*np.sum(self.wt*np.exp(-2*self.u)/(1 - np.exp(-2*self.u)))
        A += I*2*a*ng*(-(EULER + np.log(np.pi)) - np.log(1 - np.exp(-4*a)))
        if pole:
            A += pole*2*self.integrate(lambda u: np.cosh(u/2))
        return A
