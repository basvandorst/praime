# The sine process (bulk limit of GUE, density 1) by Nystrom discretisation of the Fredholm determinant (Bornemann 2010).
# K(x, y) = sin(pi (x - y)) / (pi (x - y)).
#   E(0; s) = det(I - K on (0, s));  spacing density p1(s) = d^2/ds^2 E(0; s) = Janossy density of points at 0 and s
#   with none between;  joint density of two consecutive spacings p2(s1, s2) = Janossy density of points at 0, s1,
#   s1 + s2 with no other point in (0, s1 + s2):  det(I - K_J) * det[R_J(z_i, z_j)],  R_J = K (I - K_J)^-1 (resolvent).
import numpy as np
def K(x, y):
    d = np.subtract.outer(x, y); return np.sinc(d)                      # np.sinc(t) = sin(pi t)/(pi t)
def janossy(z, m=40):
    """density of points exactly at z (sorted, z[0] = 0) with no other point in (0, z[-1])."""
    L = z[-1]; x, w = np.polynomial.legendre.leggauss(m); t = (x + 1)*L/2; w = w*L/2; sw = np.sqrt(w)
    A = sw[:, None]*K(t, t)*sw[None, :]; M = np.eye(m) - A
    det = np.linalg.det(M); z = np.asarray(z, float)
    Kzt = K(z, t)*sw[None, :]
    R = K(z, z) + Kzt @ np.linalg.solve(M, Kzt.T)
    return det*np.linalg.det(R)
def p1(s, m=40): return janossy([0.0, s], m)
def p2(s1, s2, m=40): return janossy([0.0, s1, s1 + s2], m)
if __name__ == '__main__':
    # checks: normalisation, mean spacing 1, variance 0.180 (GUE: 0.1800), and the marginal of p2
    s, ws = np.polynomial.legendre.leggauss(120); s = (s + 1)*3.5; ws = ws*3.5                 # (0, 7)
    P1 = np.array([p1(v) for v in s])
    print('int p1 = %.8f, mean = %.8f, variance = %.6f' % (ws @ P1, ws @ (P1*s), ws @ (P1*s*s) - (ws @ (P1*s))**2))
    h1 = -ws @ (P1*np.log(np.maximum(P1, 1e-300)))
    print('entropy of one spacing h1 = %.6f nats = %.6f bits (exponential: 1 nat)' % (h1, h1/np.log(2)))
    g, wg = np.polynomial.legendre.leggauss(64); g = (g + 1)*2.25; wg = wg*2.25                # (0, 4.5)^2
    P2 = np.array([[p2(a, b) for b in g] for a in g]); P1g = np.array([p1(v) for v in g])
    print('int p2 = %.8f; marginal error max %.2e' % (wg @ P2 @ wg, np.abs(P2 @ wg - P1g).max()))
    mean = wg @ (P1g*g); var = wg @ (P1g*g*g) - mean**2
    cov = (wg*g) @ P2 @ (wg*g) - mean**2
    mi = np.einsum('i,j,ij->', wg, wg, P2*np.log(np.maximum(P2, 1e-300)/np.maximum(np.outer(P1g, P1g), 1e-300)))
    print('correlation of consecutive spacings = %.5f; mutual information = %.6f nats = %.6f bits' % (cov/var, mi, mi/np.log(2)))
