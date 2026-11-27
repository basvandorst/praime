# Dual certificate (float64 check; the ball-arithmetic version is a separate step).
# From the LP multipliers y_i >= 0 of the cuts: P_H = sum y_i v_i v_i' (cuts of H), P_G = sum y_i u_i u_i' (cuts of G_T).
# They define ONE test function g with ghat(gamma) = sum |phihat_i|^2 + (gamma^2 - T^2) sum |psihat_j|^2, which is >= 0
# outside (-T, T).  For an L-function of the class with no zero in (-T, T):  0 <= sum_gamma ghat(gamma) = W(g), and
#   W(g)/d = sum_kappa nu(kappa) F(kappa) + l * c0 + sum_p <P_p(theta)>,
# with F(kappa) = Phi(A_kappa), c0 = Phi(2a I), P_p(theta) = -sum_m cos(m theta) Phi(S_pm).  So
#   W(g)/d <= sup_kappa F(kappa) + max(c0, 0)*l + sum_p max_theta P_p(theta).
# If c0 <= 0 and  B = sup_kappa F + sum_p max_theta P_p < 0  the class has no member with a zero-free (-T, T).
# Here the suprema are taken on fine grids (kappa far beyond the LP grid, 20001 angles).
# Usage: cert.py <real|complex> <X> <T> [out.npz]
import sys, numpy as np
from gap import Gap, shift_grid, EULER_LOGPI
def build(cls, X, T):
    g = Gap(X, shift_grid(cls)); lam, t, rd = g.solve(T)
    y = g.dual; n = g.n
    PH = np.zeros((n, n)); PG = np.zeros((n - 1, n - 1))
    for (kind, v), yi in zip(g.cuts, y):
        if yi <= 0: continue
        if kind == 'H': PH += yi*np.outer(v, v)
        else: PG += yi*np.outer(v, v)
    return g, t, PH, PG
def Phi(g, T, PH, PG):
    W = g.Z @ PG @ g.Z.T; KWK = g.Kd[:, None]*W*g.Kd[None, :]
    D = PH + KWK - T*T*W                                    # Phi(M) = <D, M>
    return lambda M: np.einsum('jk,...jk->...', D, M)
def check(g, T, PH, PG, cls, verbose=True):
    phi = Phi(g, T, PH, PG); a = g.a; B = g.B
    c0 = float(phi(2*a*g.I))
    # F(kappa) on a fine grid, from the scalar s(u) = Phi(S(u)) at the quadrature nodes
    s = np.concatenate([phi(B.S(B.u[i:i + 400])) for i in range(0, len(B.u), 400)])
    den = 1 - np.exp(-2*B.u)
    const = c0*(2*np.sum(B.wt*np.exp(-2*B.u)/den) - EULER_LOGPI - np.log(1 - np.exp(-4*a)))
    def F(kap):
        kap = np.atleast_1d(np.asarray(kap, complex))
        out = np.empty(len(kap)); wsd = B.wt*s/den
        for i in range(0, len(kap), 2000):
            k = kap[i:i + 2000]
            out[i:i + 2000] = const - (np.exp(-(0.5 + k.real)[:, None]*B.u[None, :])*np.cos(k.imag[:, None]*B.u[None, :])) @ wsd
        return out
    if cls == 'real':
        ks = np.arange(0, 200.0001, 0.005); Fv = F(ks); i = int(np.argmax(Fv)); supF = float(Fv[i]); arg = '%g' % ks[i]
        tail = float(F([1e3, 1e4, 1e6]).max())
    else:
        best = (-np.inf, None)
        for re in np.arange(0, 12.01, 0.05):
            im = np.arange(0, 80.0001, 0.02); Fv = F(re + 1j*im); i = int(np.argmax(Fv))
            if Fv[i] > best[0]: best = (float(Fv[i]), '%g%+gi' % (re, im[i]))
        supF, arg = best; tail = float(max(F(np.array([50, 200, 1e4]) + 0j).max(), F(1j*np.array([200.0, 1e3, 1e4])).max()))
    th = np.linspace(0, np.pi, 20001); tot = supF; per = []
    for ip, p in enumerate(g.primes):
        idx = [i for i, (q, m) in enumerate(g.pp) if q == p]
        Pp = -sum(np.cos(g.pp[i][1]*th)*float(phi(g.S[i])) for i in idx); per.append((p, float(Pp.max()), float(th[int(np.argmax(Pp))])))
        tot += Pp.max()
    if verbose:
        print('certificate %s X %g T %g: rank P_H %d, rank P_G %d; g(0) coefficient c0 = %.3e (must be <= 0)' % (cls, g.X, T, np.linalg.matrix_rank(PH, 1e-12), np.linalg.matrix_rank(PG, 1e-12), c0))
        print('  sup over shifts F = %.6e at kappa = %s; far shifts: %.3e' % (supF, arg, tail))
        print('  primes: ' + ' '.join('%d: max %.4e at theta %.4f (cos %.4f)' % (p, v, t_, np.cos(t_)) for p, v, t_ in per))
        print('  B = sup F + sum of the prime maxima = %.6e  (certificate holds if B < 0 and c0 <= 0)' % tot, flush=True)
    return c0, supF, per, tot
if __name__ == '__main__':
    cls = sys.argv[1]; X = float(sys.argv[2]); T = float(sys.argv[3])
    g, t, PH, PG = build(cls, X, T)
    print('LP value t* = %.6e (normalisation: sum of multipliers times scale = %.6f)' % (t, sum(yi*(2*g.a if k == 'H' else 2*g.a*T*T) for (k, v), yi in zip(g.cuts, g.dual))))
    c0, supF, per, tot = check(g, T, PH, PG, cls)
    if len(sys.argv) > 4: np.savez(sys.argv[4], PH=PH, PG=PG, X=X, T=T, N=g.N, c0=c0, supF=supF, B=tot)
