"""Determinantal (kernel-change) closure: the gap ratio and the spacing variance under a kernel change
K = sinc + c * K1.  (a) K1 = (pi x/6) sin(pi x) is the 1/N^2 term of CUE(N) (BBLM eq. 6); we check that <r~> has no
first-order change (the CUE identity) and that dVar/d(1/N^2) = -(1 + V0)/6.  (b) K1 = x^2 cos(pi x) is the kernel
that carries the Q-term of the two-point function, -(4 pi Q/L^3) r sin(2 pi r) = -2 sinc * dK -> dK = (4 pi^2 Q/L^3) x^2 cos(pi x)
(the Nishigaki-type kernel correction): its d<r~>/dc and dVar/dc."""
import sys, os, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/lib')
import numpy as np, fred

def stats(Kfun):
    s, E = fred.E0_grid(Kfun, smax=5.0, ds=0.005, m=32)
    v = fred.var_from_E(s, E); r, nrm = fred.ratio_mean(Kfun)
    return v, r, nrm

v0, r0, n0 = stats(fred.sinc_kernel)
print('sine kernel: Var(s) = %.7f   <r~> = %.7f   (Janossy norm %.7f)' % (v0, r0, n0))
out = {'sine': {'var': v0, 'r': r0}}
for name, K1 in (('cue_1/N^2', lambda x: (np.pi * x / 6) * np.sin(np.pi * x)), ('Q_kernel', lambda x: x ** 2 * np.cos(np.pi * x))):
    res = []
    for c in (-0.02, -0.01, 0.01, 0.02):
        v, r, nrm = stats(lambda x, c=c: np.sinc(x) + c * K1(x))
        res.append((c, v, r, nrm))
        print('  %-10s c=%+.3f  Var %.7f  <r~> %.7f  norm %.6f' % (name, c, v, r, nrm))
    c = np.array([x[0] for x in res]); V = np.array([x[1] for x in res]); R = np.array([x[2] for x in res])
    # derivative by central differences (Richardson with h = 0.01, 0.02)
    dV = (8 * (V[2] - V[1]) - (V[3] - V[0])) / (12 * 0.01); dR = (8 * (R[2] - R[1]) - (R[3] - R[0])) / (12 * 0.01)
    print('  %-10s dVar/dc = %+.5f   d<r~>/dc = %+.5f' % (name, dV, dR))
    out[name] = {'dVar_dc': dV, 'dR_dc': dR, 'rows': res}
print('check: CUE 1/N^2 coefficient of Var predicted -(1 + V0)/6 = %.5f' % (-(1 + v0) / 6))
# CUE(N) exact (Toeplitz, BBLM eq. 22) for comparison
def cue_E(s, N):
    j = np.arange(N); d = j[:, None] - j[None, :]
    with np.errstate(invalid='ignore', divide='ignore'):
        M = np.where(d == 0, s / N, np.sin(np.pi * s * d / N) / (np.pi * d))
    return np.linalg.det(np.eye(N) - M)
print('CUE(N) exact by the Toeplitz determinant: Var(s), and N^2 (Var - V0)')
cue = {}
for N in (2, 3, 4, 5, 6, 8, 10, 12, 16, 24, 32):
    s = np.arange(0, N + 1e-9, 0.002); E = np.array([cue_E(v, N) for v in s])
    V = 2 * np.trapezoid(E, s) - 1; cue[N] = V
    print('  N=%2d  Var %.6f   N^2 (Var - V0) = %+.5f' % (N, V, N * N * (V - v0)))
out['cue_var'] = cue
json.dump(out, open('kernels.json', 'w'), indent=1)
