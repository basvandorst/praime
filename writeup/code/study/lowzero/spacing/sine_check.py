"""Sine-kernel (GUE limit) reference values by Fredholm determinants: Var(s), and <r~> from the Janossy density."""
import sys, os, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/lib')
import numpy as np, fred
t = time.time()
for m in (24, 32, 48):
    s, E = fred.E0_grid(fred.sinc_kernel, smax=5.0, ds=0.005, m=m)
    print('m=%d  Var(s) = %.7f   E(0;1)=%.7f  E(0;2)=%.3e' % (m, fred.var_from_E(s, E), E[200], E[400]))
# normalisation of the Janossy density and <r~>: integrate over s1, s2 in (0, 4)
def rtilde(Kfun, m=28, n=48, smax=4.2):
    # variables: total T in (0, smax), fraction f = s1/T in (0, 1); jac = T
    xT, wT = fred.gl(n, 0, smax); xf, wf = fred.gl(n, 0, 1)
    norm = 0.0; mr = 0.0; m1 = 0.0
    for T, wt in zip(xT, wT):
        for f, wfv in zip(xf, wf):
            s1, s2 = f * T, (1 - f) * T
            j = fred.janossy3(Kfun, s1, s2, m)
            ww = wt * wfv * T * j
            norm += ww; mr += ww * min(s1, s2) / max(s1, s2); m1 += ww * s1
    return norm, mr / norm, m1 / norm
for m, n in ((20, 32), (28, 48)):
    nrm, r, m1 = rtilde(fred.sinc_kernel, m, n)
    print('Janossy m=%d n=%d: norm %.6f  <s1> %.6f  <r~> %.6f' % (m, n, nrm, m1, r))
print('time %.1f s' % (time.time() - t))
