"""Per-prime changes of Var(s) and <r~> (p missing minus p present), the wave closure on the sine base, exact in the
wave; against the second-order formula with Berry's structure function and against the slow-wave (dilation) limit."""
import sys, os, math, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/lib')
import numpy as np
import closure, gapcalc
t0 = time.time()
B = gapcalc.Base(); print('tables ready %.0f s' % (time.time() - t0))
m1, v0, r0, nrm = gapcalc.stats_with(B, None)
print('sine base from the tables: <s> %.6f Var %.6f <r~> %.6f norm %.6f  (exact: 1, 0.179998, 0.599750)' % (m1, v0, r0, nrm))

def berry_v(p, L):
    lp = math.log(p)
    def v(u):
        u = np.asarray(u, float); s = 0
        m = 1
        while p ** m < 1e16:
            s = s + np.sin(math.pi * m * u * lp / L) ** 2 / (m * m * p ** m); m += 1
        return (2 / math.pi ** 2) * s
    return v

out = {'base': {'mean': m1, 'var': v0, 'r': r0}, 'cells': []}
PR = (2, 3, 5, 7, 11, 13, 17, 19, 23)
print('\n  L     p | dVar: closure  2nd-order  slow(dilation)  BBLM 1/L^2  BBLM+Q | d<r~>: closure  2nd-order   (missing minus present)')
for L in (8.0, 9.4, 10.2, 10.7, 11.0):
    for p in PR:
        w = closure.Wave(p, L)
        mw, vw, rw, nw = gapcalc.stats_with(B, w)
        dv2, dr2 = gapcalc.second_order(B, berry_v(p, L))
        # slow limit: dilation by the local density rho(0) (exact distribution of the wave's density): <s^2> = <u^2> E[1/rho]
        ph = 2 * math.pi * np.arange(4096) / 4096; rho = w.rho(0.0, ph)
        vslow = (1 + v0) * np.mean(1 / rho) - 1
        lam = math.log(p) ** 2 / (p - 1); q = math.log(p) ** 3 / (p - 1) ** 2
        bblm = -(2 * lam / L ** 2) * (1 + v0); bblmq = bblm - 6 * q * (1 + v0) / L ** 3
        cell = dict(L=L, p=p, dvar_closure=v0 - vw, dvar_2nd=-dv2, dvar_slow=v0 - vslow, dvar_bblm=bblm, dvar_bblmq=bblmq,
                    dr_closure=r0 - rw, dr_2nd=-dr2, mean_with=mw, norm_with=nw)
        out['cells'].append(cell)
        print('%5.1f %3d | %+.5f  %+.5f  %+.5f  %+.5f  %+.5f | %+.5f  %+.5f   (<s> %.5f norm %.5f)' % (L, p, cell['dvar_closure'], cell['dvar_2nd'], cell['dvar_slow'], bblm, bblmq, cell['dr_closure'], cell['dr_2nd'], mw, nw))
json.dump(out, open('percell.json', 'w'), indent=1)
print('time %.0f s' % (time.time() - t0))
