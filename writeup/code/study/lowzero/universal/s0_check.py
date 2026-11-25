"""build every L-function, print parameters, root number, FE check, the phase identity at 3 points, lfuneuler timing."""
import sys, os, time, math, cmath, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, cypari2
from ulib import *
pari = cypari2.Pari(); pari.allocatemem(900 * 10 ** 6); pari.default('realprecision', 38)
out = {}
for key, spec, T2, fam in CATALOGUE:
    t0 = time.time()
    try:
        L = build(pari, spec)
        Q, k, mu = params(pari, L)
        rr = pari.lfunrootres(L)
        eps = to_c(rr[2])
        feq = float(pari.lfuncheckfeq(L))
        errs = []
        for t in (7.3, 23.1, 51.7):
            v = to_c(pari.lfun(L, k / 2 + 1j * t))
            th = float(theta(t, Q, mu))
            errs.append(abs(v - eps * cmath.exp(-2j * th) * v.conjugate()) / max(abs(v), 1e-30))
        t1 = time.time(); loc = Local(pari, L, spec, Q, k, 3000); P2 = loc.poly(2)
        for p in (3, 5, 7, 11, 13, 101, 1009): loc.poly(p)
        tl = (time.time() - t1)
        out[key] = dict(Q=Q, k=k, mu=mu, eps=[eps.real, eps.imag], feq=feq, phase_err=max(errs))
        print('%-5s Q=%-7g k=%-4g mu=%-28s eps=%-22s feq=%-6.0f phase_err=%.1e  euler %.3fs  P2=%s  (%.1fs)' % (
            key, Q, k, [round(m, 2) for m in mu], '%.6f%+.6fi' % (eps.real, eps.imag), feq, max(errs), tl,
            [complex(round(c.real, 4), round(c.imag, 4)) for c in P2], time.time() - t0), flush=True)
    except Exception as e:
        print(key, 'FAILED', repr(e)[:300], flush=True)
json.dump(out, open('s0_check.json', 'w'), indent=1)
