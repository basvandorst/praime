"""fresh.py: frozen predictions for the variance of log|L'(rho, chi_D)| - ln L at the zeros 600 <= gamma <= 1000 of the
fresh sets C (the discriminants of r134_C.npz) and Bp (those of r138_Bp.npz). Only the discriminants are read from the
zero files (the keys), never a zero height. For every L-function: the zero-weighted band average of
  cum2_3pt  the three-point (ratios-conjecture) change of the variance by its bad primes (assemble.Theory),
  cum2_ldm  the registered round-142 heuristic (leading order + LDM 1/L term),   cum2_lead  leading order only,
  cum2_none 0,   cum1_3pt  the mean shift (two-point law; secondary),
and mean_lnL (zero-weighted mean of ln L), as lprime_test.design uses it. Output JSON + its sha256.
usage: fresh.py OUT.json
"""
import sys, os, json, math, hashlib, importlib.util
sys.dont_write_bytecode = True
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from assemble import Theory
from retro import load_mod, LP, REPO
LO, HI = 600.0, 1000.0
ZDIR = REPO + '/data/riemann/lzeros/'
X_, W_ = np.polynomial.legendre.leggauss(32)


def keys_of(f):
    with np.load(f) as z: return sorted((int(k[1:]) for k in z.files), key=lambda d: (abs(d), d))


def main(out):
    th = Theory()
    eul = {}
    J = {'version': 'varlp-predictions-1', 'heights': [LO, HI], 'statistic': 'per-L-function sample variance (ddof 1) of log|L\'(rho)| - ln(q gamma/2pi) over 600 <= gamma <= 1000',
         'design': 'lprime_test.design: intercept, set dummy, mean_lnL, indicators of 2,3,5,7,11,13 (centred)',
         'theory': 'derivation.md (three-point Palm formula; ratios-conjecture T2, T3 with the bad primes deleted; base S = {} + pair corrections)'}
    for name, f in (('C', 'r134_C.npz'), ('Bp', 'r138_Bp.npz')):
        rows = []
        for D in keys_of(ZDIR + f):
            q = abs(D); a = 0 if D > 0 else 1; bad = LP.primes_of(q)
            t = LO + (HI - LO) * (X_ + 1) / 2; w = (HI - LO) / 2 * W_ * LP.density(t, q, a); n = float(w.sum()); w = w / w.sum()
            Lt = np.log(q * t / (2 * math.pi))
            parts = [th.dvar(bad, L, parts=True) for L in Lt]
            dv = float(sum(wi * p_['total'] for wi, p_ in zip(w, parts))); dm = float(sum(wi * p_['dmu'] for wi, p_ in zip(w, parts)))
            inv = float((w / Lt).sum()); lnL = float((w * np.log(Lt)).sum())
            lead = 0.0; ldm = 0.0
            for p in bad:
                if p not in eul: eul[p] = LP.euler_cumulants(p)
                K0, K1 = eul[p]; lead -= K0[1]; ldm -= K1[1] * math.log(p) * inv
            rows.append(dict(D=D, bad=bad, n_zeros_expected=n, mean_lnL=lnL, mean_invL=inv, Lmin=float(Lt.min()), Lmax=float(Lt.max()),
                             cum2_3pt=dv, cum2_3pt_numunc=th.unc(bad, math.exp(lnL)), cum2_ldm=lead + ldm, cum2_lead=lead, cum2_none=0.0, cum1_3pt=dm))
        J[name] = rows
        print('%s: %d L-functions, expected zeros %.0f, L %.2f .. %.2f' % (name, len(rows), sum(r['n_zeros_expected'] for r in rows),
                                                                         min(r['Lmin'] for r in rows), max(r['Lmax'] for r in rows)), flush=True)
    s = json.dumps(J, indent=0, sort_keys=True); open(out, 'w').write(s)
    print('sha256', hashlib.sha256(s.encode()).hexdigest())


if __name__ == '__main__':
    main(sys.argv[1])
