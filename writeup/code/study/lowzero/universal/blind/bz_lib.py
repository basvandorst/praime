"""blind test: L-function builders and the registered zero routine (completeness rule).
Spec kinds: ['char', q, m] | ['mf', k, N, idx] | ['sym', ainv, m] | ['g2', f, h] | ['nf', poly] | ['artin', poly] | ['zeta'].
Completeness rule (registered): lfunzeros at default subdivision on (0.01, T2+1]; then up to 3 rounds re-searching the 30
largest unfolded gaps (divz 64; 256 for gaps > 2.2 mean spacings) and the first interval; duplicates (< 1e-7) removed;
m0 = lfunorderzero at the centre; r_n = n - theta(g_n)/pi; COMPLETE iff |drift| <= 0.2 (mean r over the last 15 % minus
the first 15 %) and |mean r - c_exp| <= 0.25, c_exp = 1/2 + (poles at s = 1) - m0/2 - arg(eps)/2pi (arg taken 0 for real
eps); otherwise one more pass over every gap > 1.8 spacings with divz 256, then re-check; if still failing: REJECTED."""
import math, time, json, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from ulib import theta, to_c

def build(pari, spec):
    k = spec[0]
    if k == 'char':
        return pari('lfuncreate([G=znstar(%d,1), znconreychar(G,%d)])' % (spec[1], spec[2]))
    if k == 'mf':
        mf = pari.mfinit([spec[2], spec[1]], 0); return pari.lfunmf(mf, pari.mfeigenbasis(mf)[spec[3]])
    if k == 'sym':
        return pari.lfunsympow(pari.ellinit(spec[1]), spec[2])
    if k == 'g2':
        return pari('lfungenus2([%s,%s])' % (spec[1], spec[2]))
    if k == 'nf':
        return pari('lfuncreate(%s)' % spec[1])
    if k == 'artin':
        return pari('lfundiv(lfuncreate(%s), lfuncreate(1))' % spec[1])
    if k == 'zeta':
        return pari.lfuncreate(1)
    raise ValueError(spec)

def params(pari, L):
    N, k, V = pari.lfunparams(L); k = float(k)
    return float(N), k, [float(v) + (k - 1) / 2 for v in V]

def unfold(g, Q, mu): return theta(g, Q, mu) / math.pi

def check(g, Q, mu, cexp):
    r = np.arange(1, len(g) + 1) - unfold(g, Q, mu); m = max(5, int(0.15 * len(g)))
    drift = float(r[-m:].mean() - r[:m].mean()); rmean = float(r.mean())
    return drift, rmean, (abs(drift) <= 0.2 and abs(rmean - cexp) <= 0.25)

def zeros(pari, spec, T2, poles=0):
    t0 = time.time()
    L = build(pari, spec)
    Q, k, mu = params(pari, L)
    eps = to_c(pari.lfunrootres(L)[2])
    phi = 0.0 if abs(eps.imag) < 1e-9 else math.atan2(eps.imag, eps.real)
    Li = pari.lfuninit(L, [1, int(T2) + 3])
    m0 = int(pari.lfunorderzero(L)) if spec[0] != "zeta" else 0
    cexp = 0.5 + poles - m0 / 2 - phi / (2 * math.pi)
    g = np.sort(np.array([float(x) for x in pari.lfunzeros(Li, [0.01, T2 + 1.0])]))
    n0 = len(g); added = 0
    for it in range(3):
        u = unfold(g, Q, mu); gaps = np.diff(u); new = []
        for i in np.argsort(-gaps)[:30]:
            new += [float(x) for x in pari.lfunzeros(Li, [g[i] + 1e-7, g[i + 1] - 1e-7], 256 if gaps[i] > 2.2 else 64)]
        new += [float(x) for x in pari.lfunzeros(Li, [0.01, g[0] - 1e-7], 64)]
        if not new: break
        g = np.sort(np.concatenate([g, new])); added += len(new)
    keep = np.concatenate([[True], np.diff(g) > 1e-7]); ndup = int((~keep).sum()); g = g[keep]
    drift, rmean, ok = check(g, Q, mu, cexp)
    extra = 0
    if not ok:
        u = unfold(g, Q, mu); gaps = np.diff(u); new = []
        for i in np.nonzero(gaps > 1.8)[0]:
            new += [float(x) for x in pari.lfunzeros(Li, [g[i] + 1e-7, g[i + 1] - 1e-7], 256)]
        if new:
            g = np.sort(np.concatenate([g, new])); keep = np.concatenate([[True], np.diff(g) > 1e-7]); g = g[keep]; extra = len(new)
        drift, rmean, ok = check(g, Q, mu, cexp)
    info = dict(spec=spec, T2=T2, Q=Q, mu=mu, eps=[eps.real, eps.imag], m0=m0, poles=poles, c_exp=cexp, n_first=n0,
                added=added, added_extra=extra, duplicates=ndup, n=len(g), drift=drift, r_mean=rmean,
                status='COMPLETE' if ok else 'REJECTED', secs=round(time.time() - t0, 1), pari=str(pari('version()')))
    return g, info
