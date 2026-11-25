"""degree-4 follow-up: builders and the registered zero routine with the CORRECTED completeness rule.
Spec kinds: ['sym', ainv, m] | ['g2', f, h] | ['rs', ainv1, ainv2, eps] | ['char', q, m] (checks only).
Rankin-Selberg L(E1 x E2), gcd(N1, N2) = 1: Dirichlet series from local factors (arithmetic, weight 3, centre 3/2):
 good p: 1 - ab X + (p a^2 + p b^2 - 2 p^2) X^2 - p^2 ab X^3 + p^4 X^4; p | N1: 1 - ab X + a^2 p X^2; p | N2: 1 - ab X + b^2 p X^2;
 Gamma_C(s) Gamma_C(s-1) (Vga [-1,0,0,1]), conductor (N1 N2)^2; eps fixed by the functional equation (the other sign fails).
Completeness rule (registered; corrected after round 151): as round 151, except the offset
 c_exp = 1/2 + (poles at s = 1) - m0/2 - A0/pi, A0 = arg of L(1/2 + delta)/delta^m0 followed continuously along the real axis
 from sigma = 12 down to sigma = 1/2 + delta (delta = 1e-6; 150 points, phase unwrapped) for complex L; for real L
 (real coefficients and eps) A0 = 0 under GRH, checked by the sign of L(1/2 + 1e-3) (one evaluation). For real L this is the
 round-151 rule; for complex L it replaces -arg(eps)/2pi, which can be on the wrong branch."""
import math, time, json, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from ulib import theta, to_c

RSDEF = ('rsloc(E1,E2,N1,N2,p) = my(a = ellap(E1,p), b = ellap(E2,p)); if(N1%p==0, 1/(1 - a*b*X + a^2*p*X^2), '
         'if(N2%p==0, 1/(1 - a*b*X + b^2*p*X^2), 1/(1 - a*b*X + (a^2*p + b^2*p - 2*p^2)*X^2 - a*b*p^2*X^3 + p^4*X^4)))')

def build(pari, spec):
    k = spec[0]
    if k == 'sym':
        return pari.lfunsympow(pari.ellinit(spec[1]), spec[2])
    if k == 'g2':
        return pari('lfungenus2([%s,%s])' % (spec[1], spec[2]))
    if k == 'char':
        return pari('lfuncreate([G=znstar(%d,1), znconreychar(G,%d)])' % (spec[1], spec[2]))
    if k == 'rs':
        pari(RSDEF)
        return pari('my(E1 = ellinit(%s), E2 = ellinit(%s), N1 = ellglobalred(E1)[1], N2 = ellglobalred(E2)[1]); '
                    'lfuncreate([n -> direuler(p=2, n, rsloc(E1,E2,N1,N2,p)), 0, [-1,0,0,1], 3, (N1*N2)^2, %d])' % (spec[1], spec[2], spec[3]))
    raise ValueError(spec)

def params(pari, L):
    N, k, V = pari.lfunparams(L); k = float(k)
    return float(N), k, [float(v) + (k - 1) / 2 for v in V]

def unfold(g, Q, mu): return theta(g, Q, mu) / math.pi

def arg_centre(pari, L, k, m0, real=None):
    """A0: continuous argument of L(s)/(s-1/2)^m0 along s = sigma (analytic) from 12 down to 1/2 + 1e-6.
    Real L (real coefficients, real eps): L(sigma) is real and, under GRH, has no zero right of 1/2, so A0 = 0 when
    L(1/2 + 1e-3) / 1e-3^m0 > 0 (one evaluation; A0 = pi and a flag otherwise). Complex L: 150-point unwrap."""
    if real:
        v = to_c(pari.lfun(L, 0.5 + 1e-3 + (k - 1) / 2))
        return (0.0 if v.real > 0 else math.pi), abs(v)
    sig = 0.5 + np.concatenate([np.geomspace(11.5, 1e-2, 100), np.geomspace(1e-2, 1e-6, 50)[1:]])
    ph = []
    for s in sig:
        v = to_c(pari.lfun(L, s + (k - 1) / 2))
        ph.append(math.atan2(v.imag, v.real))
    ph = np.unwrap(np.array(ph))
    return float(ph[-1]), float(abs(to_c(pari.lfun(L, sig[-1] + (k - 1) / 2))))

def check(g, Q, mu, cexp):
    r = np.arange(1, len(g) + 1) - unfold(g, Q, mu); m = max(5, int(0.15 * len(g)))
    drift = float(r[-m:].mean() - r[:m].mean()); rmean = float(r.mean())
    return drift, rmean, (abs(drift) <= 0.2 and abs(rmean - cexp) <= 0.25)

def zeros(pari, spec, T2, poles=0):
    t0 = time.time()
    L = build(pari, spec)
    Q, k, mu = params(pari, L)
    eps = to_c(pari.lfunrootres(L)[2])
    m0 = int(pari.lfunorderzero(L))
    A0, _ = arg_centre(pari, L, k, m0, real=(abs(eps.imag) < 1e-12 and spec[0] != 'char'))
    cexp = 0.5 + poles - m0 / 2 - A0 / math.pi
    Li = pari.lfuninit(L, [1, int(T2) + 3])
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
    info = dict(spec=spec, T2=T2, Q=Q, mu=mu, eps=[eps.real, eps.imag], m0=m0, A0=A0, poles=poles, c_exp=cexp, n_first=n0,
                added=added, added_extra=extra, duplicates=ndup, n=len(g), drift=drift, r_mean=rmean,
                status='COMPLETE' if ok else 'REJECTED', secs=round(time.time() - t0, 1), pari=str(pari('version()')))
    return g, info
