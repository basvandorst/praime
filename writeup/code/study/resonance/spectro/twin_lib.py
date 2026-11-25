"""twin_lib.py: chunk sums for the constellation tests (twin primes, labels on shifted primes).

For a chunk [L, U) (L >= 0) and a Hann window [v0, v1] in v = sqrt(n), labels in LABELS, statistics (complex sums):
  twinlo   sum_{p, p+2 prime}  lambda(p)   ln p ln(p+2) p^{-1/4} W(sqrt p) e(-alpha 2 sqrt(p/N))            -> 2C2 c MT
  twinhi   sum_{p, p+2 prime}  lambda(p+2) ln p ln(p+2) (p+2)^{-1/4} W(sqrt(p+2)) e(-alpha 2 sqrt((p+2)/N))  -> 2C2 c MT
  single   sum_p lambda(p) ln p p^{-1/4} W e(...)                                                           -> c MT
  shift2   sum_p lambda(p+2) ln p (p+2)^{-1/4} W(sqrt(p+2)) e(-alpha 2 sqrt((p+2)/N))   (lattice labels only)  -> c_shift MT
  sg       sum_{p, 2p+1 prime} lambda(p) ln p ln(2p+1) p^{-1/4} W e(...)                                     -> 2C2 c MT
for alpha = 1 (the line) and the 8 control factors CTRL; plus the noise sums sum w^2 and the sharp sums (W = 1, weight
without p^{-1/4}) at alpha = 1 for the running plots.  MT = main_term_unit (constant-1 law) of spectro_lib.
The chunk covers primes p with L <= p < U (for twinhi / shift2 the label sits at p + 2)."""
import math
import numpy as np
import spectro_lib as S

CTRL = [0.6713, 0.7919, 0.8837, 0.9431, 1.0617, 1.1237, 1.2389, 1.3711]
ALPHAS = [1.0] + CTRL
LABELS = ['EC:11a1', 'EC:15a1', 'A4', 'C23']
LATTICE = ['A4', 'C23']
STATS = ['twinlo', 'twinhi', 'single', 'shift2', 'shift6', 'sg']

APTAB = {}   # name -> (p, a) tables (discovery range); otherwise PARI ellap at the needed primes
def lam_on(label, nvals, L, U, lattice_cache):
    """lambda at the integers nvals (lattice labels: inside [L, U + 4); curves: nvals are primes)."""
    if label in LATTICE:
        arr = lattice_cache[label]; return arr[nvals - L]
    name = label[3:]
    if len(nvals) == 0: return np.zeros(0)
    if name in APTAB and APTAB[name][0][-1] >= nvals.max():
        p, a = APTAB[name]; k = np.searchsorted(p, nvals); assert np.all(p[k] == nvals); a = a[k]
    else:
        a = S.ap_at(name, nvals)
    return a / np.sqrt(nvals.astype(np.float64))

def chunk_sums(L, U, v0, v1, labels=LABELS, stats=STATS, lines=None):
    """lines: optional {lattice label: [alpha, ...]}: also sum_p lambda(p) ln p Hann(sqrt p; v0, v1) e(-alpha sqrt p) and the
    noise sum (returned under out[label]['lines'])."""
    isp = S.is_prime_seg(L, U + 4)            # p and p + 2
    pr = np.nonzero(isp[:U - L])[0].astype(np.int64) + L
    tw = pr[isp[pr + 2 - L]]
    del isp
    sgp = np.zeros(0, np.int64)
    if 'sg' in stats:
        isq = S.is_prime_seg(2 * L, 2 * U + 2)   # 2p + 1
        sgp = pr[isq[2 * pr + 1 - 2 * L]]
        del isq
    lattice_cache = {}
    for lab in labels:
        if lab == 'A4': lattice_cache[lab] = S.lam_gauss_seg(L, U + 8, 4)
        if lab == 'C23': lattice_cache[lab] = S.lam_c23_seg(L, U + 8)
    out = {}
    for lab in labels:
        N = S.label_info(lab)['N']; rec = {}
        for st in stats:
            if st in ('shift2', 'shift6') and lab not in LATTICE: continue
            if st == 'twinlo': n = tw; wt = np.log(tw) * np.log(tw + 2.0); lam = lam_on(lab, n, L, U, lattice_cache)
            elif st == 'twinhi': n = tw + 2; wt = np.log(tw) * np.log(tw + 2.0); lam = lam_on(lab, n, L, U, lattice_cache)
            elif st == 'single': n = pr; wt = np.log(pr.astype(np.float64)); lam = lam_on(lab, n, L, U, lattice_cache)
            elif st == 'shift2': n = pr + 2; wt = np.log(pr.astype(np.float64)); lam = lattice_cache[lab][n - L]
            elif st == 'shift6': n = pr + 6; wt = np.log(pr.astype(np.float64)); lam = lattice_cache[lab][n - L]
            elif st == 'sg': n = sgp; wt = np.log(sgp.astype(np.float64)) * np.log(2.0 * sgp + 1); lam = lam_on(lab, n, L, U, lattice_cache)
            nf = n.astype(np.float64); v = np.sqrt(nf)
            w = lam * wt * nf ** -0.25 * S.hann(v, v0, v1)
            ph = S.TWO_PI * 2 * np.sqrt(nf / N)
            sm = [complex(np.dot(w, np.exp(-1j * a * ph))) for a in ALPHAS]
            ws = lam * wt
            sharp = complex(np.dot(ws, np.exp(-1j * ph)))
            rec[st] = dict(smooth=[[z.real, z.imag] for z in sm], noise=float(np.dot(w, w)), sharp=[sharp.real, sharp.imag],
                           noise_sharp=float(np.dot(ws, ws)), count=int(len(n)))
        if lines and lab in lines:
            n = pr; nf = n.astype(np.float64); v = np.sqrt(nf); lam = lattice_cache[lab][n - L]
            w = lam * np.log(nf) * S.hann(v, v0, v1); k = w != 0; w = w[k]; v = v[k]
            ls = [complex(np.dot(w, np.exp(-1j * S.TWO_PI * a * v))) for a in lines[lab]]
            rec['lines'] = dict(S=[[z.real, z.imag] for z in ls], noise=float(np.dot(w, w)))
        out[lab] = rec
    return out

def add_rec(a, b):
    if a is None: return b
    o = {}
    for lab in b:
        o[lab] = {}
        for st in b[lab]:
            x, y = a[lab][st], b[lab][st]
            o[lab][st] = dict(smooth=[[p[0] + q[0], p[1] + q[1]] for p, q in zip(x['smooth'], y['smooth'])], noise=x['noise'] + y['noise'],
                              sharp=[x['sharp'][0] + y['sharp'][0], x['sharp'][1] + y['sharp'][1]], noise_sharp=x['noise_sharp'] + y['noise_sharp'],
                              count=x['count'] + y['count'])
    return o

# ------------------------------------------------------------------ predicted constants
TWO_C2 = 1.3203236316937391   # 2 prod_{p>2} (1 - 1/(p-1)^2)

def ftab(t, l):
    x = np.arange(l)
    if t in ('twinlo', 'sg_lo'):
        if l == 2: return 4.0 * (x == 1)
        bad2 = (l - 2) % l if t == 'twinlo' else ((l - 1) * pow(2, -1, l)) % l   # n + 2 = 0  /  2n + 1 = 0
        return (l / (l - 1)) ** 2 * ((x != 0) & (x != bad2))
    if t == 'twinhi':
        if l == 2: return 4.0 * (x == 1)
        return (l / (l - 1)) ** 2 * ((x != 0) & (x != 2 % l))
    if t == 'shift2':
        if l == 2: return 2.0 * (x == 1)
        return (l / (l - 1)) * (x != 2 % l)
    if t == 'shift6':
        if l == 2: return 2.0 * (x == 1)
        return (l / (l - 1)) * (x != 6 % l)
    if t == 'single': return (l / (l - 1)) * (x != 0)
    raise ValueError(t)

def c_const(label, t, pmax=10 ** 6):
    """predicted constant of statistic t for a label (spectro.md, section 1): prod_l L_l(t), L_l = (1/l) sum_c f_l(c) R_l(c).
    For t in (twinlo, twinhi, sg_lo, single) it is computed as c_exact * prod_{l <= pmax} L_l(t)/L_l(single) (the ratio
    converges absolutely; c_exact = B/L(1, Ad) from curves_c.json / round 153); shift2 directly (absolutely convergent)."""
    loc = S.local_lams(label, pmax); prod = 1.0
    def Lfac(tt, l, lp, chi):
        if l <= 50 or (label == 'C23' and l == 23):
            return S.local_factor(label, l, lp, chi, ftab(tt, l))
        rho = S.rho_p(label, l, lp, chi)          # R uniform over the nonzero classes: L = mean (1 - rho) + rho f(0)
        if tt == 'single': mean, f0 = 1.0, 0.0
        elif tt in ('shift2', 'shift6'): mean, f0 = 1.0, l / (l - 1.0)
        else: mean, f0 = l * (l - 2.0) / (l - 1.0) ** 2, 0.0
        return mean * (1 - rho) + rho * f0
    for l, (lp, chi) in loc.items():
        Lt = Lfac(t, l, lp, chi)
        if t in ('shift2', 'shift6'): prod *= Lt
        else: prod *= Lt / Lfac('single', l, lp, chi)
    return prod if t in ('shift2', 'shift6') else S.label_info(label)['c'] * prod

def c_shift2(label, pmax=10 ** 6):
    return c_const(label, 'shift2', pmax)
