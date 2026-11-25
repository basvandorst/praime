"""Readings beyond the registered ones (own code; frozen recipe coefficients hm / hm_no3 / g / c4, verified on 20 triples
by myrecipe.py).  Sections:
 A  off-diagonal fraction against H_diag, and H_diag under other AFE truncations (smooth weight V^2, weight V,
    unbalanced sharp splits)
 B  lower-order terms: amplitudes of the |m| = 1, 2, 3 parts; 'pairs-only' constructions
 C  thirds of [100, 1000]: level, CV level, gamma_no3, singles, fourth moments, pairs
 D  subgroups: conductor size, shared primes, signs; skewness"""
import json, math, os, itertools
import numpy as np
from scipy.special import erfc
from scipy.integrate import quad
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
import sys; sys.path.insert(0, HERE)
from stats import jk, slope, P, R, tri, pairs, TT, PP, SS, m
EG = float(mp.euler)
EDGES = [(100.0, 1000.0), (100.0, 400.0), (400.0, 700.0), (700.0, 1000.0)]
xg, wg = np.polynomial.legendre.leggauss(40)


def gl(f, a, b):
    u = 0.5 * (math.log(b) - math.log(a)) * xg + 0.5 * (math.log(b) + math.log(a)); t = np.exp(u)
    return float(np.sum(wg * f(t) * t) * 0.5 * (math.log(b) - math.log(a)))


MONOS3 = list(itertools.product((0, 1), repeat=3)); MONOS2 = list(itertools.product((0, 1), repeat=2)); MONOS1 = [(0,), (1,)]


def Pdens(Ds, hm, t, monos, order=None):
    """recipe-type density from multilinear coefficients; order = set of |m| to keep (None: all)."""
    Ls = [np.log(abs(D) * t / (2 * math.pi)) + 2 * EG for D in Ds]
    out = 0.0
    for c, mm in zip(hm, monos):
        if order is not None and sum(mm) not in order:
            continue
        term = c
        for j, e in enumerate(mm):
            if e == 0:
                term = term * Ls[j]
        out = out + term
    return out


ST = [float(mp.stieltjes(n)) for n in range(4)]


def szeta_pow(k):
    zs = [1.0] + [(-1) ** n * ST[n] / math.factorial(n) for n in range(k)]
    zk = [1.0] + [0.0] * k
    for _ in range(k):
        new = [0.0] * (k + 1)
        for a in range(k + 1):
            for b in range(k + 1 - a):
                new[a + b] += zk[a] * zs[b]
        zk = new
    return zk


Z3 = szeta_pow(3)


def diag_dens(Ds, g, t, W=(0, 0, 0), split=None):
    """product-AFE diagonal: 2 Res zeta(1+s)^3 g(s) C^{s/2} Wt(s) (W: Laurent coeffs w0,w1,w2 of Wt - 1/s);
    split = theta: sharp cutoffs C^theta and C^(1-theta) instead of sqrt(C) twice."""
    k = 3
    c = [sum(Z3[a] * g[n - a] for a in range(n + 1)) for n in range(k + 1)]
    lnC = sum(np.log(abs(D) * t / (2 * math.pi)) for D in Ds)
    Wc = [1.0] + list(W)

    def coef(ell):
        E = [ell ** n / math.factorial(n) for n in range(k + 1)]
        tot = 0.0
        for a in range(k + 1):
            for b in range(k + 1 - a):
                tot = tot + c[a] * E[b] * Wc[k - a - b]
        return tot
    if split is None:
        return 2 * coef(lnC / 2)
    return coef(split * lnC) + coef((1 - split) * lnC)


def laurent_W(Vfun):
    """w0, w1, w2 of the Mellin transform of W(x) = Vfun(u = ln x) (minus 1/s)."""
    out = []
    for n in range(3):
        a = quad(lambda u: (Vfun(u) - 1) * u ** n, -60, 0, limit=400)[0]
        b = quad(lambda u: Vfun(u) * u ** n, 0, 60, limit=400)[0]
        out.append((a + b) / math.factorial(n))
    return out


def main():
    rp = {p: m['pair_trap'][i] / PP[p]['full'] for i, p in enumerate(pairs)}
    M3 = m['triple_trap']
    pf = np.array([TT[t]['full'] for t in tri]); r3 = M3 / pf
    cpair = np.array([rp[(a, b)] + rp[(a, c)] + rp[(b, c)] - 3 for a, b, c in tri])
    zlog = np.array([math.log(rp[(a, b)]) + math.log(rp[(a, c)]) + math.log(rp[(b, c)]) for a, b, c in tri])
    dratio = np.array([TT[t]['diag'] / TT[t]['full'] for t in tri])
    print('=== A. off-diagonal fraction (lambda = 1: recipe, 0: H_diag sharp) ===')
    for nm, y in (('CV', r3 - cpair), ('raw', r3)):
        lam = jk(lambda k: 1 - (y[k].mean() - 1) / (dratio[k].mean() - 1), tri)
        print('  %-3s lambda = %.3f +- %.3f   (mean P_diag/P_full %.4f)' % (nm, lam[0], lam[1], dratio.mean()))
    # H_diag variants: recompute the sharp one from g as a check, then smooth weights and unbalanced splits
    variants = {'sharp (registered)': dict()}
    for d in (0.5, 1.0, 2.0, 4.0):
        Vf = lambda u, d=d: (0.5 * erfc(d * u / 2)) ** 2
        variants['smooth V^2, G=exp(w^2/%.1f^2)' % d] = dict(W=laurent_W(Vf))
        Vf1 = lambda u, d=d: 0.5 * erfc(d * u / 2)
        variants['weight V (not squared), delta %.1f' % d] = dict(W=laurent_W(Vf1))
    for th in (0.45, 0.4, 0.3):
        variants['sharp unbalanced theta %.2f' % th] = dict(split=th)
    for name, kw in variants.items():
        rat = []
        for t in tri:
            r = TT[t]; Ds = r['key']
            Id = gl(lambda tt: diag_dens(Ds, r['g'], tt, **kw), 100, 1000)
            rat.append(Id / r['full'])
        rat = np.array(rat)
        if name.startswith('sharp (reg'):
            chk = np.max(np.abs(rat / dratio - 1))
        lvl = jk(lambda k: (r3[k] - cpair[k]).mean() / rat[k].mean(), tri)
        print('  %-42s mean P_c/P_full %.4f [%.3f, %.3f]  W=%s   data/H_c (CV) %.4f +- %.4f' % (
            name, rat.mean(), rat.min(), rat.max(), np.round(kw.get('W', []), 4).tolist(), lvl[0], lvl[1]))
    print('  check: my sharp H_diag from frozen g vs frozen diag integral: max rel diff %.1e' % chk)

    print('=== B. lower-order terms ===')
    parts = {k: [] for k in range(4)}
    for t in tri:
        r = TT[t]; Ds = r['key']
        for k in range(4):
            parts[k].append(gl(lambda tt: Pdens(Ds, r['hm'], tt, MONOS3, {k}), 100, 1000))
    parts = {k: np.array(v) for k, v in parts.items()}
    tot = sum(parts.values()); print('  check: sum of |m| parts / frozen full: max rel %.1e' % np.max(np.abs(tot / pf - 1)))
    for k in range(4):
        sh = parts[k] / pf
        print('  |m| = %d share of int P: median %+.3f  [%+.3f, %+.3f]  sd %.3f' % (k, np.median(sh), sh.min(), sh.max(), sh.std()))
    y = np.log(M3) - zlog
    for k in (1, 2, 3):
        x = parts[k] / pf
        # amplitude lambda_k: ln M - z = ln(P_full + (lam-1) part_k) ~ ln P_full + (lam-1) x ; centred slope -> lam - 1
        s = jk(lambda kk: slope(x[kk], y[kk] - np.log(pf[kk])), tri)
        # level version (no intercept): mean of (M - cpair Pfull - (Pfull - part_k)) / part_k is unstable; use ratio of means
        lv = jk(lambda kk: 1 - ((r3[kk] - cpair[kk]).mean() - 1) / x[kk].mean(), tri)
        print('  amplitude of |m| = %d terms: shape (centred slope, CV) %.3f +- %.3f ; from the level (CV) %.3f +- %.3f' % (
            k, 1 + s[0], s[1], lv[0], lv[1]))
    # alternative constructions with the same h(0) (incl. A_3): kirk * A3ratio, no3 * A3ratio, only m = {} term, drop |m|>=2
    A3 = np.array([TT[t]['A3_over_Aij'] for t in tri])
    alts = {'kirk x A3/prodAij (pairs-only lower order)': np.array([TT[t]['kirk'] for t in tri]) * A3,
            'no3 x A3/prodAij (no three-body derivatives)': np.array([TT[t]['no3'] for t in tri]) * A3,
            'h(0) prod(L_j + 2 gamma) (|m| = 0 only)': parts[0],
            'drop |m| >= 2': parts[0] + parts[1],
            'lead h(0) prod L_j': np.array([TT[t]['lead'] for t in tri])}
    lm = np.log(M3)
    for name, pc in alts.items():
        xf = np.log(pf); xc = np.log(pc)
        g = jk(lambda k: slope(xf[k] - xc[k], lm[k] - zlog[k] - xc[k]), tri)
        lv = jk(lambda k: ((M3[k] - cpair[k] * pf[k]) / pc[k]).mean(), tri)
        lv0 = jk(lambda k: (M3[k] / pc[k]).mean(), tri)
        print('  %-46s P_c/P_full mean %.4f sd(ln) %.4f | gamma (1 recipe, 0 alt) %.3f +- %.3f | level CV %.4f +- %.4f raw %.4f +- %.4f' % (
            name, np.mean(pc / pf), np.std(np.log(pc / pf)), g[0], g[1], lv[0], lv[1], lv0[0], lv0[1]))

    print('=== C. thirds of [100, 1000] ===')
    M3p = m['triple_parts']; M2p = m['pair_parts']; M1p = m['single_parts']; M4p = m['single4_parts']
    for ip, (a, b) in enumerate(EDGES[1:]):
        p3 = np.array([gl(lambda tt: Pdens(TT[t]['key'], TT[t]['hm'], tt, MONOS3), a, b) for t in tri])
        n3 = np.array([gl(lambda tt: Pdens(TT[t]['key'], TT[t]['hm_no3'], tt, MONOS3), a, b) for t in tri])
        p2 = {p: gl(lambda tt: Pdens(p, PP[p]['hm'], tt, MONOS2), a, b) for p in pairs}
        rpp = {p: M2p[i, ip] / p2[p] for i, p in enumerate(pairs)}
        r2 = np.array([rpp[p] for p in pairs])
        p1 = np.array([gl(lambda tt: Pdens((D,), SS[D]['hm1'], tt, MONOS1), a, b) for D in R])
        p4 = np.array([gl(lambda tt: sum(c * np.log(tt) ** kk for kk, c in enumerate(SS[D]['c4'])), a, b) for D in R])
        r1 = M1p[:, ip] / p1; r4 = M4p[:, ip] / p4
        rr = M3p[:, ip] / p3
        cp = np.array([rpp[(x, y)] + rpp[(x, z)] + rpp[(y, z)] - 3 for x, y, z in tri])
        zl = np.array([math.log(rpp[(x, y)] * rpp[(x, z)] * rpp[(y, z)]) for x, y, z in tri])
        raw = jk(lambda k: rr[k].mean(), tri); cv = jk(lambda k: (rr[k] - cp[k]).mean(), tri)
        pr = jk(lambda k: r2[k].mean(), pairs)
        xf = np.log(p3); xc = np.log(n3); lmm = np.log(M3p[:, ip])
        gno3 = jk(lambda k: slope(xf[k] - xc[k], lmm[k] - zl[k] - xc[k]), tri)
        print('  [%4.0f,%4.0f] triples raw %.4f +- %.4f  CV %.4f +- %.4f  gamma_no3 %.3f +- %.3f | pairs %.4f +- %.4f | singles %.4f +- %.4f | fourth %.4f +- %.4f' % (
            a, b, raw[0], raw[1], cv[0], cv[1], gno3[0], gno3[1], pr[0], pr[1], r1.mean(), r1.std(ddof=1) / math.sqrt(56),
            r4.mean(), r4.std(ddof=1) / math.sqrt(56)))

    print('=== D. subgroups and skewness ===')
    from math import gcd
    q3 = np.array([abs(t[0] * t[1] * t[2]) for t in tri])
    nsh = np.array([(gcd(abs(t[0]), abs(t[1])) > 1) + (gcd(abs(t[0]), abs(t[2])) > 1) + (gcd(abs(t[1]), abs(t[2])) > 1) for t in tri])
    nneg = np.array([sum(D < 0 for D in t) for t in tri])
    cv = r3 - cpair
    def grp(name, lab):
        for g_ in sorted(set(lab)):
            k = np.nonzero(lab == g_)[0]
            sub = [tri[i] for i in k]
            a1 = jk(lambda kk: r3[k][kk].mean(), sub); a2 = jk(lambda kk: cv[k][kk].mean(), sub)
            print('  %-10s %-6s n %4d  raw %.4f +- %.4f  CV %.4f +- %.4f' % (name, g_, len(k), a1[0], a1[1], a2[0], a2[1]))
    terc = np.digitize(q3, np.quantile(q3, [1 / 3, 2 / 3]))
    grp('q1q2q3', terc); grp('n shared', nsh); grp('n neg D', nneg)
    from scipy.stats import skew
    print('  r3: mean %.4f median %.4f sd %.4f skew %.2f ; CV quantity: mean %.4f median %.4f skew %.2f' % (
        r3.mean(), np.median(r3), r3.std(), skew(r3), cv.mean(), np.median(cv), skew(cv)))
    r2 = np.array([rp[p] for p in pairs]); r1 = m['single_trap'] / np.array([SS[D]['P1'] for D in R])
    print('  pairs: mean %.4f median %.4f skew %.2f ; singles: mean %.4f median %.4f skew %.2f' % (
        r2.mean(), np.median(r2), skew(r2), r1.mean(), np.median(r1), skew(r1)))
    # influence of characters on the level (jackknife pseudo-values)
    reps = []
    for c in R:
        k = np.array([i for i, t in enumerate(tri) if c not in t]); reps.append((c, r3[k].mean(), cv[k].mean()))
    reps.sort(key=lambda x: x[1])
    print('  deleting one character moves raw level to [%.4f (del %d), %.4f (del %d)]; CV [%.4f, %.4f]' % (
        reps[0][1], reps[0][0], reps[-1][1], reps[-1][0], min(x[2] for x in reps), max(x[2] for x in reps)))


if __name__ == '__main__':
    main()
