"""moment_test.py: the readings of the registered test (round_draft.txt) of the mean square of products of two and three
real primitive Dirichlet L-functions on [100, 1000] against the frozen recipe predictions.

usage:  moment_test.py PRED.json VALDIR [--zeta ZETADIR]
        PRED.json = predictions_triple.json (frozen; its sha256 is checked against predictions_triple.sha256 if present)
                    or predictions_dev.json (dry run on the development set)
        VALDIR    = output of moment_values.py for every character of PRED's set (absL_<D>.npy);
        ZETADIR   = directory holding absL_1.npy (zeta) for the zeta * L control (default: VALDIR)
Measured integrals: trapezoid rule on the grid t = 100 + k/20 (k = 0..18000) of |L_1 ... L_k|^2 (quadrature error
< 2e-5 relative: h = 0.05 vs 0.1 differs by 5e-5 and the error is O(h^2)).
Jackknife: delete one character (all tuples containing it), n = number of characters.
"""
import sys, os, json, math, hashlib
import numpy as np

H = 0.05


def trap(y):
    return H * (y.sum() - 0.5 * (y[0] + y[-1]))


def jackknife(chars, tuples, stat):
    """delete-one-character jackknife of stat(index list)."""
    full = stat(list(range(len(tuples))))
    reps = []
    for c in chars:
        keep = [i for i, t in enumerate(tuples) if c not in t]
        reps.append(stat(keep))
    reps = np.array(reps); n = len(chars)
    se = math.sqrt((n - 1) / n * np.sum((reps - reps.mean(axis=0)) ** 2))
    return full, se


def ols_slope(x, y):
    x = np.asarray(x) - np.mean(x); y = np.asarray(y) - np.mean(y)
    return float(np.sum(x * y) / np.sum(x * x))


def main(argv):
    pf = argv[1]; vd = argv[2]
    zd = argv[argv.index('--zeta') + 1] if '--zeta' in argv else vd
    shaf = os.path.join(os.path.dirname(os.path.abspath(pf)), 'predictions_triple.sha256')
    if os.path.basename(pf) == 'predictions_triple.json' and os.path.exists(shaf):
        h = hashlib.sha256(open(pf, 'rb').read()).hexdigest(); reg = open(shaf).read().split()[0]
        print('prediction file sha256 %s: %s' % (h, 'MATCHES the registered hash' if h == reg else 'DOES NOT MATCH ' + reg))
    P = json.load(open(pf))
    chars = P['sets']['reserved']
    V = {}
    for D in chars:
        v = np.load(os.path.join(vd, 'absL_%d.npy' % D))
        assert v.shape == (18001,) and np.all(np.isfinite(v)), D
        V[D] = v ** 2
    have_zeta = os.path.exists(os.path.join(zd, 'absL_1.npy'))
    if have_zeta:
        Z = np.load(os.path.join(zd, 'absL_1.npy')) ** 2
    print('characters %d, pairs %d, triples %d; zeta control %s' % (len(chars), len(P['pairs']), len(P['triples']), 'yes' if have_zeta else 'no'))

    # ---------------------------------------------------------------- controls (theorems)
    S = {r['key'][0]: r for r in P['singles']}
    r1 = np.array([trap(V[D]) / S[D]['P1'] for D in chars])
    r4 = np.array([trap(V[D] ** 2) / S[D]['P4'] for D in chars])
    print('\nC1 single L, second moment (theorem):   mean ratio %.4f +- %.4f  (sd %.4f, min %.4f, max %.4f)' % (r1.mean(), r1.std(ddof=1) / math.sqrt(len(r1)), r1.std(ddof=1), r1.min(), r1.max()))
    print('C2 single L, fourth moment (theorem):   mean ratio %.4f +- %.4f  (sd %.4f)' % (r4.mean(), r4.std(ddof=1) / math.sqrt(len(r4)), r4.std(ddof=1)))
    if have_zeta:
        mz = np.array([trap(Z * V[D]) for D in chars])
        fz = np.array([S[D]['zetaL']['full'] for D in chars]); dz = np.array([S[D]['zetaL']['diag'] for D in chars])
        rz = mz / fz
        gz = float(np.sum(np.log(mz / dz) * np.log(fz / dz)) / np.sum(np.log(fz / dz) ** 2))
        print('C3 zeta * L(chi) (theorem, mixed swaps up to 26 %%): mean ratio to recipe %.4f (sd %.4f); to AFE diagonal %.4f; '
              'gamma_diag %.3f (1 = recipe, 0 = diagonal; naive se %.3f; zeta is common to all, so the se is optimistic)' % (
                  rz.mean(), rz.std(ddof=1), (mz / dz).mean(), gz,
                  np.std(np.log(rz)) / math.sqrt(np.sum(np.log(fz / dz) ** 2))))

    # ---------------------------------------------------------------- pairs (control)
    pairs = [tuple(r['key']) for r in P['pairs']]
    PP = {tuple(r['key']): r for r in P['pairs']}
    M2 = {p: trap(V[p[0]] * V[p[1]]) for p in pairs}
    r2 = np.array([M2[p] / PP[p]['full'] for p in pairs])
    m, se = jackknife(chars, pairs, lambda k: r2[k].mean())
    print('\nP1 pairs: pooled ratio measured/recipe %.4f +- %.4f (jackknife over characters); sd per pair %.4f' % (m, se, r2.std(ddof=1)))
    lm2 = np.log([M2[p] for p in pairs]); lp2 = np.log([PP[p]['full'] for p in pairs])
    s, se = jackknife(chars, pairs, lambda k: ols_slope(lp2[k], lm2[k]))
    print('P2 pairs: slope of ln M on ln P_recipe %.4f +- %.4f' % (s, se))
    for c in ['diag', 'lead']:
        rc = np.array([M2[p] / PP[p][c] for p in pairs]); m, se = jackknife(chars, pairs, lambda k: rc[k].mean())
        print('    comparator %-5s pooled ratio %.4f +- %.4f (H_%s would read 1; recipe predicts %.4f)' % (c, m, se, c, np.mean([PP[p]['full'] / PP[p][c] for p in pairs])))

    # ---------------------------------------------------------------- triples (test)
    tri = [tuple(r['key']) for r in P['triples']]
    TT = {tuple(r['key']): r for r in P['triples']}
    M3 = np.array([trap(V[a] * V[b] * V[c]) for a, b, c in tri])
    pf3 = np.array([TT[t]['full'] for t in tri])
    r3 = M3 / pf3
    cpair = np.array([M2[(a, b)] / PP[(a, b)]['full'] + M2[(a, c)] / PP[(a, c)]['full'] + M2[(b, c)] / PP[(b, c)]['full'] - 3 for a, b, c in tri])
    zlog = np.array([math.log(M2[(a, b)] / PP[(a, b)]['full']) + math.log(M2[(a, c)] / PP[(a, c)]['full']) + math.log(M2[(b, c)] / PP[(b, c)]['full']) for a, b, c in tri])
    comps = ['full', 'diag', 'kirk', 'no3', 'ind', 'lead']
    PC = {c: np.array([TT[t][c] for t in tri]) for c in comps}
    print('\nT1 triples, level (primary for the polynomial):')
    m, se = jackknife(chars, tri, lambda k: r3[k].mean())
    mc, sec = jackknife(chars, tri, lambda k: (r3[k] - cpair[k]).mean())
    print('   raw pooled ratio measured/recipe          %.4f +- %.4f   (sd per triple %.3f)' % (m, se, r3.std(ddof=1)))
    print('   control-variate (minus sum of pair ratio excesses, beta = 1)  %.4f +- %.4f' % (mc, sec))
    print('   expected under H_c (mean P_c/P_full):  ' + '  '.join('%s %.4f' % (c, np.mean(PC[c] / pf3)) for c in comps))
    print('\nT2 three-body Euler factor (primary for the arithmetic): gamma_c = centred OLS slope of (ln M - z - ln P_c) on (ln P_full - ln P_c), z = sum of pair log-ratios')
    lm3 = np.log(M3)
    for c in ['no3', 'kirk', 'diag', 'ind', 'lead']:
        xc = np.log(PC[c]); xf = np.log(pf3)
        g, se = jackknife(chars, tri, lambda k: ols_slope(xf[k] - xc[k], lm3[k] - zlog[k] - xc[k]))
        g0, se0 = jackknife(chars, tri, lambda k: ols_slope(xf[k] - xc[k], lm3[k] - xc[k]))
        print('   vs %-5s gamma (CV) %.4f +- %.4f   raw %.4f +- %.4f   [H_full 1, H_%s 0; sd of (x_f - x_c) %.3f]' % (c, g, se, g0, se0, c, np.std(xf - xc)))
    print('\nT3 slope of ln M on ln P_full across triples:')
    s, se = jackknife(chars, tri, lambda k: ols_slope(np.log(pf3[k]), lm3[k]))
    s2, se2 = jackknife(chars, tri, lambda k: ols_slope(np.log(pf3[k]), lm3[k] - zlog[k]))
    print('   raw %.4f +- %.4f   CV %.4f +- %.4f   (H_full 1; under H_c the expected slope is the slope of ln P_c on ln P_full: %s)' % (
        s, se, s2, se2, ' '.join('%s %.3f' % (c, ols_slope(np.log(pf3), np.log(PC[c]))) for c in comps[1:])))
    print('\nT4 comparator levels: mean of M/P_c (raw | CV) and sd of ln(M/P_c):')
    for c in comps:
        rc = M3 / PC[c]
        m, se = jackknife(chars, tri, lambda k: rc[k].mean())
        mcv, secv = jackknife(chars, tri, lambda k: (rc[k] - cpair[k] * PC['full'][k] / PC[c][k]).mean())
        print('   %-5s %.4f +- %.4f | %.4f +- %.4f   sd(ln) %.3f' % (c, m, se, mcv, secv, np.std(np.log(rc), ddof=1)))
    print('\nT5 (secondary) by the p = 2 configuration of the triple (m: chi(2) = -1, p: +1, r: 2 | D):')
    p2 = [TT[t]['p2'] for t in tri]
    for g in sorted(set(p2)):
        k = [i for i, x in enumerate(p2) if x == g]
        if len(k) < 15:
            continue
        pred = np.mean(np.log(pf3[k] / PC['no3'][k])); meas = np.mean(lm3[k] - zlog[k] - np.log(PC['no3'][k]))
        print('   %s n = %4d  predicted mean ln(P_full/P_no3) %+.3f   measured mean (ln M - z - ln P_no3) %+.3f' % (g, len(k), pred, meas))


if __name__ == '__main__':
    main(sys.argv)
