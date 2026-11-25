"""Spacing variance and gap ratio of the zeros of L(s, chi_D) per prime dividing D: the registered test.
Self-contained (numpy, scipy).  Reads the frozen predictions predictions_spacing.json (sha256 checked) and the zeros.

usage:
  spacing_test.py TEST [--zeros FILE [--no-completeness]] [--noise-only] [--out RESULT.json]   TEST in B Dn Fn C Bp R120
  spacing_test.py --pool RESULT.json ...                                       pooled readings over several tests
  spacing_test.py TEST --make-synthetic {cue,cuewave} --seed N --out FILE.npz synthetic zero sets for dry runs:
        cue     : GUE bulk spectra (the sine process locally) placed at the heights of the smooth zero count of chi_D
        cuewave : the same, re-mapped by the explicit-formula waves of the good primes among 2, 3, 5, 7, 11, 13
                  (random phase per prime and set): the wave law holds exactly for these primes, the level is GUE's.

Statistic per L-function (zeros lo <= gamma <= hi of the TEST's band): u = (gamma/2pi) ln(q gamma/(2 pi e)), s = diff(u),
Var = mean((s - mean s)^2), r~ = mean of min(s_n, s_n+1)/max(s_n, s_n+1) over consecutive pairs in the band.
Completeness (decided before any statistic; a flagged set is excluded and reported):
  (a) |N(t) - (theta(t) - theta(0))/pi| <= 1.5 at 200 equally spaced heights of [lo, hi] (ends included);
  (b) largest unfolded gap in the band <= 3.0;
  (c) the mean of N(t) - (theta(t) - theta(0))/pi over each eighth of the band (25 heights each) within +-0.4
      (on the used set R120 these means have sd 0.05; removing one zero moves the later eighths by -1, a close pair by -2);
  (d) smallest unfolded gap > 1e-6 (no duplicate).
"""
import sys, os, json, math, hashlib, argparse
import numpy as np
from scipy.special import loggamma

H = os.path.dirname(os.path.abspath(__file__))
REG_SHA = '0c6501d00771fc91467d45562db555da561e1cd5e28d591408a76bb006f0adee'
LZ = os.path.join(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..')), 'data', 'riemann', 'lzeros') + '/'
FILES = {'B': ['r134_B.npz'], 'Dn': ['r136_Dn.npz'], 'Fn': ['r138_Fn.npz'], 'C': ['r134_C.npz'], 'Bp': ['r138_Bp.npz'],
         'R120': ['r117_part0.npz', 'r117_part1.npz', 'r120_part0.npz', 'r120_part1.npz']}
PR = (2, 3, 5, 7, 11, 13); ALT = ('law', 'law2', 'slow', 'bblm', 'ignored'); NBOOT = 2000


def theta(t, q, a):
    return np.imag(loggamma(0.25 + a / 2 + 0.5j * np.asarray(t, float))) + np.asarray(t, float) / 2 * math.log(q / math.pi)


def completeness(g, D, lo, hi):
    q = abs(D); a = 0 if D > 0 else 1; t = np.linspace(lo, hi, 200)
    dev = np.searchsorted(g, t) - (theta(t, q, a) - theta(0.0, q, a)) / math.pi
    b = g[(g >= lo) & (g <= hi)]; u = b / (2 * math.pi) * np.log(q * b / (2 * math.pi * math.e)); d = np.diff(u)
    octm = dev.reshape(8, 25).mean(1)
    flags = []
    if np.abs(dev).max() > 1.5: flags.append('count %.2f' % np.abs(dev).max())
    if d.max() > 3.0: flags.append('gap %.2f' % d.max())
    if np.abs(octm).max() > 0.4: flags.append('eighth-mean %.2f' % octm[np.argmax(np.abs(octm))])
    if d.min() <= 1e-6: flags.append('duplicate')
    return flags, dict(maxdev=float(np.abs(dev).max()), maxgap=float(d.max()), mingap=float(d.min()), maxoct=float(np.abs(octm).max()), n_small=int((d < 0.1).sum()))


def stats(g, D, lo, hi):
    q = abs(D); b = g[(g >= lo) & (g <= hi)]
    u = b / (2 * math.pi) * np.log(q * b / (2 * math.pi * math.e)); s = np.diff(u)
    r = np.minimum(s[1:], s[:-1]) / np.maximum(s[1:], s[:-1])
    return float(np.mean((s - s.mean()) ** 2)), float(r.mean()), len(s)


def load_pred():
    raw = open(os.path.join(H, 'predictions_spacing.json')).read()
    h = hashlib.sha256(raw.encode()).hexdigest()
    if h != REG_SHA: sys.exit('predictions_spacing.json does not match the registered sha256 (%s)' % h)
    return json.loads(raw)


def load_zeros(test, zfile):
    Z = {}
    for f in ([zfile] if zfile else [LZ + x for x in FILES[test]]):
        z = np.load(f); Z.update({k: np.sort(z[k]) for k in z.files})
    return Z


def design(D): return np.array([[1.0] + [float(abs(d) % p == 0) for p in PR] for d in D])


# ------------------------------------------------------------------------------------------------ synthetic sets
def make_synthetic(P, test, kind, seed, out):
    rng = np.random.default_rng(seed); S = P[test]; lo, hi = S['height']; Z = {}
    for D in S['D']:
        q = abs(D); a = 0 if D > 0 else 1
        tt = np.linspace(25.0, hi + 4.0, 200001); c = (theta(tt, q, a) - theta(0.0, q, a)) / math.pi
        n0 = int(math.floor(c[0])); ntot = int(math.ceil(c[-1] - c[0])) + 2
        # GUE bulk (Dumitriu-Edelman tridiagonal model, beta = 2; scipy eigvalsh_tridiagonal), unfolded by the semicircle
        # law with the radius from the second moment; the central eigenvalues |lambda| < R/2 are kept (about 61 %)
        from scipy.linalg import eigvalsh_tridiagonal
        Ng = int(ntot / 0.6) + 50
        dg = rng.standard_normal(Ng) * math.sqrt(2); eg = np.sqrt(rng.chisquare(2 * np.arange(Ng - 1, 0, -1)))
        lam = eigvalsh_tridiagonal(dg / math.sqrt(2), eg / math.sqrt(2)); Rr = 2 * math.sqrt(np.mean(lam ** 2))
        lam = lam[np.abs(lam) < 0.5 * Rr]
        F = 0.5 + (lam * np.sqrt(Rr ** 2 - lam ** 2)) / (math.pi * Rr ** 2) + np.arcsin(lam / Rr) / math.pi
        y = np.sort(Ng * F); y = (y - y[0])[:ntot] + c[0]
        if kind == 'cuewave':
            good = [p for p in PR if q % p]; ph = {p: rng.uniform(0, 2 * math.pi) for p in good}
            x = y.copy()
            for _ in range(40):
                t = np.interp(x, c, tt)
                dN = sum(np.angle(1 - p ** -0.5 * np.exp(1j * (t * math.log(p) + ph[p]))) / math.pi for p in good)
                x = y - dN
            y = x
        g = np.interp(y[(y > c[0]) & (y < c[-1])], c, tt)
        dummy = np.linspace(1.0, 24.0, n0) if n0 > 0 else np.zeros(0)       # placeholders below the band (count check)
        Z['D%d' % D] = np.sort(np.concatenate([dummy, g]))
    np.savez(out, **Z); print('wrote', out, len(Z), 'sets')


# ------------------------------------------------------------------------------------------------ readings
def run(P, test, zfile, noise_only, outfile, nocomp=False):
    S = P[test]; lo, hi = S['height']; Z = load_zeros(test, zfile)
    keep, excl, V, R, N = [], [], [], [], []; comp = []
    for i, D in enumerate(S['D']):
        g = Z['D%d' % D]; fl, info = completeness(g, D, lo, hi); comp.append(info)
        if fl and not nocomp: excl.append((D, fl)); continue
        if fl: excl.append((D, 'flagged, kept (--no-completeness)'))
        v, r, n = stats(g, D, lo, hi); keep.append(i); V.append(v); R.append(r); N.append(n)
    V, R, N = map(np.array, (V, R, N)); D = [S['D'][i] for i in keep]; X = design(D); n = len(D)
    print('TEST %s (%s): %d L-functions, %d kept, heights %g-%g, mean spacings %.0f' % (test, zfile or 'registered zero files', len(S['D']), n, lo, hi, N.mean()))
    print('completeness: max |N - theta/pi| %.2f, max gap %.2f, max eighth-mean %.2f, min gap %.4f; excluded: %s'
          % (max(c['maxdev'] for c in comp), max(c['maxgap'] for c in comp), max(c['maxoct'] for c in comp), min(c['mingap'] for c in comp), ('%d flagged' % len(excl)) if nocomp else (excl or 'none')))
    A = np.linalg.inv(X.T @ X); Hm = A @ X.T
    if noise_only:
        for nm, y in (('Var', V), ('r~', R)):
            res = y - X @ (Hm @ y); print('NOISE ONLY  %-3s residual sd after the indicator regression: %.5f (dof %d)' % (nm, math.sqrt(res @ res / (n - X.shape[1])), n - X.shape[1]))
        rv = V - X @ (Hm @ V); rr = R - X @ (Hm @ R); print('NOISE ONLY  correlation of the residuals (Var, r~): %+.3f' % (np.corrcoef(rv, rr)[0, 1]))
        return
    pred = {v: {st: np.array(S['pred'][v][st])[keep] for st in ('var', 'r')} for v in ALT}
    rng = np.random.default_rng(1371); BI = [rng.integers(0, n, n) for _ in range(NBOOT)]
    res = {'test': test, 'n': n, 'excluded': excl, 'meanL': float(np.mean(np.array(S['meanL'])[keep]))}
    for st, y in (('var', V), ('r', R)):
        yy = y - pred['ignored'][st]; c = Hm @ yy
        cp = {v: Hm @ (pred[v][st] - pred['ignored'][st]) for v in ALT}
        # measurement error: bootstrap of the deviation from the law, d = OLS(stat - law); the regression-matched predictions
        # are computed on the same (full) sample, so their resampling spread is not an error of the reading
        dd = y - pred['law'][st]
        cb = np.array([np.linalg.lstsq(X[b], dd[b], rcond=None)[0] for b in BI])
        se = cb.std(0); C3 = np.cov(cb[:, 1:4].T)
        law = cp['law'][1:4]; Wd = 1 / se[1:4] ** 2
        gam = float((Wd * c[1:4] * law).sum() / (Wd * law * law).sum())
        gb = np.array([(Wd * x[1:4] * law).sum() / (Wd * law * law).sum() for x in cb])
        cb = cb + cp['law'][None, :]
        gexp = {v: float((Wd * cp[v][1:4] * law).sum() / (Wd * law * law).sum()) for v in ALT}
        chi = {v: float((c[1:4] - cp[v][1:4]) @ np.linalg.solve(C3, c[1:4] - cp[v][1:4])) for v in ALT}
        lev = float(np.mean(y - pred['law'][st])); levb = np.array([np.mean((y - pred['law'][st])[b]) for b in BI]).std()
        rat = (c[2] / c[1], c[3] / c[1]); ratb = np.array([(x[2] / x[1], x[3] / x[1]) for x in cb]).std(0)
        res[st] = dict(coef=c[1:].tolist(), se=se[1:].tolist(), pred={v: cp[v][1:].tolist() for v in ALT}, gamma=gam, gamma_se=float(gb.std()),
                       gamma_expected=gexp, chi2=chi, level_minus_law=lev, level_se=float(levb), level_pred={v: float(np.mean(pred[v][st] - pred['law'][st])) for v in ALT},
                       ratio_3_2=float(rat[0]), ratio_5_2=float(rat[1]), ratio_se=ratb.tolist(), mean=float(y.mean()))
        print('\n[%s] per-prime coefficients (y = stat - ignored; OLS on [1, I(2|D) .. I(13|D)]; bootstrap se, %d resamples)' % ('Var(s)' if st == 'var' else '<r~>', NBOOT))
        print('   p             ' + ' '.join('%9d' % p for p in PR))
        print('   measured      ' + ' '.join('%+9.5f' % x for x in c[1:]))
        print('   se            ' + ' '.join('%9.5f' % x for x in se[1:]))
        for v in ALT:
            print('   %-12s  ' % v + ' '.join('%+9.5f' % x for x in cp[v][1:]) + '   chi2(p=2,3,5) %8.1f   gamma expected %.3f' % (chi[v], gexp[v]))
        print('   gamma (amplitude on the law pattern, p = 2,3,5) = %.3f +- %.3f;  ratios c3/c2 %.3f, c5/c2 %.3f (se %.3f, %.3f)' % (gam, gb.std(), rat[0], rat[1], ratb[0], ratb[1]))
        print('   level: mean(stat - law) = %+.5f +- %.5f   (mean stat %.5f)' % (lev, levb, y.mean()))
    decide(res)
    if outfile: json.dump(res, open(outfile, 'w'), indent=1, default=float)
    return res


def decide(res):
    """the registered decisions (round_draft.txt), on one test"""
    V, R = res['var'], res['r']; z = lambda a, b, s: (a - b) / s
    print('\nDECISIONS (%s)' % res['test'])
    gv, sv = V['gamma'], V['gamma_se']; gr, sr = R['gamma'], R['gamma_se']
    print(' D1 Var amplitude: gamma_V = %.3f +- %.3f; law passes by measurement if |gamma_V - 1| <= 4 se: %s; within the theory band [0.85, 1.05]: %s'
          % (gv, sv, abs(gv - 1) <= 4 * sv, 0.85 <= gv <= 1.05))
    print('    slow-wave (dilation) comparator rejected if gamma_V < (1 + gamma_slow)/2 = %.3f and gamma_slow - gamma_V >= 5 se: %s'
          % ((1 + V['gamma_expected']['slow']) / 2, gv < (1 + V['gamma_expected']['slow']) / 2 and (V['gamma_expected']['slow'] - gv) >= 5 * sv))
    print('    bad primes ignored rejected if gamma_V >= 5 se: %s' % (gv >= 5 * sv))
    print(' D2 gap ratio: gamma_R = %.3f +- %.3f; "no change" (dilation / effective CUE size / kernel change / ignored) rejected if gamma_R >= 5 se and > 0.5: %s;'
          ' law passes by measurement if |gamma_R - 1| <= 4 se: %s; theory band [0.85, 1.25]: %s; negative sign (scout) rejected if gamma_R >= 5 se: %s'
          % (gr, sr, gr >= 5 * sr and gr > 0.5, abs(gr - 1) <= 4 * sr, 0.85 <= gr <= 1.25, gr >= 5 * sr))
    gb_ = V['gamma_expected']['bblm']
    print(' D3 effective-matrix-size law with the bad-prime constants (bblm) rejected if D2 rejects "no change", or if gamma_V < (1 + gamma_bblm)/2 = %.3f'
          ' and gamma_bblm - gamma_V >= 5 se: %s' % ((1 + gb_) / 2, (gr >= 5 * sr and gr > 0.5) or (gv < (1 + gb_) / 2 and gb_ - gv >= 5 * sv)))
    print(' D4 levels (secondary): Var - law %+.5f +- %.5f (expected +0.002..+0.004 from the formula\'s known overshoot; falsifier |.| > 0.008): %s;'
          ' r~ - law %+.5f +- %.5f (expected -0.001..-0.004; falsifier |.| > 0.008): %s'
          % (V['level_minus_law'], V['level_se'], abs(V['level_minus_law']) <= 0.008, R['level_minus_law'], R['level_se'], abs(R['level_minus_law']) <= 0.008))


def pool(files):
    rs = [json.load(open(f)) for f in files]
    print('POOLED over %s' % ', '.join('%s(L %.2f)' % (r['test'], r['meanL']) for r in rs))
    for st in ('var', 'r'):
        g = np.array([r[st]['gamma'] for r in rs]); s = np.array([r[st]['gamma_se'] for r in rs]); w = 1 / s ** 2
        gp = (w * g).sum() / w.sum(); sp = 1 / math.sqrt(w.sum()); chi = ((g - gp) ** 2 * w).sum()
        print('  %-3s gamma pooled %.3f +- %.3f; per test %s; heterogeneity chi2 %.1f on %d' % (st, gp, sp, ' '.join('%.3f' % x for x in g), chi, len(g) - 1))
    byname = {r['test']: r for r in rs}
    if 'B' in byname and 'Bp' in byname:
        for st, flat in (('var', 1.16), ('r', 1.31)):
            d = byname['Bp'][st]['gamma'] - byname['B'][st]['gamma']; sd = math.hypot(byname['Bp'][st]['gamma_se'], byname['B'][st]['gamma_se'])
            print('  L-dependence (same discriminants, B at heights 100-600 -> Bp at 600-1000), %s: gamma(Bp) - gamma(B) = %+.3f +- %.3f (law 0; no L-dependence about %+.2f)' % (st, d, sd, flat - 1))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('test', nargs='?'); ap.add_argument('--zeros'); ap.add_argument('--noise-only', action='store_true'); ap.add_argument('--no-completeness', action='store_true')
    ap.add_argument('--out'); ap.add_argument('--pool', nargs='+'); ap.add_argument('--make-synthetic'); ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    if a.pool: pool(a.pool); sys.exit()
    P = load_pred()
    if a.make_synthetic: make_synthetic(P, a.test, a.make_synthetic, a.seed, a.out); sys.exit()
    if a.no_completeness and not a.zeros: sys.exit('--no-completeness is for synthetic sets only')
    run(P, a.test, a.zeros, a.noise_only, a.out, a.no_completeness)
