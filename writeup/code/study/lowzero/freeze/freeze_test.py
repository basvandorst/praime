"""freeze_test.py: the registered test (self-contained; numpy only).
Per-prime readings of the maximum and the partition function of |L(1/2+it, chi_D)| on intervals of length 2 pi,
real primitive chi_D, against frozen predictions (predictions_freeze.json) and against the two surgeries.

usage: freeze_test.py VALUES_DIR SETS.json KEY PREDICTIONS.json OUT.json [--sha SHA256] [--boot 2000] [--label TXT]
                     [--model]   (--model: synthetic values, skips the mean-square completeness check)
  VALUES_DIR holds absL_<D>.npy (float64, 18001 values of |L(1/2 + i t)| at t = 100 + k/20, k = 0..18000), as written
  by moment_values.py; if VALUES_DIR/manifest.json exists, every file's sha256 must match it.
  SETS.json[KEY] is the list of discriminants.  --sha: required sha256 of PREDICTIONS.json (as committed).

Definitions (all fixed before any test statistic is computed):
 intervals  I_j = [100 + 2 pi j, 100 + 2 pi (j+1)), j = 0..142.
 per interval: M = max ln|L| (grid maximum refined by the parabola through it and its neighbours, vertex kept only
   inside I_j; or the value at an exact end of I_j, quadratic interpolation through the 3 nearest grid points); Z_b = (1/2 pi) int_{I_j} |L|^(2b) dt (trapezoid, linear interpolation in the two partial cells);
   mean_ln = (1/2 pi) int ln max(|L|, 1e-12).
 per character: means over j of M, ln Z_b (b = 0.5, 1, 1.5, 2, 3), mean_ln;  slope32 = mean(ln Z_3 - ln Z_2);
   varM = variance (ddof 1) of M over j;  lnEZ_b = ln mean_j Z_b (b = 1, 2; annealed).
 regression: OLS of each statistic on [1, ln q - mean, I(p|D) for p = 2, 3, 5, 7, 11, 13]; c_p = coefficient of
   I(p|D) = 'p missing minus p present'.
 surgery (E_p(t) = 1 - chi(p) p^(-1/2-it)):
   removal  r_p = mean over characters with p not | D of [stat(|L| |E_p|) - stat(|L|)];
   addition a_p = mean over characters with p | D of [stat(|L|) - stat(|L| / |1 - e p^(-1/2-it)|)], e = +1, -1 averaged.
 bootstrap: B resamples of characters (numpy default_rng(149150)); every reading is recomputed on each resample;
   se = bootstrap standard deviation.
Readings: see the registration text; the code below computes C1-C4, T1-T4 and the decisions.
"""
import sys, os, json, hashlib
import numpy as np

T0, H, NPT = 100.0, 0.05, 18001
THETA = 2 * np.pi
NINT = 143
BETAS = (0.5, 1.0, 1.5, 2.0, 3.0)
PR = (2, 3, 5, 7)
REGP = (2, 3, 5, 7, 11, 13)
EG = 0.5772156649015329
STATS = ['M', 'lnZ_0.5', 'lnZ_1', 'lnZ_1.5', 'lnZ_2', 'lnZ_3', 'slope32', 'mean_ln', 'varM', 'lnEZ_1', 'lnEZ_2']
REL = ['M', 'lnZ_3', 'varM', 'lnZ_1.5', 'lnZ_2', 'slope32', 'lnZ_0.5', 'lnZ_1', 'lnEZ_2']
# registered theory bands (the theory's own uncertainty; NOT measurement bands)
BAND = {'crit_M': (0.92, 1.08), 'crit_lnZ_3': (0.85, 1.10), 'frozen_slope32': (0.90, 1.20), 'W': (-0.15, 0.15)}


# ---------------------------------------------------------------- observables
def interval_weights(a, b):
    ka = (a - T0) / H; kb = (b - T0) / H
    i0 = int(np.floor(ka)); i1 = int(np.floor(kb))
    idx = np.arange(i0, i1 + 2); w = np.zeros(len(idx))
    for c in range(i0, i1 + 1):
        lo = max(ka, c); hi = min(kb, c + 1)
        if hi <= lo:
            continue
        s0, s1 = lo - c, hi - c
        w[c - i0] += ((s1 - s0) - (s1 ** 2 - s0 ** 2) / 2) * H
        w[c + 1 - i0] += ((s1 ** 2 - s0 ** 2) / 2) * H
    return idx, w


WTS = [interval_weights(T0 + j * THETA, T0 + (j + 1) * THETA) for j in range(NINT)]


def observables(absL):
    lnL = np.log(np.maximum(absL, 1e-12))
    out = {'M': np.zeros(NINT), 'mean_ln': np.zeros(NINT)}
    for b in BETAS:
        out['lnZ_%g' % b] = np.zeros(NINT)
    for j, (idx, w) in enumerate(WTS):
        a = T0 + j * THETA; e = a + THETA
        ic = np.clip(idx, 0, NPT - 1); vals = absL[ic]
        for b in BETAS:
            out['lnZ_%g' % b][j] = np.log(np.dot(w, vals ** (2 * b)) / THETA)
        out['mean_ln'][j] = np.dot(w, lnL[ic]) / THETA
        k0 = int(np.ceil((a - T0) / H - 1e-9)); k1 = min(int(np.floor((e - T0) / H + 1e-9)), NPT - 1)
        km = k0 + int(np.argmax(lnL[k0:k1 + 1]))
        v = lnL[km]
        if 1 <= km <= NPT - 2:
            y0, y1, y2 = lnL[km - 1], lnL[km], lnL[km + 1]
            den = y0 - 2 * y1 + y2
            if den < 0:
                d = max(-1.0, min(1.0, 0.5 * (y0 - y2) / den))
                tv = T0 + (km + d) * H
                if a <= tv <= e:
                    v = y1 - 0.25 * (y0 - y2) * d
        for x in (a, e):          # exact interval ends: quadratic interpolation through the 3 nearest grid points
            u = (x - T0) / H; kc = int(min(max(round(u), 1), NPT - 2)); s = u - kc
            ve = lnL[kc] + 0.5 * s * (lnL[kc + 1] - lnL[kc - 1]) + 0.5 * s * s * (lnL[kc + 1] - 2 * lnL[kc] + lnL[kc - 1])
            v = max(v, ve)
        out['M'][j] = v
    return out


def char_stats(absL):
    o = observables(absL)
    s = {k: float(o[k].mean()) for k in ['M', 'lnZ_0.5', 'lnZ_1', 'lnZ_1.5', 'lnZ_2', 'lnZ_3', 'mean_ln']}
    s['slope32'] = s['lnZ_3'] - s['lnZ_2']
    s['varM'] = float(o['M'].var(ddof=1))
    s['lnEZ_1'] = float(np.log(np.mean(np.exp(o['lnZ_1']))))
    s['lnEZ_2'] = float(np.log(np.mean(np.exp(o['lnZ_2']))))
    return s


def kron(D, p):
    if D % p == 0:
        return 0
    if p == 2:
        return 1 if D % 8 in (1, 7) else -1
    return 1 if pow(D % p, (p - 1) // 2, p) == 1 else -1


def primes_of(n):
    n = abs(n); out = []; d = 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def twisted_ratio(D, p):
    """theorem (twisted second moment): mean(|L E_p|^2)/mean(|L|^2) over [100, 1000], main terms."""
    t = T0 + H * np.arange(NPT)
    Lp = np.log(abs(D) * t / (2 * np.pi)) + 2 * EG + 2 * sum(np.log(l) / (l - 1) for l in primes_of(D))
    f = (1 + 1 / p) * Lp - (2 / p) * (Lp - np.log(p))
    return float(np.trapezoid(f, t) / np.trapezoid(Lp, t))


def second_moment_ratio(D, absL):
    t = T0 + H * np.arange(NPT)
    q = abs(D)
    main = (np.prod([1 - 1 / l for l in primes_of(q)]) *
            (np.log(q * t / (2 * np.pi)) + 2 * EG + 2 * sum(np.log(l) / (l - 1) for l in primes_of(q))))
    return float(np.trapezoid(absL ** 2, t) / np.trapezoid(main, t))


def per_character(D, absL):
    g = T0 + H * np.arange(NPT)
    s0 = char_stats(absL)
    rec = {'D': D, 'stats': s0, 'surg': {}, 'tw': {}, 'm2ratio': second_moment_ratio(D, absL)}
    for p in PR:
        c = kron(D, p)
        if c:
            s = char_stats(absL * np.abs(1 - c * p ** -0.5 * np.exp(-1j * g * np.log(p))))
            rec['surg'][str(p)] = {'kind': 'rem', **{k: s[k] - s0[k] for k in STATS}}
            rec['tw'][str(p)] = twisted_ratio(D, p)
        else:
            acc = {k: 0.0 for k in STATS}
            for e in (1, -1):
                s = char_stats(absL / np.abs(1 - e * p ** -0.5 * np.exp(-1j * g * np.log(p))))
                for k in STATS:
                    acc[k] += 0.5 * (s0[k] - s[k])
            rec['surg'][str(p)] = {'kind': 'add', **acc}
    return rec


# ---------------------------------------------------------------- readings
def design(Ds):
    lnq = np.log(np.abs(np.array(Ds, float)))
    return np.column_stack([np.ones(len(Ds)), lnq - lnq.mean()] +
                           [[1.0 if D % p == 0 else 0.0 for D in Ds] for p in REGP])


def core(recs, idx):
    """c, r, a per statistic and prime on the resample idx; plus control quantities."""
    R = [recs[i] for i in idx]
    Ds = [r['D'] for r in R]
    X = design(Ds)
    XtXi = np.linalg.pinv(X.T @ X)
    out = {}
    for s in STATS:
        y = np.array([r['stats'][s] for r in R])
        b = XtXi @ (X.T @ y)
        c = b[2:6]
        rr, aa = [], []
        for p in PR:
            vr = [r['surg'][str(p)][s] for r in R if r['surg'][str(p)]['kind'] == 'rem']
            va = [r['surg'][str(p)][s] for r in R if r['surg'][str(p)]['kind'] == 'add']
            rr.append(np.mean(vr) if vr else np.nan); aa.append(np.mean(va) if va else np.nan)
        out[s] = (c, np.array(rr), np.array(aa))
    ctl = {}
    for p in PR:
        tr = [np.exp(r['surg'][str(p)]['lnEZ_1']) / r['tw'][str(p)] for r in R if r['surg'][str(p)]['kind'] == 'rem']
        ad = [np.exp(-r['surg'][str(p)]['lnEZ_1']) * (1 - 1 / p) for r in R if r['surg'][str(p)]['kind'] == 'add']
        ctl[p] = (np.mean(tr), np.mean(ad))
    return out, ctl


def amp(c, l, w):
    return float(np.sum(w * c * l) / np.sum(w * l * l))


def readings(core_out, ctl, pred, W):
    """all scalar readings from one (re)sample; W = dict of weight vectors (fixed from the bootstrap)."""
    rd = {}
    co = pred['coef']
    lc = lambda h, s: np.array([co[h][s][str(p)] for p in PR])
    # C1: annealed mean square against the theorem's per-prime factor
    rd['C1_A_EZ1'] = amp(core_out['lnEZ_1'][0], lc('crit', 'lnEZ_1'), W['c_lnEZ_1'])
    # C2, C4: removal ratio / twisted value; addition ratio x (1 - 1/p)
    for p in PR:
        rd['C2_tw_%d' % p] = float(ctl[p][0]); rd['C4_add_%d' % p] = float(ctl[p][1])
    # C3: mean of ln|L|
    for k, p in enumerate(PR):
        rd['C3_mean_ln_%d' % p] = float(core_out['mean_ln'][0][k])
    # T1, T2: relational gamma = sum w d e / sum w e^2, d = c - r, e = a - r
    for s in REL:
        c, r, a = core_out[s]
        d, e = c - r, a - r
        rd['G_%s' % s] = float(np.sum(W['d_' + s] * d * e) / np.sum(W['d_' + s] * e * e))
        for k, p in enumerate(PR):
            rd['c_%s_%d' % (s, p)] = float(c[k]); rd['r_%s_%d' % (s, p)] = float(r[k])
            rd['a_%s_%d' % (s, p)] = float(a[k])
    # T3: absolute amplitudes on the crit pattern; T4: frozen slope
    for s in ('M', 'lnZ_3', 'lnZ_1.5', 'lnZ_2', 'slope32'):
        rd['A_crit_%s' % s] = amp(core_out[s][0], lc('crit', s), W['c_' + s])
        rd['A_crit_rem_%s' % s] = amp(core_out[s][1], lc('crit', s), W['c_' + s])
    return rd


def expected_amplitudes(pred, W):
    co = pred['coef']; e = {}
    for s in ('M', 'lnZ_3', 'lnZ_1.5', 'lnZ_2', 'slope32', 'lnEZ_1'):
        l = np.array([co['crit'][s][str(p)] for p in PR])
        for h in ('lead', 'naive', 'spread', 'zero'):
            if co[h][s] is None:
                continue
            m = np.array([co[h][s][str(p)] for p in PR])
            e['%s_%s' % (h, s)] = amp(m, l, W['c_' + s] if 'c_' + s in W else np.ones(4))
    return e


def main(argv):
    argv = list(argv)
    opt = {}
    for key in ('--sha', '--boot', '--label'):
        if key in argv:
            i = argv.index(key); opt[key] = argv[i + 1]; del argv[i:i + 2]
    model = '--model' in argv          # dry runs on the random model: no mean-square completeness check
    if model:
        argv.remove('--model')
    vdir, setsf, key, predf, outf = argv[1:6]
    B = int(opt.get('--boot', 2000))
    raw = open(predf, 'rb').read()
    sha = hashlib.sha256(raw).hexdigest()
    if '--sha' in opt and opt['--sha'] != sha:
        sys.exit('predictions sha256 mismatch: %s' % sha)
    pred = json.loads(raw)
    Ds = json.load(open(setsf))[key]
    if [int(x) for x in pred['test_set']] != [int(x) for x in Ds] and opt.get('--label', '') != 'dryrun':
        sys.exit('set differs from the predictions file')
    man = None
    if os.path.exists(os.path.join(vdir, 'manifest.json')):
        man = json.load(open(os.path.join(vdir, 'manifest.json')))['files']
    # completeness
    recs, excluded = [], []
    for D in Ds:
        f = os.path.join(vdir, 'absL_%d.npy' % D)
        why = None
        if not os.path.exists(f):
            why = 'missing'
        else:
            if man is not None and man.get('absL_%d.npy' % D) != hashlib.sha256(open(f, 'rb').read()).hexdigest():
                why = 'sha256 differs from manifest'
            v = np.load(f)
            if why is None and (v.shape != (NPT,) or not np.all(np.isfinite(v)) or v.min() < 0 or v.max() > 1e6):
                why = 'bad values'
        if why is None:
            rec = per_character(D, v)
            if not model and not (0.85 <= rec['m2ratio'] <= 1.15):
                why = 'mean square / theorem = %.3f outside [0.85, 1.15]' % rec['m2ratio']
        if why:
            excluded.append([D, why]); continue
        recs.append(rec)
    n = len(recs)
    full, ctl = core(recs, np.arange(n))
    rng = np.random.default_rng(149150)
    boots = [core(recs, rng.integers(0, n, n)) for _ in range(B)]
    # weights from the bootstrap variances
    W = {}
    for s in STATS:
        cb = np.array([b[0][s][0] for b in boots]); rb = np.array([b[0][s][1] for b in boots])
        W['c_' + s] = 1 / np.var(cb, axis=0)
        W['d_' + s] = 1 / np.var(cb - rb, axis=0)
    rd = readings(full, ctl, pred, W)
    bd = [readings(b[0], b[1], pred, W) for b in boots]
    se = {k: float(np.std([x[k] for x in bd])) for k in rd}
    ex = expected_amplitudes(pred, W)
    dec = decisions(rd, se, ex)
    res = {'label': opt.get('--label', ''), 'n': n, 'excluded': excluded, 'predictions_sha256': sha,
           'groups': {p: sum(1 for r in recs if r['D'] % p == 0) for p in REGP},
           'readings': rd, 'se': se, 'expected': ex, 'decisions': dec,
           'm2ratio': [r['m2ratio'] for r in recs], 'chars': recs}
    json.dump(res, open(outf, 'w'), indent=1)
    report(res)


def decisions(rd, se, ex):
    d = {}
    # controls
    d['C1'] = 'pass' if abs(rd['C1_A_EZ1'] - 1) <= 4 * se['C1_A_EZ1'] else 'FAIL'
    d['C2'] = 'pass' if all(abs(rd['C2_tw_%d' % p] - 1) <= max(4 * se['C2_tw_%d' % p], 0.01) for p in PR) else 'FAIL'
    d['C3'] = 'pass' if all(abs(rd['C3_mean_ln_%d' % p]) <= 4 * se['C3_mean_ln_%d' % p] for p in PR) else 'FAIL'
    d['C4'] = 'pass' if all(abs(rd['C4_add_%d' % p] - 1) <= max(4 * se['C4_add_%d' % p], 0.01) for p in PR) else 'FAIL'
    # T1/T2 relational: H_W gamma = 0, H_G gamma = 1
    for s in ('M', 'lnZ_3', 'varM', 'lnZ_1.5', 'lnZ_2'):
        g, e = rd['G_' + s], se['G_' + s]
        w = 'pass' if abs(g) <= 4 * e else ('consistent within theory band' if BAND['W'][0] <= g <= BAND['W'][1]
                                             else 'H_W FALSIFIED' if abs(g) > 5 * e else 'undecided')
        gg = 'H_G rejected' if (1 - g) >= 5 * e and g < 0.5 else ('H_G pass' if abs(g - 1) <= 4 * e else 'H_G not rejected')
        d['T_rel_' + s] = '%s; %s' % ('H_W ' + w if not w.startswith('H_W') else w, gg)
    # T3 absolute on the crit pattern
    for s in ('M', 'lnZ_3'):
        A, e = rd['A_crit_' + s], se['A_crit_' + s]
        lo, hi = BAND['crit_' + s]
        crit = 'pass' if abs(A - 1) <= 4 * e else ('consistent at the theory accuracy' if lo <= A <= hi else
                                                    'FALSIFIED' if min(abs(A - lo), abs(A - hi)) > 4 * e else 'undecided')
        El = ex['lead_' + s]
        lead = 'rejected' if (El - A) >= 5 * e and A < (1 + El) / 2 else 'not rejected'
        z = 'rejected' if A >= 5 * e else 'not rejected'
        d['T3_' + s] = 'H_crit %s; H_lead (%.3f) %s; H_0 %s' % (crit, El, lead, z)
    # T4 freezing: slope32 on the crit pattern
    A, e = rd['A_crit_slope32'], se['A_crit_slope32']
    lo, hi = BAND['frozen_slope32']
    fr = 'pass' if (lo - 4 * e) <= A <= (hi + 4 * e) else 'FALSIFIED'
    out = []
    for h in ('naive', 'spread'):
        Eh = ex['%s_slope32' % h]
        out.append('H_%s (%.2f) %s' % (h, Eh, 'rejected' if (Eh - A) >= 5 * e and A < (hi + Eh) / 2 else 'not rejected'))
    d['T4_slope32'] = 'frozen %s; %s' % (fr, '; '.join(out))
    return d


def report(res):
    rd, se = res['readings'], res['se']
    print('%s  n = %d  excluded %d %s  groups %s' % (res['label'], res['n'], len(res['excluded']), res['excluded'][:5],
                                                     res['groups']))
    print('mean square / theorem: mean %.4f sd %.4f' % (np.mean(res['m2ratio']), np.std(res['m2ratio'])))
    print('C1 annealed mean-square amplitude on ln R_p: %.4f +- %.4f  (leading-only factor would give %.3f)' % (
        rd['C1_A_EZ1'], se['C1_A_EZ1'], res['expected'].get('lead_lnEZ_1', float('nan'))))
    print('C2 removal / twisted: ' + ' '.join('%d: %.4f+-%.4f' % (p, rd['C2_tw_%d' % p], se['C2_tw_%d' % p]) for p in PR))
    print('C4 addition x (1-1/p): ' + ' '.join('%d: %.4f+-%.4f' % (p, rd['C4_add_%d' % p], se['C4_add_%d' % p]) for p in PR))
    print('C3 mean ln|L| coef: ' + ' '.join('%d: %+.4f+-%.4f' % (p, rd['C3_mean_ln_%d' % p], se['C3_mean_ln_%d' % p]) for p in PR))
    print('\nper prime (p = 2, 3, 5, 7): regression c | removal r | addition a')
    for s in REL:
        print('%-8s c %s | r %s | a %s | gamma %.3f +- %.3f' % (
            s, ' '.join('%+.3f(%.3f)' % (rd['c_%s_%d' % (s, p)], se['c_%s_%d' % (s, p)]) for p in PR),
            ' '.join('%+.3f' % rd['r_%s_%d' % (s, p)] for p in PR),
            ' '.join('%+.3f' % rd['a_%s_%d' % (s, p)] for p in PR), rd['G_' + s], se['G_' + s]))
    print('\namplitudes on the crit pattern (data | removal):')
    for s in ('M', 'lnZ_1.5', 'lnZ_2', 'lnZ_3', 'slope32'):
        print('  %-8s %.4f +- %.4f | %.4f +- %.4f   expected: %s' % (
            s, rd['A_crit_' + s], se['A_crit_' + s], rd['A_crit_rem_' + s], se['A_crit_rem_' + s],
            ' '.join('%s %.3f' % (k.split('_')[0], v) for k, v in res['expected'].items() if k.endswith('_' + s))))
    print('\ndecisions:')
    for k, v in res['decisions'].items():
        print('  %-14s %s' % (k, v))


if __name__ == '__main__':
    main(sys.argv)
