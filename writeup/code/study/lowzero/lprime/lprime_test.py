"""lprime_test.py: the registered readings of log|L'(rho, chi_D)| at the zeros 100 <= gamma <= 600 (round_draft.txt).

usage:  lprime_test.py PRED.json SET:VALDIR:ZEROFILE[,ZEROFILE...] [SET:VALDIR:ZEROFILES ...] [--noise-only] [--boot 2000]
        SET in B R120 Dn Fn; VALDIR = output directory of lprime_values.py for that set; ZEROFILE(s) = the zero file(s) the
        values were computed from (used for the registered completeness rule). With several SETs the readings are made on
        the pooled L-functions with one nuisance intercept per set.
--noise-only: prints only noise (bootstrap se of every reading, residual scatter) and the power of each reading against its
        alternatives (|difference of predictions|/se); no fitted coefficient, level or comparison with a prediction is shown.
Self-contained: numpy, scipy.
"""
import sys, os, json, math
import numpy as np
from scipy import stats
from scipy.special import loggamma

PR = (2, 3, 5, 7, 11, 13); RD = (2, 3, 5)          # design primes; primes with a decision
LO, HI = 100.0, 600.0; YC = 0.2                      # band; truncation |L'| floor for the k = -1 reading
KS = (-0.5, 0.5)                                     # light-tailed moments read by OLS of log moments
WP = (2, 3, 5, 7)                                    # primes in the within-function reading


def theta(t, q, a): return np.imag(loggamma(0.25 + a / 2 + 0.5j * t)) + t / 2 * math.log(q / math.pi)


def complete(g, D):
    """Registered rule (round 138): |N - theta/pi| <= 1.5 on 200 heights of the band, largest unfolded gap <= 3.0."""
    q = abs(D); a = 0 if D > 0 else 1; t = np.linspace(LO, HI, 200)
    s = np.searchsorted(g, t) - (theta(t, q, a) - theta(0.0, q, a)) / math.pi
    b = g[(g >= LO) & (g <= HI)]; u = b / (2 * math.pi) * np.log(q * b / (2 * math.pi * math.e))
    return float(np.abs(s).max()), float(np.diff(u).max())


def kron(D, p):
    """Kronecker symbol (D/p) for a prime p."""
    if D % p == 0: return 0
    if p == 2: return 1 if D % 8 in (1, 7) else -1
    r = pow(D % p, (p - 1) // 2, p); return 1 if r == 1 else -1


def load(pred, spec):
    name, vdir, zf = spec.split(':'); Z = {}
    for f in zf.split(','):
        z = np.load(f); Z.update({k: np.sort(z[k]) for k in z.files})
    rows, excl = [], []
    for e in pred[name]:
        D = e['D']; k = 'D%d' % D; fn = os.path.join(vdir, k + '.npz')
        if not os.path.exists(fn): excl.append((D, 'no values')); continue
        v = np.load(fn); g = Z[k]; b = g[(g >= LO) & (g <= HI)]
        s, gap = complete(g, D)
        if s > 1.5 or gap > 3.0: excl.append((D, 'completeness %.2f %.2f' % (s, gap))); continue
        if len(v['gamma']) != len(b) or np.abs(v['gamma'] - b).max() > 1e-9: excl.append((D, 'values do not match the zeros')); continue
        chk = float((v['absL'] / np.exp(v['logabs'])).max())
        if chk > 1e-8: excl.append((D, 'not a zero: |L|/|L\'| = %.1e' % chk)); continue
        rows.append(dict(e, set=name, gamma=v['gamma'], logabs=v['logabs'], chk=chk))
    return rows, excl


def stats_per_lfunction(r):
    q = abs(r['D']); g = r['gamma']; Lj = np.log(q * g / (2 * math.pi)); y = r['logabs'] - np.log(Lj); m = y.mean(); c = y - m
    m2 = (c ** 2).mean(); out = {'mu': m, 'var': y.var(ddof=1), 'k3': (c ** 3).mean(), 'k4': (c ** 4).mean() - 3 * m2 ** 2, 'n': len(g)}
    a2 = np.exp(2 * r['logabs']); out['S1'] = float(a2.sum())
    inv = np.exp(-2 * r['logabs']); out['Sm1'] = float(inv.sum()); out['Sm1c'] = float(np.minimum(inv, YC ** -2).sum())
    out['n_trunc'] = int((r['logabs'] < math.log(YC)).sum())
    for k in KS: out['lm%+.1f' % k] = math.log(np.mean(np.exp(2 * k * r['logabs']) / Lj ** (k * (k + 2))))
    out['tail'] = float(np.mean(r['logabs'] - np.log(Lj) < math.log(0.05)))
    return out


def design(rows):
    sets = sorted(set(r['set'] for r in rows)); X = []
    for r in rows:
        x = [1.0] + [1.0 * (r['set'] == s) for s in sets[1:]] + [r['mean_lnL']] + [1.0 * (p in r['bad']) for p in PR]
        X.append(x)
    X = np.array(X); X[:, 1:] -= X[:, 1:].mean(0)        # centred nuisance and indicator columns
    return X, len(sets)


def ols(X, y):
    A = np.linalg.pinv(X.T @ X); b = A @ X.T @ y; res = y - X @ b; s2 = res @ res / (len(y) - X.shape[1])
    return b, np.sqrt(np.diag(A) * s2), res


def glm(X, S, P, it=60):
    """Quasi-Poisson ratio regression E S_i = P_i exp(X_i beta): solves sum X_i (S_i - P_i e^{X_i beta}) = 0."""
    b = np.zeros(X.shape[1]); b[0] = math.log(S.sum() / P.sum())
    for _ in range(it):
        mu = P * np.exp(X @ b); H = X.T @ (X * mu[:, None]); b = b + np.linalg.solve(H, X.T @ (S - mu))
    return b


def amp(fit, se, pref, palt):
    """Amplitude of the predicted contrast: sum w (fit - alt)(pref - alt) / sum w (pref - alt)^2 over the decision primes."""
    d = pref - palt; w = 1 / se ** 2; return float((w * (fit - palt) * d).sum() / (w * d * d).sum())


def within_parts(rows):
    """Per L-function normal equations of the pooled within-L-function regression of y = log|L'| - ln L on the
    Euler-factor variable A_p and the LDM term E_p (good p) and on cos, sin(gamma ln p) (bad p), one intercept each."""
    parts = []
    for r in rows:
        q = abs(r['D']); g = r['gamma']; Lj = np.log(q * g / (2 * math.pi)); y = r['logabs'] - np.log(Lj); cols = []
        for p in WP:
            chi = kron(r['D'], p); ph = g * math.log(p)
            if chi != 0:
                w = chi * p ** -0.5 * np.exp(-1j * ph); A = -np.log(np.abs(1 - w)); E = -2 * math.log(p) / Lj * np.real(w / (1 - w))
                cols += [A, E, 0 * g, 0 * g]
            else:
                cols += [0 * g, 0 * g, np.cos(ph), np.sin(ph)]
        M = np.array(cols).T; M = M - M.mean(0); yc = y - y.mean()
        amp2 = []                       # per-L-function squared amplitude of y at frequency ln p and at a control frequency
        for p in WP:
            for om in (math.log(p), math.log(p) + (0.1 if p < 7 else -0.1)):
                C = np.array([np.cos(om * g), np.sin(om * g)]).T; C = C - C.mean(0); cs = np.linalg.lstsq(C, yc, rcond=None)[0]
                amp2.append(float(cs @ cs))
        parts.append((M.T @ M, M.T @ yc, np.array(amp2), np.array([kron(r['D'], p) == 0 for p in WP])))
    return parts


def within(parts, idx):
    XtX = sum(parts[i][0] for i in idx); Xty = sum(parts[i][1] for i in idx)
    return np.linalg.lstsq(XtX, Xty, rcond=None)[0]


def amplitudes(parts, idx):
    """For each p in WP: mean squared amplitude at ln p over the L-functions with p | D (bad) and p not dividing D (good),
    and the noise floor = mean squared amplitude at the control frequency over the L-functions with p | D (no signal at or
    near ln p there). Returns bad/floor, good/floor."""
    A = np.array([parts[i][2] for i in idx]); Bd = np.array([parts[i][3] for i in idx]); out = []
    for j, p in enumerate(WP):
        a, c = A[:, 2 * j], A[:, 2 * j + 1]; fl = c[Bd[:, j]].mean()
        out += [a[Bd[:, j]].mean() / fl, a[~Bd[:, j]].mean() / fl]
    return np.array(out)


def readings(rows, st, idx, pred, wp):
    R = {}; X, ns = design([rows[i] for i in idx]); ip = ns + 1      # columns of the prime indicators start at ip
    sub = lambda key: np.array([st[i][key] for i in idx]); pr = lambda key: np.array([rows[i][key] for i in idx])
    H = np.linalg.pinv(X.T @ X) @ X.T
    match = lambda v: (H @ v)[ip:ip + len(PR)]
    for stat, key in (('mean', 'cum1'), ('var', 'cum2'), ('k3', 'cum3'), ('k4', 'cum4')):
        y = sub({'mean': 'mu'}.get(stat, stat)); b, se, res = ols(X, y)
        R[stat] = {'fit': b[ip:ip + len(PR)], 'se_ols': se[ip:ip + len(PR)], 'res_sd': float(res.std()),
                   'ldm': match(pr(key + '_ldm')), 'lead': match(pr(key + '_lead')), 'none': match(pr(key + '_none'))}
    for k in KS:
        y = sub('lm%+.1f' % k); b, se, res = ols(X, y)
        R['m%+.1f' % k] = {'fit': b[ip:ip + len(PR)], 'se_ols': se[ip:ip + len(PR)], 'res_sd': float(res.std()),
                           'ldm': match(pr('logmom%+.1f_ldm' % k)), 'lead': match(pr('logmom%+.1f_lead' % k)), 'none': 0 * match(pr('logmom%+.1f_ldm' % k))}
    # k = 1, exact recipe (absolute): beta = 0 predicted, including the level
    S1 = sub('S1'); P1 = pr('S1_recipe'); b = glm(X, S1, P1)
    R['S1'] = {'level': float(b[0]), 'fit': b[ip:ip + len(PR)], 'pooled_ratio': float(S1.sum() / P1.sum()),
               'none': (H @ np.log(pr('S1_none') / P1))[ip:ip + len(PR)], 'lead': (H @ np.log(pr('S1_lead') / P1))[ip:ip + len(PR)],
               'lead_level': float(np.log(pr('S1_lead').sum() / P1.sum())), 'ldm': np.zeros(len(PR))}
    # k = -1 (Pearce-Crump), truncated at |L'| >= YC: tail term 2 C_q y_c per zero (leading-order tail theory)
    Pc = pr('Sm1_conj') - 2 * YC * pr('tail_C_weighted'); Pn = pr('Sm1_none') - 2 * YC * pr('tail_C_weighted') * np.exp(-pr('tail_logfactor_lead'))
    Sc = sub('Sm1c'); b = glm(X, Sc, Pc)
    R['Sm1'] = {'level': float(b[0]), 'fit': b[ip:ip + len(PR)], 'pooled_ratio_trunc': float(Sc.sum() / Pc.sum()),
                'pooled_ratio_untrunc': float(sub('Sm1').sum() / pr('Sm1_conj').sum()), 'none': (H @ np.log(Pn / Pc))[ip:ip + len(PR)],
                'ldm': np.zeros(len(PR)), 'lead': np.zeros(len(PR)), 'n_trunc': int(sub('n_trunc').sum())}
    R['tail'] = {'fit': ols(X, sub('tail'))[0][ip:ip + len(PR)]}
    R['within'] = within(wp, idx); R['amp'] = amplitudes(wp, idx)
    return R


def main(argv):
    noise = '--noise-only' in argv; nb = int(argv[argv.index('--boot') + 1]) if '--boot' in argv else 2000
    pred = json.load(open(argv[1])); specs = [a for a in argv[2:] if ':' in a]
    rows, excl = [], []
    for s in specs:
        r, e = load(pred, s); rows += r; excl += e
    st = [stats_per_lfunction(r) for r in rows]; n = len(rows); idx = np.arange(n)
    print('sets %s: %d L-functions kept; excluded: %s' % ([s.split(':')[0] for s in specs], n, excl if excl else 'none'))
    print('zeros used %d; max |L(rho)|/|L\'(rho)| %.1e; groups (2,3,5,7,11,13): %s' % (sum(x['n'] for x in st), max(r['chk'] for r in rows), [sum(p in r['bad'] for r in rows) for p in PR]))
    wp = within_parts(rows); R = readings(rows, st, idx, pred, wp); rng = np.random.default_rng(20261125)
    B = [readings(rows, st, rng.integers(0, n, n), pred, wp) for _ in range(nb)]
    def bse(key, sub='fit'):
        return np.std([b[key][sub] for b in B], axis=0)
    out = {}
    print('\nreading            p:      ' + '      '.join('%-6d' % p for p in PR))
    for key in ('mean', 'var', 'k3', 'k4', 'm-0.5', 'm+0.5', 'S1', 'Sm1'):
        se = bse(key); out[key] = {'se': se}
        print('%-6s bootstrap se        ' % key + ' '.join('%.5f' % v for v in se))
        if not noise:
            print('%-6s fit                 ' % key + ' '.join('%+.5f' % v for v in R[key]['fit']))
        for h in ('ldm', 'lead', 'none'):
            if h in R[key]:
                print('%-6s pred %-4s (z vs fit)' % (key, h) + ' '.join('%+.5f' % v for v in R[key][h]) +
                      ('' if noise else '   z: ' + ' '.join('%+.1f' % v for v in (R[key]['fit'] - R[key][h]) / se)))
    # amplitude readings over the decision primes
    J = [PR.index(p) for p in RD]; print(); AM = {}
    for key, ref, alts in (('mean', 'ldm', ('none',)), ('var', 'ldm', ('lead', 'none')), ('S1', 'ldm', ('lead', 'none')), ('Sm1', 'ldm', ('none',)),
                           ('m-0.5', 'ldm', ('lead', 'none')), ('m+0.5', 'ldm', ('lead', 'none')), ('k3', 'ldm', ('lead', 'none')), ('k4', 'ldm', ('lead', 'none'))):
        se = out[key]['se'][J]
        for alt in alts:
            pr_, pa = R[key][ref][J], R[key][alt][J]
            if np.allclose(pr_, pa): continue
            g = amp(R[key]['fit'][J], se, pr_, pa); gs = float(np.std([amp(b[key]['fit'][J], se, b[key][ref][J], b[key][alt][J]) for b in B]))
            AM[(key, alt)] = (g, gs)
            # power: under H_ref the amplitude is 1, under H_alt 0 -> separation 1/gs
            print('amplitude %-6s %s vs %-4s: se %.4f, separation %.1f se' % (key, ref, alt, gs, 1 / gs) + ('' if noise else ', reading %.4f' % g))
    for key in ('S1', 'Sm1'):
        lv = float(np.std([b[key]['level'] for b in B]))
        print('level %-4s (log of data/prediction, centred design): se %.4f' % (key, lv) + ('' if noise else ', reading %+.4f, pooled ratio %s' % (R[key]['level'], R[key].get('pooled_ratio', R[key].get('pooled_ratio_trunc')))))
    if 'S1' in R: print('S1 leading-order-only level alternative: log ratio %+.3f (%.0f se)' % (R['S1']['lead_level'], abs(R['S1']['lead_level']) / float(np.std([b['S1']['level'] for b in B]))))
    print('k = -1 truncation: %d zeros below |L\'| = %.2f' % (R['Sm1']['n_trunc'], YC))
    big = max((float(np.exp(-2 * r['logabs']).max()), r['D']) for r in rows); tot = sum(x['Sm1'] for x in st)
    print('k = -1 largest single term |L\'|^-2 = %.1f (D = %d), %.2f%% of the pooled sum; min |L\'| = %.3f' % (big[0], big[1], 100 * big[0] / tot, big[0] ** -0.5))
    W = np.array([b['within'] for b in B]); wse = W.std(0); lab = []
    for p in WP: lab += ['A_%d' % p, 'E_%d' % p, 'cos_%d' % p, 'sin_%d' % p]
    print('\nwithin-function reading (prediction: A_p 1, E_p 1 [LDM; 0 if P-only], cos/sin of a bad prime 0; orbit kept: amplitude p^-1/2):')
    for j, l in enumerate(lab):
        print('  %-6s se %.4f' % (l, wse[j]) + ('' if noise else '  fit %+.4f  z(vs pred) %+.1f' % (R['within'][j], (R['within'][j] - (1 if l[0] in 'AE' else 0)) / wse[j])))
    Am = np.array([b['amp'] for b in B]); ase = Am.std(0)
    print('\nphase-amplitude reading: mean squared amplitude of y at frequency ln p / noise floor (control frequency);')
    print('  prediction: bad p 1 (no orbit); good p 1 + p^-1 (1 - 2 ln p/L)^2 / floor; orbit kept: bad p as good p')
    for j, p in enumerate(WP):
        print('  p=%d  bad/floor se %.3f   good/floor se %.3f' % (p, ase[2 * j], ase[2 * j + 1]) + ('' if noise else '   readings %.3f  %.3f' % (R['amp'][2 * j], R['amp'][2 * j + 1])))
    if not noise:
        print('\nDECISIONS (round_draft.txt; se = bootstrap se on these data)')
        f1, s1 = R['S1']['fit'][J], out['S1']['se'][J]; lv1 = float(np.std([b['S1']['level'] for b in B]))
        print('  D1 k = 1 recipe: every |beta_p| <= 4 se + 0.01 and |level| <= 4 se + 0.02: %s;  no-conductor-prime alternative rejected (gamma >= 5 se, >= 0.5): %s' % (
            bool(np.all(np.abs(f1) <= 4 * s1 + 0.01) and abs(R['S1']['level']) <= 4 * lv1 + 0.02), bool(AM[('S1', 'none')][0] >= max(5 * AM[('S1', 'none')][1], 0.5))))
        g, gs = AM[('Sm1', 'none')]
        print('  D2 Pearce-Crump factor p/(p+1): gamma_PC = %.3f +- %.3f; passes (|gamma - 1| <= 4 se + 0.15): %s; no factor rejected (gamma >= 5 se, >= 0.5): %s; falsified (gamma < 0.5 or > 1.5 by 4 se): %s' % (
            g, gs, abs(g - 1) <= 4 * gs + 0.15, g >= max(5 * gs, 0.5), (g + 4 * gs < 0.5) or (g - 4 * gs > 1.5)))
        g, gs = AM[('mean', 'none')]
        print('  D3 mean shift ln p ln(p/(p-1))/L: gamma_M = %.3f +- %.3f; passes (|gamma - 1| <= 4 se + 0.2): %s; no shift rejected (>= 5 se, >= 0.5): %s; falsified (gamma < 0.5 or > 1.5 by 4 se): %s' % (
            g, gs, abs(g - 1) <= 4 * gs + 0.2, g >= max(5 * gs, 0.5), (g + 4 * gs < 0.5) or (g - 4 * gs > 1.5)))
        gv, gvs = AM[('var', 'none')]; g1, g1s = AM[('var', 'lead')]
        print('  D4 variance: gamma_V = %.3f +- %.3f (no-dependence rejected: %s; missing-factor law falsified if gamma_V + 4 se < 0.5: %s); gamma_1 = %.3f +- %.3f: LDM passes (|gamma_1 - 1| <= 4 se + 0.35): %s; leading-order-only rejected (gamma_1 >= 5 se, >= 0.5): %s' % (
            gv, gvs, gv >= max(5 * gvs, 0.5), gv + 4 * gvs < 0.5, g1, g1s, abs(g1 - 1) <= 4 * g1s + 0.35, g1 >= max(5 * g1s, 0.5)))
        w = R['within']; ok_A = all(abs(w[lab.index('A_%d' % p)] - 1) <= 4 * wse[lab.index('A_%d' % p)] + 0.1 for p in (2, 3))
        ok_E = all(abs(w[lab.index('E_%d' % p)] - 1) <= 4 * wse[lab.index('E_%d' % p)] + 0.35 for p in (2, 3))
        rej_P = all(w[lab.index('E_%d' % p)] >= 5 * wse[lab.index('E_%d' % p)] for p in (2, 3))
        ok_ph = all(abs(w[lab.index('%s_%d' % (c, p))]) <= 4 * wse[lab.index('%s_%d' % (c, p))] for p in (2, 3, 5) for c in ('cos', 'sin'))
        rej_kept = all(R['amp'][2 * j] <= 1 + 4 * ase[2 * j] + 0.1 and R['amp'][2 * j + 1] >= 10 for j, p in enumerate(WP) if p in (2, 3))
        print('  D5 within-function: A_p = 1 (p = 2, 3): %s; E_p = 1 (LDM): %s; P-only rejected: %s; no phase dependence for bad p = 2, 3, 5: %s; orbit kept rejected: %s' % (ok_A, ok_E, rej_P, ok_ph, rej_kept))


if __name__ == '__main__':
    if len(sys.argv) < 3: sys.exit(__doc__)
    main(sys.argv)
