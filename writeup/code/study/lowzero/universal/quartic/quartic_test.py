"""quartic_test.py: frozen blind test of the murmuration density of S4 quartic fields (degree-3 Artin L-functions zeta_K/zeta).

  counts --range disc|hold --r2 0|1|2 --budget SEC   field list (PARI nflist) and per-field bin sums; resumable (data/<range>_r<r2>.npz)
  read   --range disc|hold                            statistics and decisions; writes result_<range>.json
Ranges: disc = 5000 <= |D| < 20000 (discovery), hold = 20000 <= |D| < 60000 (holdout). Primes p <= 8|D|, p not dividing D.
Statistic: coarse bins in u = (p/|D|)^(1/3), edges 0, 1/24, 3/24, ..., 47/24; m_j = sum ln p sqrt|D| a_K(p) / sum ln p over fields and primes;
errors: delete-one jackknife over fields (full covariance). Refuses to run (read; counts on the holdout) unless predictions_quartic.json has the
registered sha256."""
import sys, os, json, time, math, hashlib, argparse
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
PRED_FILE = os.path.join(HERE, 'predictions_quartic.json')
EXPECTED_SHA = '148c1231be5767a32a2f4c71d06ad219f7f9e5ad5dfe276e2773263e1ace19c1'
RANGES = {'disc': (5000, 20000), 'hold': (20000, 60000)}

def load_pred():
    raw = open(PRED_FILE, 'rb').read(); sha = hashlib.sha256(raw).hexdigest()
    if sha != EXPECTED_SHA: sys.exit('REFUSED: %s has sha256 %s, registered %s' % (PRED_FILE, sha, EXPECTED_SHA))
    return json.loads(raw)

def stem(rng, r2): return os.path.join(HERE, 'data', '%s_r%d' % (rng, r2))

def counts(rng, r2, budget):
    if rng == 'hold': load_pred()          # the holdout is listed and counted only after the predictions are frozen
    import qcounts as Q
    os.makedirs(os.path.join(HERE, 'data'), exist_ok=True); lo, hi = RANGES[rng]
    Q.counts(stem(rng, r2), lo, hi, r2, budget, log=lambda s: print(s, flush=True))

def coarse(x, fi): return np.add.reduceat(x, fi[:-1], axis=-1)[..., :len(fi) - 1]

def jk(num, den):
    n = len(num); N = num.sum(0); Dn = den.sum(0); m = N / Dn; mj = (N[None] - num) / (Dn[None] - den)
    dm = mj - mj.mean(0); return m, (n - 1) / n * dm.T @ dm, mj

def gls(m, Ci, P):
    a = float(P @ Ci @ m / (P @ Ci @ P)); r = m - a * P; return a, float(1 / math.sqrt(P @ Ci @ P)), float(r @ Ci @ r)

def diag_amp(m, mj, P, w):
    n = len(mj); s = float(np.sum(P * m * w) / np.sum(P * P * w)); sj = (mj * (P * w)[None]).sum(1) / np.sum(P * P * w)
    return s, float(math.sqrt((n - 1) / n * np.sum((sj - sj.mean())**2)))

def cond_model(base, types, small_l):
    """finite-X conditional reading: Bhargava local factors at l <= 13 replaced by the observed type frequencies (marginals) of the same fields."""
    import pred4 as P
    from collections import Counter
    w = base.c + 1j * base.t; Dv = base.D.copy(); r3 = base.poles[1][1]; n = len(types)
    for k, l in enumerate(small_l):
        fr = {nm: c / n for nm, c in Counter(types[:, k]).items()}; X = np.exp(-w * math.log(l))
        Dv = Dv * P.E_local(l, X, fr) / P.E_local(l, X); x3 = l**(-1 / 3); r3 *= P.E_local(l, x3, fr) / P.E_local(l, x3)
    return P.Model(base.c, base.t, Dv, [base.poles[0], (1 / 3, r3)])

def read(rng):
    pred = load_pred(); import pred4 as P
    ue = np.array(pred['u_edges']); ye = np.array(pred['y_edges']); fi = np.round(ue * 96).astype(int); h0, h1 = pred['hi_bins']; hi = slice(h0, h1)
    lo = slice(1, h0); dec = pred['decisions']; res = {'range': rng, 'sig': {}}; ok = True
    base = P.quartic_model(T=dec['cond_T'])
    for r2 in (1, 2, 0):
        key = str(r2); z = np.load(stem(rng, r2) + '.npz')
        if not z['done'].all(): sys.exit('counts incomplete for %s r2=%d: %d of %d fields' % (rng, r2, z['done'].sum(), len(z['done'])))
        num, den = coarse(z['num'], fi), coarse(z['den'], fi); n = len(num); S = pred['sig'][key]; Pf = np.array(S['pred'])
        m, C, mj = jk(num[:, hi], den[:, hi]); sig = np.sqrt(np.diag(C)); w = 1 / sig**2
        Ci = np.linalg.inv(C) if n >= dec['min_fields_full_cov'] else np.diag(w)
        a_g, a_ge, c2 = gls(m, Ci, Pf[hi]); s_d, s_de = diag_amp(m, mj, Pf[hi], w)
        comps = {'zero': float(m @ Ci @ m) - c2}
        for nm, Cv in S['comparators'].items(): comps[nm] = gls(m, Ci, np.array(Cv)[hi])[2] - c2
        corr = float(np.corrcoef(m, Pf[hi])[0, 1]); corr_c = {nm: float(np.corrcoef(m, np.array(Cv)[hi])[0, 1]) for nm, Cv in S['comparators'].items()}
        mL, CL, mjL = jk(num[:, lo], den[:, lo]); sL, sLe = diag_amp(mL, mjL, Pf[lo], 1 / np.diag(CL))
        lo_, hi_ = RANGES[rng]; f_count = pred['bhargava_c'][key] * (hi_ - lo_) / n
        cm = cond_model(base, z['types'], tuple(dec['cond_l'])); Pc = P.bin_means(cm, tuple(S['mu']), ye)
        s_c, s_ce = diag_amp(m, mj, Pc[hi], w); a_cg, a_cge, c2c = gls(m, Ci, Pc[hi])
        D = {'D1_comparators': bool(all(v >= dec['delta_min'] for v in comps.values())), 'D2_correlation': bool(corr >= dec['corr_min']),
             'S1_shape_chi2': bool(c2 <= dec['shape_chi2_max']), 'S2_limit_amplitude_1': bool(abs(s_d - 1) <= 3 * s_de),
             'S3_conditional_amplitude_1': bool(abs(s_c - 1) <= 3 * s_ce), 'S4_smally_count_factor': bool(abs(sL / f_count - 1) <= 3 * sLe / f_count)}
        prim = key in dec['primary_sigs']
        if prim: ok &= D['D1_comparators'] and D['D2_correlation']
        r = {'n_fields': n, 'gls_amplitude': a_g, 'gls_amplitude_err': a_ge, 'shape_chi2_gls': c2, 'diag_amplitude': s_d, 'diag_amplitude_err': s_de,
             'delta_chi2': comps, 'corr': corr, 'corr_comparators': corr_c, 's_lo': sL, 's_lo_err': sLe, 'count_factor_f': f_count,
             'cond_diag_amplitude': s_c, 'cond_diag_amplitude_err': s_ce, 'cond_gls_amplitude': a_cg, 'cond_shape_chi2_gls': c2c, 'decisions': D,
             'm': m.tolist(), 'sig': sig.tolist(), 'pred': Pf[hi].tolist(), 'pred_cond': Pc[hi].tolist()}
        if rng == 'hold': r['sigma_ratio_obs_over_expected'] = float(np.median(sig / np.array(S['sigma_hold_expected'])))
        res['sig'][key] = r
        print('%s r2=%d (%s): %d fields' % (rng, r2, {0: 'totally real', 1: 'one complex place, D<0', 2: 'totally complex'}[r2], n))
        print('   amplitude: GLS %.3f +- %.3f (shape chi2 %.1f / 20 bins), diagonal %.3f +- %.3f; correlation %.3f' % (a_g, a_ge, c2, s_d, s_de, corr))
        print('   delta chi2: ' + ', '.join('%s %.0f' % (k_, v) for k_, v in comps.items()) + '; corr of comparators: ' + ', '.join('%s %.2f' % kv for kv in corr_c.items()))
        print('   small y (u < 7/24): s = %.3f +- %.3f, count factor N_main/N = %.3f; conditional (l<=13 observed types): diag %.3f +- %.3f, GLS %.3f (chi2 %.1f)' % (sL, sLe, f_count, s_c, s_ce, a_cg, c2c))
        print('   decisions: ' + ', '.join('%s %s' % (k_, 'PASS' if v else 'FAIL') for k_, v in D.items()) + ('   [D1, D2 primary]' if prim else '   [reported only]'))
    res['verdict'] = 'PRIMARY CONFIRMED' if ok else 'PRIMARY FAILED'
    print('verdict (%s): %s' % (rng + (', dry run on discovery data' if rng == 'disc' else ''), res['verdict']))
    json.dump(res, open(os.path.join(HERE, 'result_%s.json' % rng), 'w'), indent=1)

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('mode', choices=['counts', 'read']); ap.add_argument('--range', choices=['disc', 'hold'], required=True)
    ap.add_argument('--r2', type=int, choices=[0, 1, 2]); ap.add_argument('--budget', type=float, default=540.0)
    A = ap.parse_args()
    if A.mode == 'counts':
        if A.r2 is None: sys.exit('counts needs --r2')
        counts(A.range, A.r2, A.budget)
    else:
        read(A.range)
