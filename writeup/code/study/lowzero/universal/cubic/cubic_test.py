"""cubic_test.py: frozen test of the predicted murmuration density of S3 cubic fields.

  counts --range disc|hold --sign -1|1 --budget SEC   build the field list (PARI nflist, in chunks) and the per-field
                                                       bin sums; resumable (state in cubic_<range>_<sign>.npz)
  read   --range disc|hold [--sign -1|1]               statistics and decisions; writes cubic_result_<range>.json
Ranges: disc = 5000 <= |D| < 10000 (discovery), hold = 20000 <= |D| < 40000 (holdout). Primes p <= 3|D|, p not dividing D.
Per field and bin j (y = p/|D|): num_j = sum ln p * sqrt|D| * a_K(p), den_j = sum ln p; estimate m_j = sum num / sum den.
The script refuses to run if predictions_cubic.json does not have the registered sha256.
"""
import sys, os, json, time, math, hashlib, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PRED_FILE = os.path.join(HERE, 'predictions_cubic.json')
EXPECTED_SHA = '301e02f06efc0c07c390dc271d9c6faace8cfcffca287138e2a4477a8fed07c1'
RANGES = {'disc': (5000, 10000), 'hold': (20000, 40000)}
YMAX = 3.0
EDGES = np.round(np.linspace(0.0, YMAX, 61), 10)   # 60 bins of width 0.05 in y = p/|D|
CHUNK = 2500        # width of the |D| chunks given to nflist

def load_pred():
    raw = open(PRED_FILE, 'rb').read(); sha = hashlib.sha256(raw).hexdigest()
    if sha != EXPECTED_SHA:
        sys.exit('REFUSED: %s has sha256 %s, registered %s' % (PRED_FILE, sha, EXPECTED_SHA))
    return json.loads(raw)

def pari_init():
    import cypari2
    pari = cypari2.Pari(); pari.allocatemem(200_000_000); pari.default('parisizemax', 400_000_000)
    pari('acnt(f, P) = apply(p -> #polrootsmod(f, p) - 1, P)')
    pari('acnt1(f, p) = my(nf = nfinit(f)); #select(q -> q.f == 1, idealprimedec(nf, p)) - 1')
    return pari

def field_list(pari, rng, sgn):
    fn = os.path.join(HERE, 'fields_%s_%+d.txt' % (rng, sgn)); lo, hi = RANGES[rng]
    done = {}
    if os.path.exists(fn):
        for line in open(fn):
            if line.startswith('#chunk'): done[int(line.split()[1])] = True
    with open(fn, 'a') as out:
        for c0 in range(lo, hi, CHUNK):
            if c0 in done: continue
            t0 = time.time(); c1 = min(hi, c0 + CHUNK) - 1
            fl = pari('nflist("S3", [%d, %d], %d)' % (c0, c1, 1 if sgn < 0 else 0))
            for f in fl:
                D = int(pari.nfdisc(f)); out.write('%d %s\n' % (abs(D) * (1 if D > 0 else -1), str(f).replace(' ', '')))
            out.write('#chunk %d %d fields %.1f s\n' % (c0, len(fl), time.time() - t0)); out.flush()
            print('  nflist chunk [%d, %d]: %d fields, %.1f s' % (c0, c1, len(fl), time.time() - t0), flush=True)
    rows = [l.split() for l in open(fn) if not l.startswith('#')]
    rows = sorted(set((int(D), f) for D, f in rows), key=lambda r: (abs(r[0]), r[1]))
    for D, f in rows: assert (D > 0) == (sgn > 0) and lo <= abs(D) < hi
    return rows

def field_sums(pari, f, D, edges):
    A = abs(D); P = pari.primes([2, int(YMAX * A)])
    a = np.array([int(v) for v in pari('acnt')(f, P)], float); p = np.array([int(v) for v in P], np.int64)
    keep = (A % p != 0)
    disc_f = abs(int(pari.poldisc(f))); ind = disc_f // A
    for i in np.nonzero(keep & (ind % p == 0))[0]:      # index divisors: count the degree-1 primes directly
        a[i] = int(pari('acnt1')(f, int(p[i])))
    p = p[keep]; a = a[keep]; y = p / A
    j = np.searchsorted(edges, y, side='left') - 1     # bin j: edges[j] < y <= edges[j+1]
    ok = (j >= 0) & (j < len(edges) - 1)
    lp = np.log(p[ok].astype(float))
    num = np.bincount(j[ok], weights=lp * math.sqrt(A) * a[ok], minlength=len(edges) - 1)
    den = np.bincount(j[ok], weights=lp, minlength=len(edges) - 1)
    return num, den, len(p)

def counts(rng, sgn, budget):
    if rng == 'hold': load_pred()          # the holdout is counted only after the predictions are frozen
    edges = EDGES; pari = pari_init(); t0 = time.time()
    rows = field_list(pari, rng, sgn)
    st = os.path.join(HERE, 'cubic_%s_%+d.npz' % (rng, sgn))
    if os.path.exists(st):
        z = np.load(st); num, den, D, npr, done = z['num'], z['den'], z['D'], z['npr'], z['done']
        assert len(D) == len(rows) and all(int(D[i]) == rows[i][0] for i in range(len(rows)))
    else:
        nb = len(edges) - 1; num = np.zeros((len(rows), nb)); den = np.zeros((len(rows), nb))
        D = np.array([r[0] for r in rows], np.int64); npr = np.zeros(len(rows), np.int64); done = np.zeros(len(rows), bool)
    todo = np.nonzero(~done)[0]; nrun = 0; tr = time.time()
    for i in todo:
        if time.time() - t0 > budget: break
        D_, f = rows[i]; num[i], den[i], npr[i] = field_sums(pari, pari(f), D_, edges); done[i] = True; nrun += 1
        if nrun % 100 == 0: np.savez(st, num=num, den=den, D=D, npr=npr, done=done)
    np.savez(st, num=num, den=den, D=D, npr=npr, done=done)
    dt = time.time() - tr
    print('%s sign %+d: %d fields in list, %d done (%d this run, %.1f s, %.2f us per root count)%s' %
          (rng, sgn, len(rows), int(done.sum()), nrun, dt, 1e6 * dt / max(1, int(npr[todo[:nrun]].sum()) if nrun else 1),
           '  DONE' if done.all() else '  (run again to continue)'))

# ---------------------------------------------------------------------------------------------------------------
def jack(num, den, P=None, w=None):  # w: fixed weights 1/sigma^2 (full-sample jackknife)
    """delete-one jackknife over fields: estimate, covariance; amplitude s = <P, m>_w / <P, P>_w and its error."""
    n = num.shape[0]; N = num.sum(0); Dn = den.sum(0); m = N / Dn
    mj = (N[None, :] - num) / (Dn[None, :] - den)
    dm = mj - mj.mean(0); cov = (n - 1) / n * dm.T @ dm
    out = {'m': m, 'cov': cov, 'sig': np.sqrt(np.diag(cov))}
    if P is not None:
        s = np.sum(P * m * w) / np.sum(P * P * w); sj = (mj * (P * w)[None, :]).sum(1) / np.sum(P * P * w)
        out['s'] = s; out['s_err'] = math.sqrt((n - 1) / n * np.sum((sj - sj.mean())**2))
    return out

def coarse(x, g):
    return x.reshape(x.shape[:-1] + (x.shape[-1] // g, g)).sum(-1)

def chi2(r, cov):
    return float(r @ np.linalg.solve(cov, r))

def read(rng, signs, dry):
    pred = load_pred(); edges = np.array(pred['edges'], float); G = pred['coarse_group']
    assert np.allclose(edges, EDGES)
    hi = slice(*pred['range_hi_bins']); lo = slice(*pred['range_lo_bins']); res = {'range': rng, 'signs': {}}; primary_ok = True
    for sgn in signs:
        key = '%+d' % sgn; st = os.path.join(HERE, 'cubic_%s_%+d.npz' % (rng, sgn))
        z = np.load(st); num, den, done = z['num'], z['den'], z['done']
        if not done.all(): sys.exit('counts incomplete for %s %+d: %d of %d fields' % (rng, sgn, done.sum(), len(done)))
        P = np.array(pred['pred'][key]); d = pred['decisions'][key]; f = pred['finiteX_factor'][key]['disc' if rng == 'disc' else 'hold']
        Jf = jack(num, den); w = 1 / Jf['sig']**2
        Jh = jack(num[:, hi], den[:, hi], P[hi], w[hi]); Jl = jack(num[:, lo], den[:, lo], P[lo], w[lo])
        Jc = jack(coarse(num[:, hi], G), coarse(den[:, hi], G)); Pc = coarse(P[hi], G) / G
        c2_pred = chi2(Jc['m'] - Pc, Jc['cov']); c2_shape = chi2(Jc['m'] - Jh['s'] * Pc, Jc['cov'])
        comps = {n: chi2(Jc['m'] - coarse(np.array(C)[hi], G) / G, Jc['cov']) - c2_pred for n, C in pred['comparators'][key].items()}
        dec = {'A_amplitude': bool(abs(Jh['s'] - 1) <= d['amp_z'] * Jh['s_err'] and Jh['s'] / Jh['s_err'] >= d['amp_signal']),
               'B_shape': bool(c2_shape <= d['shape_chi2_max']),
               'C_comparators': bool(all(v >= d['delta_min'] for v in comps.values())),
               'D_smally_finiteX': bool(abs(Jl['s'] / f - 1) <= d['lo_z'] * Jl['s_err'] / f and (Jl['s'] - 1) >= d['lo_excess_z'] * Jl['s_err'])}
        r = {'n_fields': int(len(done)), 's_hi': Jh['s'], 's_hi_err': Jh['s_err'], 'chi2_coarse_pred': c2_pred, 'chi2_coarse_shape': c2_shape,
             's_lo': Jl['s'], 's_lo_err': Jl['s_err'], 'finiteX_factor': f, 'delta_chi2_comparators': comps, 'decisions': dec,
             'm': Jf['m'].tolist(), 'sig': Jf['sig'].tolist()}
        if rng == 'hold': r['sigma_ratio_obs_over_registered'] = float(np.median(Jf['sig'] / np.array(pred['sigma_hold'][key])))
        res['signs'][key] = r
        primary_ok &= (dec['A_amplitude'] and dec['B_shape'] and dec['C_comparators'] and dec['D_smally_finiteX']) if d['primary'] else dec['D_smally_finiteX']
        print('%s D %s 0: %d fields' % (rng, '<' if sgn < 0 else '>', r['n_fields']))
        print('   y in (0.4,3]: amplitude s = %.4f +- %.4f; coarse chi2 (13 bins) pred %.1f, shape (s fitted) %.1f [max %.2f]' % (Jh['s'], Jh['s_err'], c2_pred, c2_shape, d['shape_chi2_max']))
        print('   delta chi2 of the comparators: ' + ', '.join('%s %.0f' % (n, v) for n, v in comps.items()))
        print('   y in (0,0.2]: s = %.4f +- %.4f; finite-X factor 1/(1+eps) = %.4f' % (Jl['s'], Jl['s_err'], f))
        print('   decisions: ' + ', '.join('%s %s' % (n, 'PASS' if v else 'FAIL') for n, v in dec.items()) + ('   [A-C primary]' if d['primary'] else '   [A-C secondary, D primary]'))
    res['verdict'] = 'PRIMARY CONFIRMED' if primary_ok else 'PRIMARY FAILED'
    print('verdict (%s): %s' % (rng + (', dry run on discovery data' if dry else ''), res['verdict']))
    json.dump(res, open(os.path.join(HERE, 'cubic_result_%s.json' % rng), 'w'), indent=1)

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('mode', choices=['counts', 'read']); ap.add_argument('--range', choices=['disc', 'hold'], required=True)
    ap.add_argument('--sign', type=int, choices=[-1, 1]); ap.add_argument('--budget', type=float, default=500.0)
    A = ap.parse_args()
    if A.mode == 'counts':
        if A.sign is None: sys.exit('counts needs --sign')
        counts(A.range, A.sign, A.budget)
    else:
        read(A.range, [A.sign] if A.sign else [-1, 1], A.range == 'disc')
