"""ps_freeze.py: the frozen predictions of the prime-side test (predictions_primeside.json) and its sha256.

For every series and both runs (dry = seen range, fresh = registered range) it writes:
  K (complex), kappa, c, eps; the main term MT of the sharp increment over (x0, x1] and of the smooth sums (main Hann window
  and 8 sub-windows), including the first-order kernel factor (1 - i kappa/t); the prime-power term PP (lambda-Lambda rule:
  the primes-only sum is predicted as MT - sum_{p^k in range, k >= 2} lambda(p^k) ln p e(...) * weight); the alternatives
  'constant 1' (MT/c) and 'no effect' (0); and for Q(i) the ray-model predictions of the Pythagorean zone sums.
usage: ps_freeze.py   (needs constants.json, small_lams.npz from ps_constants.py)"""
import sys, os, json, math, hashlib, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ps_lib as L
import ps_sums as S

HERE = os.path.dirname(os.path.abspath(__file__))
const = json.load(open(os.path.join(HERE, 'constants.json')))
sm = np.load(os.path.join(HERE, 'small_lams.npz'))

def cplx(z): return [float(z.real), float(z.imag)]

def lam_powers(key, p, lam, kmax):
    """lambda(p^k), k = 0..kmax (rows), for arrays p, lam."""
    d = L.SERIES[key][1]; out = np.zeros((kmax + 1, len(p))); out[0] = 1.0; out[1] = lam
    if d == 2:
        chi = L.chi_neben(key, p).astype(float)
        for k in range(2, kmax + 1): out[k] = lam * out[k - 1] - chi * out[k - 2]
    else:   # Sym^2 of 11a: generating function 1/((1-X)(1-(s)X+X^2)), s+1 = lam_S; bad p = 11: 11^-k
        s1 = lam
        for k in range(2, kmax + 1):
            out[k] = s1 * out[k - 1] - s1 * out[k - 2] + (out[k - 3] if k >= 3 else 0.0)
        bad = p == 11
        for k in range(kmax + 1): out[k][bad] = 11.0 ** (-k)
    return out

def pp_terms(key, x0, x1, W):
    """prime-power term: sum over p^k in (x0, x1], k >= 2, of lambda(p^k) ln p e(-d (p^k/Q)^(1/d)) [sharp] and with the smooth weights."""
    grp, d, Q, mu, k = L.SERIES[key]
    p = sm[key + '_p']; lam = sm[key + '_lam']
    if key[0] in 'AB':      # small_lams holds all primes (inert with lam = 0, specials included)
        pass
    kmax = int(math.log(x1) / math.log(2)) + 1
    lp = lam_powers(key, p.astype(np.int64), lam, kmax)
    sharp = 0j; smooth = np.zeros(len(W), complex); n = 0
    for kk in range(2, kmax + 1):
        n_ = p.astype(np.float64) ** kk; m = (n_ > x0) & (n_ <= x1)
        if not m.any(): continue
        nn = n_[m]; wt = lp[kk][m] * np.log(p[m].astype(np.float64)) * np.exp(-2j * math.pi * d * (nn / Q) ** (1 / d))
        sharp += wt.sum(); n += int(m.sum())
        v = nn ** (1 / d); g = nn ** (-0.5 + 0.5 / d)
        for i, (a, b) in enumerate(W): smooth[i] += (wt * g * L.hann(v, a, b)).sum()
    assert p.max() ** 2 > x1, 'small_lams does not reach sqrt(x1)'
    return sharp, smooth, n

def predictions():
    out = {'what': 'prime-side test of the one-point law (fact 2055): frozen predictions', 'created': time.strftime('%Y-%m-%d %H:%M:%S'),
           'config': {'runs': S.RUNS, 'chunk': S.CHUNK, 'ctrl': L.CTRL, 'nsub': 8, 'beta': L.BETA, 'zone_dirs': S.ZONE_DIRS, 'zone_ms': list(S.ZONE_MS)},
           'series': {}, 'zones': {}}
    for key, (grp, d, Q, mu, k) in L.SERIES.items():
        cs = const[key]; c = cs['c_exact']; eps = complex(*cs['eps'])
        K, kappa, phi0 = L.law_constants(d, Q, mu, eps, c)
        ent = {'d': d, 'Q': Q, 'mu': mu, 'k': k, 'eps': cplx(eps), 'c': c, 'c_method': cs['c_method'], 'K': cplx(K), 'K_abs': abs(K),
               'K_phase_deg': math.degrees(math.atan2(K.imag, K.real)), 'kappa': kappa, 'phi0': phi0, 'exponent': (d + 1) / (2 * d), 'runs': {}}
        for run, cfgs in S.RUNS.items():
            cfg = cfgs['ABC' if grp[0] in 'ABC' else ('DD' if key == 'DD' else 'D')]
            x0, x1 = cfg['read']; W = L.windows(x0 ** (1 / d), x1 ** (1 / d))
            MTs = L.main_term(K, kappa, d, Q, x0, x1)
            MTw = [L.main_term(K, kappa, d, Q, x0, x1, weight=(lambda v, a=a, b=b: v ** ((1 - d) / 2) * L.hann(v, a, b))) for (a, b) in W]
            MT0 = L.main_term(K, 0.0, d, Q, x0, x1)
            PPs, PPw, npp = pp_terms(key, x0, x1, W)
            # sharp main term at the grid points from x0 (for the curve of the reading)
            grid = [g for g in cfg['grid'] if g >= x0]
            MTgrid = [L.main_term(K, kappa, d, Q, x0, g) if g > x0 else 0j for g in grid]
            PPgrid = [pp_terms(key, x0, g, W)[0] if g > x0 else 0j for g in grid]
            ent['runs'][run] = {'read': [x0, x1], 'windows': W, 'MT_sharp': cplx(MTs), 'MT_sharp_no_kappa': cplx(MT0), 'PP_sharp': cplx(PPs),
                                'MT_smooth': [cplx(z) for z in MTw], 'PP_smooth': [cplx(z) for z in PPw], 'n_prime_powers': npp,
                                'grid': grid, 'MT_grid': [cplx(z) for z in MTgrid], 'PP_grid': [cplx(z) for z in PPgrid]}
        out['series'][key] = ent
        print('%-5s c = %.6f eps = %+.0f  K = %.6f %+.6fi (|K| %.4f, %.1f deg)  kappa = %.4f   fresh: MT_sharp = %.4g, PP/MT = %.2e' % (
            key, c, eps.real, K.real, K.imag, abs(K), ent['K_phase_deg'], kappa, abs(complex(*ent['runs']['fresh']['MT_sharp'])),
            abs(complex(*ent['runs']['fresh']['PP_sharp'])) / abs(complex(*ent['runs']['fresh']['MT_sharp']))), flush=True)
    # Pythagorean zones (Q(i), primes, m = 0 and 4)
    dirs = L.pyth_directions(5000)
    for run, cfgs in S.RUNS.items():
        x0, x1 = cfgs['ABC']['read']; zr = {}
        for name, ph in S.ZONE_DIRS.items():
            zp = L.zone_prediction(ph, S.ZONE_MS, math.sqrt(x0), math.sqrt(x1), dirs)
            zr[name] = {str(m): cplx(zp[m]) for m in S.ZONE_MS}
        # single-ray scale (one full Pythagorean ray of weight 4/pi, direction-independent): for errors in units of the axis
        zr['_unit_ray'] = cplx(4 / math.pi * L.zone_F(0.0) * (2 / 3) * (x1 ** 0.75 - x0 ** 0.75))
        out['zones'][run] = zr
        print(run, 'zones:', ' '.join('%s %.0f%+.0fi' % (n, v['0'][0], v['0'][1]) for n, v in zr.items() if n[0] != '_'), flush=True)
    return out

if __name__ == '__main__':
    pr = predictions()
    pr['decisions'] = json.load(open(os.path.join(HERE, 'decisions.json')))
    fn = os.path.join(HERE, 'predictions_primeside.json')
    s = json.dumps(pr, indent=1, sort_keys=True)
    open(fn, 'w').write(s)
    h = hashlib.sha256(s.encode()).hexdigest()
    open(os.path.join(HERE, 'freeze_sha256.txt'), 'w').write(h + '  predictions_primeside.json\n')
    print('sha256', h)
