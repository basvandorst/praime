"""freeze_spectro.py: the frozen predictions for the HOLDOUT range (3e7, 1e9] (nothing there is computed here).

Everything below comes from formulas (spectro.md sections 1-2) and from the DISCOVERY outputs (only for the error model):
  constellation statistics (twin_lib): twinlo, twinhi, sg for EC:11a1, EC:15a1, A4, C23; shift2, shift6, single for A4, C23
  lines (lines_lib): the 60 strongest predicted lines of A4 and of C23 on the holdout window + registered extras.
Window: Hann in v = sqrt(n) over [sqrt(3e7), sqrt(1e9)].
Predicted errors: sigma_hold = sigma_disc * sqrt(I_hold / I_disc) * |MT_disc| / |MT_hold| (random-model integral I of the
weights; sigma_disc = max(control rms, random model) from out_twins_discovery.json, resp. the random-model sigma of
scan_lines.json).  Output: predictions_spectro.json and its sha256 in freeze_sha256.txt."""
import sys, os, json, math, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
import numpy as np
import spectro_lib as S
import twin_lib as T
import lines_lib as LL

H = S.HERE
V0, V1 = math.sqrt(3e7), math.sqrt(1e9)
D0, D1 = math.sqrt(6e5), math.sqrt(3e7)            # discovery window of twins_discovery.py
disc = json.load(open(os.path.join(H, 'out_twins_discovery.json')))

def integ(f, a, b, n=200001):
    v = np.linspace(a, b, n); y = f(v); h = v[1] - v[0]
    return float(h / 3 * (y[0] + y[-1] + 4 * y[1:-1:2].sum() + 2 * y[2:-1:2].sum()))

def noise_int(kind, a, b):
    """random-model integral of the weights squared, up to a constant common to discovery and holdout."""
    H2 = lambda v: S.hann(v, a, b) ** 2
    if kind in ('twinlo', 'twinhi'): return integ(lambda v: np.log(v * v) ** 2 * H2(v), a, b)
    if kind == 'sg': return integ(lambda v: np.log(v * v) * np.log(2 * v * v) * H2(v), a, b)
    return integ(lambda v: np.log(v * v) * H2(v), a, b)        # single, shift2, shift6

pred = {'window_v': [V0, V1], 'range': [30_000_000, 1_000_000_000], 'two_C2': T.TWO_C2, 'stats': {}, 'lines': {}}
STATS = {'EC:11a1': ['twinlo', 'twinhi', 'sg'], 'EC:15a1': ['twinlo', 'twinhi', 'sg'],
         'A4': ['twinlo', 'twinhi', 'sg', 'single', 'shift2', 'shift6'], 'C23': ['twinlo', 'twinhi', 'sg', 'single', 'shift2', 'shift6']}
for lab, sts in STATS.items():
    c = S.label_info(lab)['c']
    MTh = S.main_term_unit(lab, V0, V1); MTd = S.main_term_unit(lab, D0, D1)
    for st in sts:
        C = disc['consts'][lab][st]
        d = disc['res'][lab][st]
        sig_d = max(d['err_ctrl'], d['err_random'])
        sig_h = sig_d * math.sqrt(noise_int(st, V0, V1) / noise_int(st, D0, D1)) * abs(MTd) / abs(MTh)
        if st in ('twinlo', 'twinhi', 'sg'):
            alts = {'derived': 1.0, 'factorised_independence': T.TWO_C2 * c / C, 'no_resonance': 0.0,
                    'integer_like_label': T.TWO_C2 / C, 'no_twin_enhancement': c / C}
        elif st in ('shift2', 'shift6'):
            alts = {'derived': 1.0, 'integer_like_label': 1.0 / C, 'prime_like': c / C, 'no_resonance': 0.0}
        else:
            alts = {'derived': 1.0, 'no_resonance': 0.0, 'constant_1': 1.0 / C}
        Sp = C * MTh
        pred['stats']['%s|%s' % (lab, st)] = dict(label=lab, stat=st, const=C, MT_unit=[MTh.real, MTh.imag], S_pred=[Sp.real, Sp.imag],
                                                 reading_expected=[1.0, 0.0], sigma_pred=sig_h, alternatives=alts,
                                                 discovery_reading=d['reading'], discovery_sigma=sig_d)
# ------------------------------------------------------------------ lines (A4, C23), weight ln p Hann(sqrt p)
sl = json.load(open(os.path.join(H, 'scan_lines.json')))
U0d, U1d = math.sqrt(1e5), math.sqrt(3e7)
def mt_scan(label, r, a, b, npts=20001):
    K1, kappa = S.K1_kappa(label); N = S.label_info(label)['N']
    uu = np.linspace(a, b, npts); f = S.hann(uu, a, b) * uu ** 0.5 * (1 - 1j * kappa * math.sqrt(N) / (S.TWO_PI * uu * math.sqrt(r)))
    h = uu[1] - uu[0]; s = f[0] + f[-1] + 4 * f[1:-1:2].sum() + 2 * f[2:-1:2].sum()
    return 1.5 * K1 * r ** -0.25 * s * h / 3
lnint = lambda a, b: integ(lambda v: np.log(v * v) * S.hann(v, a, b) ** 2 * v, a, b)
for key, lab in (('a4', 'A4'), ('c23', 'C23')):
    N = S.label_info(lab)['N']; LL._LOC[lab] = S.local_lams(lab, 20000); BAD = json.load(open(LL.bad_factor_file(lab)))
    sig_h = sl[key]['sigma'] * math.sqrt(lnint(V0, V1) / lnint(U0d, U1d))
    cand = []
    for b in range(1, 37):
        if not all(e <= 2 for e in LL.factor(b).values()): continue
        for a in range(1, int(25.0 * N / 4 * b) + 1):
            if math.gcd(a, b) != 1: continue
            r = Fr(a, b); al = 2 * math.sqrt(float(r) / N)
            if al < 0.05 or al > 5: continue
            c = LL.c_line(lab, a, b, BAD)
            if c != c or abs(c) < 1e-12: continue
            cand.append((abs(c * mt_scan(lab, float(r), V0, V1, 2001)) / sig_h, r, c))
    cand.sort(key=lambda x: -x[0])
    chosen = [(r, c) for (z, r, c) in cand[:60]]
    extras = [Fr(205, 26), Fr(5, 9)] if lab == 'A4' else [Fr(926, 9), Fr(277, 2)]      # largest discovery deviations
    zeros = [Fr(1, 8), Fr(3, 8), Fr(1, 27)]                                            # predicted zero (cube in b or 2^-3)
    rows = []
    for r in [x for x, _ in chosen] + [e for e in extras if e not in dict(chosen)] + zeros:
        a, b = r.numerator, r.denominator
        c = LL.c_line(lab, a, b, BAD) if all(e <= 2 for e in LL.factor(b).values()) else 0.0
        if c != c: c = 0.0
        Sp = c * mt_scan(lab, float(r), V0, V1)
        rows.append(dict(r=str(r), alpha=2 * math.sqrt(float(r) / N), c=c, S_pred=[Sp.real, Sp.imag],
                         role=('top60' if r in dict(chosen) else ('extra' if r in extras else 'zero'))))
    den = sum(abs(complex(*x['S_pred'])) ** 2 for x in rows if x['role'] == 'top60')
    pred['lines'][lab] = dict(sigma_pred=sig_h, slope_expected=1.0, slope_sigma_pred=sig_h / math.sqrt(den), rows=rows)
    print(lab, 'lines: top z %.1f, 60th z %.1f, sigma %.1f, slope sigma %.4f' % (cand[0][0], cand[59][0], sig_h, sig_h / math.sqrt(den)))
# ------------------------------------------------------------------ registered decisions
pred['decisions'] = {
 'D1_primary': 'pooled twinlo reading (inverse-variance with sigma_pred, real parts, 4 labels) within 1 +- 3 sigma_pooled',
 'D2_models': 'over all 12 twin readings (twinlo, twinhi, sg x 4 labels; real and imaginary parts, sigma_pred): chi2(derived) < chi2 of '
              'no_resonance, integer_like_label, no_twin_enhancement, each by >= 25 where the predicted separation (sum ((alt-1)/sigma)^2) is >= 50',
 'D3_shift': 'A4 shift2 and A4 shift6 each within 3 sigma_pred of 1 (real part), and chi2(integer_like_label) - chi2(derived) >= 25 on the 4 shift readings',
 'D4_lines': 'for A4 and for C23 the pooled slope of the 60 top lines (projection on S_pred) within 1 +- max(3 sigma, 0.02); '
             'the registered extra lines each within 4 sigma_pred of S_pred; the zero lines each |S| < 4 sigma_pred',
 'D5_controls': 'the control-frequency readings (8 per statistic) are REPORTED with their rms in units of sigma_pred; no decision',
 'falsifier': 'pooled twinlo reading outside [0.80, 1.20], or A4 shift2 outside [0.90, 1.10], or a line slope outside [0.95, 1.05]',
}
blob = json.dumps(pred, indent=1, sort_keys=True)
open(os.path.join(H, 'predictions_spectro.json'), 'w').write(blob)
h = hashlib.sha256(blob.encode()).hexdigest()
open(os.path.join(H, 'freeze_sha256.txt'), 'w').write(h + '  predictions_spectro.json\n')
print('sha256', h)
for k, v in pred['stats'].items():
    print('%-16s const %.4f  sigma_pred %.3f  (disc reading %+.3f, sigma %.3f)  alts %s' % (k, v['const'], v['sigma_pred'], v['discovery_reading'][0], v['discovery_sigma'],
          {a: round(x, 3) for a, x in v['alternatives'].items()}))
