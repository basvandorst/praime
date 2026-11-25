"""scan_lines.py: every predicted line r = a/b (b cube-free <= 36) of the four automorphic labels at beta = 1/2 on the
DISCOVERY range, measured exactly (direct sums over primes in [1e5, 3e7), weight ln p * Hann(sqrt p)), against the
prediction c(r) MT(r) of lines_lib (no fitted parameter).  Lines kept: predicted |S| > 2 sigma, plus every scan peak with
P > 20.  Also: residual peaks (P > 20 at beta = 1/2 whose r has no predicted line above 2 sigma).
Output scan_lines.json, printed summary."""
import sys, os, json, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
import numpy as np
import spectro_lib as S
import lines_lib as LL

Z = np.load(os.path.join(S.HERE, 'labels_3e7.npz'))
p = Z['p']; m = p >= 100000; P = p[m].astype(np.float64); u = np.sqrt(P)
u0, u1 = math.sqrt(1e5), math.sqrt(3e7)
W = np.log(P) * S.hann(u, u0, u1)
LABS = {'l11': 'EC:11a1', 'l37': 'EC:37a1', 'c23': 'C23', 'a4': 'A4'}
peaks = json.load(open(os.path.join(S.HERE, 'scan_peaks.json')))

def mt_scan(label, r, npts=40001):
    K1, kappa = S.K1_kappa(label); N = S.label_info(label)['N']
    uu = np.linspace(u0, u1, npts); f = S.hann(uu, u0, u1) * uu ** 0.5 * (1 - 1j * kappa * math.sqrt(N) / (S.TWO_PI * uu * math.sqrt(r)))
    h = uu[1] - uu[0]; s = f[0] + f[-1] + 4 * f[1:-1:2].sum() + 2 * f[2:-1:2].sum()
    return 1.5 * K1 * r ** -0.25 * s * h / 3

def cubefree(b):
    return all(e <= 2 for e in LL.factor(b).values())

out = {}; t0 = time.time()
for key, label in LABS.items():
    N = S.label_info(label)['N']; w = Z[key][m] * W; sig = math.sqrt(float(np.dot(w, w)) / 2)
    LL._LOC[label] = S.local_lams(label, 20000); BAD = json.load(open(LL.bad_factor_file(label)))
    rmax = 25.0 * N / 4
    cands = {}
    for b in range(1, 37):
        if not cubefree(b): continue
        for a in range(1, int(rmax * b) + 1):
            if math.gcd(a, b) != 1: continue
            r = Fr(a, b)
            if 2 * math.sqrt(float(r) / N) < 0.05: continue
            c = LL.c_line(label, a, b, BAD)
            if c != c: continue
            Sp = c * mt_scan(label, float(r), 2001) if abs(c) > 1e-12 else 0j
            if abs(Sp) > 2 * sig: cands[r] = c
    # every scan peak P > 20 at beta = 1/2
    extra = []
    for pk in peaks:
        if pk['label'] != key or pk['beta'] != '1/2' or pk['P'] <= 20: continue
        r = Fr(pk['alpha'] ** 2 * N / 4).limit_denominator(60)
        if r not in cands: extra.append((r, pk))
    rows = []
    allr = sorted(set(cands) | set(r for r, _ in extra))
    for r in allr:
        a, b = r.numerator, r.denominator
        c = LL.c_line(label, a, b, BAD) if cubefree(b) else 0.0
        if c != c: c = 0.0
        Sp = c * mt_scan(label, float(r))
        al = 2 * math.sqrt(float(r) / N)
        Sm = complex(np.dot(w, np.exp(-1j * S.TWO_PI * al * u)))
        rows.append(dict(r=str(r), alpha=al, c=c, S_pred=[Sp.real, Sp.imag], S_meas=[Sm.real, Sm.imag], sigma=sig,
                         dev=abs(Sm - Sp) / sig, P_meas=abs(Sm) ** 2 / (2 * sig * sig), P_pred=abs(Sp) ** 2 / (2 * sig * sig)))
    # pooled slope (projection of measured on predicted, inverse-variance) and chi2
    num = sum((complex(*x['S_meas']) * complex(*x['S_pred']).conjugate()).real for x in rows)
    den = sum(abs(complex(*x['S_pred'])) ** 2 for x in rows)
    slope = num / den; slope_err = sig / math.sqrt(den)
    chi2 = sum(x['dev'] ** 2 for x in rows); ndof = 2 * len(rows)
    chi2_null = sum(abs(complex(*x['S_meas'])) ** 2 / sig ** 2 for x in rows)
    out[key] = dict(label=label, sigma=sig, n_lines=len(rows), slope=slope, slope_err=slope_err, chi2=chi2, ndof=ndof, chi2_no_lines=chi2_null,
                    unexplained_peaks=[dict(r=str(r), P=pk['P'], alpha=pk['alpha']) for r, pk in extra if not cubefree(r.denominator) or LL.c_line(label, r.numerator, r.denominator, BAD) == 0],
                    rows=rows)
    big = sorted(rows, key=lambda x: -x['P_pred'])[:12]
    print('== %s (%s): %d lines (pred > 2 sigma or scan peak); slope meas/pred = %.4f +- %.4f; chi2 = %.1f / %d (no lines: %.1f); %.0fs'
          % (key, label, len(rows), slope, slope_err, chi2, ndof, chi2_null, time.time() - t0))
    for x in big:
        print('   r=%-7s alpha=%.5f c(r)=%+.4f  S_meas/sig = %+6.2f %+6.2fi   S_pred/sig = %+6.2f %+6.2fi   |dev| %.2f'
              % (x['r'], x['alpha'], x['c'], x['S_meas'][0] / sig, x['S_meas'][1] / sig, x['S_pred'][0] / sig, x['S_pred'][1] / sig, x['dev']))
    worst = sorted(rows, key=lambda x: -x['dev'])[:5]
    print('   largest deviations:', ['%s: %.2f' % (x['r'], x['dev']) for x in worst])
    print('   peaks P > 20 without a predicted line:', out[key]['unexplained_peaks'])
json.dump(out, open(os.path.join(S.HERE, 'scan_lines.json'), 'w'), indent=0)
print('done %.0fs' % (time.time() - t0))
