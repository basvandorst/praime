"""Numbers for the registration draft: measurement errors at the new height (24 windows) scaled from the nine measured
windows (per-window error = 3 x pooled 9-window error, the larger of block and window-scatter errors; the per-block
within value does not depend on the number of primes per block to first order, so no density factor), total errors,
the height-contrast statistic T and its power."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import json, math, numpy as np
from law import MEAS, NEW, OLD, lam_of, pooled, lead_k1, lead_k2_wick
pred = json.load(open('data/predictions.json'))
nW = len(NEW); out = {}
print('cell           law     law_err  sig_meas24  sig_tot  | 1: z   | sieve chain: value z | old measured: value z | leading order')
for (q, k, H), (m, eb, ew) in MEAS.items():
    p = pred['%d_%d_%d' % (q, k, H)]; s9 = max(eb, ew); s24 = s9 * 3 / math.sqrt(nW)
    st = math.sqrt(s24 ** 2 + p['law_err'] ** 2)
    out['%d_%d_%d' % (q, k, H)] = dict(law=p['law'], law_err=p['law_err'], sig24=s24, sig_tot=st, MQ=p['MQ'], old=m, lead=p['lead'])
    print('(%d,%d,%.0e)  %.4f  %.4f   %.4f     %.4f  | %+5.1f | %.4f %+5.1f | %.4f %+5.1f | %.4f'
          % (q, k, H, p['law'], p['law_err'], s24, st, (1 - p['law']) / st, p['MQ'], (p['MQ'] - p['law']) / st, m, (m - p['law']) / st, p['lead']))
# height contrast T over k = 1, H = 1e4, q = 3, 5, 7
num = 0; den = 0; parts = []
for q in (3, 5, 7):
    key = '%d_1_10000' % q; p = out[key]
    m, eb, ew = MEAS[(q, 1, 10 ** 4)]; s9 = max(eb, ew)
    dlaw = p['law'] - pooled(OLD, 'old', q, 1, 10 ** 4)
    dtail = (p['law'] - pooled(NEW, 'new', q, 1, 10 ** 4, 'notail')) - (pooled(OLD, 'old', q, 1, 10 ** 4) - pooled(OLD, 'old', q, 1, 10 ** 4, 'notail'))
    sd = math.sqrt(p['sig24'] ** 2 + s9 ** 2 + (0.4 * dtail) ** 2)
    dlead = float(np.mean([lead_k1(q, lam_of(x)) for x in NEW]) - np.mean([lead_k1(q, lam_of(x)) for x in OLD]))
    parts.append((q, dlaw, sd, dlead, m)); 
sigT = 1 / math.sqrt(sum((d / s) ** 2 for q, d, s, dl, m in parts))
w = [d / s ** 2 for q, d, s, dl, m in parts]; W = sum(wi * d for wi, (q, d, s, dl, m) in zip(w, parts))
Tlead = sum(wi * dl for wi, (q, d, s, dl, m) in zip(w, parts)) / W
print('\nheight contrast T = sum_q w_q (W_new - W_old,measured) / sum_q w_q dlaw_q,  w_q = dlaw_q/sig_q^2')
for (q, d, s, dl, m), wi in zip(parts, w):
    print('  q=%d  dlaw %+.4f  sig_delta %.4f  (dlaw/sig %.1f)  leading-order delta %+.4f  weight %.1f  old measured %.4f' % (q, d, s, d / s, dl, wi, m))
print('  law: T = 1 +- %.3f ; no height dependence: T = 0 (%.1f sigma) ; leading order alone: T = %.2f' % (sigT, 1 / sigT, Tlead))
# slope reading, k = 1, q = 5 and 7, S = W(1e5) - W(1e4) at the new height, mean over the two q
S = []; sm = []; sl = []
for q in (5, 7):
    a = out['%d_1_10000' % q]; b = out['%d_1_100000' % q]
    S.append(b['law'] - a['law']); sm.append(math.hypot(a['sig24'], b['sig24']))
    tail_d = (b['law'] - pooled(NEW, 'new', q, 1, 10 ** 5, 'notail')) - (a['law'] - pooled(NEW, 'new', q, 1, 10 ** 4, 'notail'))
    sl.append(0.4 * abs(tail_d))
Sbar = np.mean(S); smeas = math.sqrt(sum(x ** 2 for x in sm)) / 2; slaw = np.mean(sl)
print('\nslope S = mean_{q=5,7} [W(1e5) - W(1e4)]: law %+.4f, sig_meas %.4f, law err (40%% of tail) %.4f; short-range-only law / sieve chain: %+.4f'
      % (Sbar, smeas, slaw, np.mean([pooled(NEW, 'new', q, 1, 10 ** 5, 'notail') - pooled(NEW, 'new', q, 1, 10 ** 4, 'notail') for q in (5, 7)])))
json.dump(dict(cells=out, T_sigma=sigT, T_lead=Tlead, T_parts=parts, S=Sbar, S_sig_meas=smeas, S_sig_law=slaw), open('data/reg_numbers.json', 'w'), indent=1)
# --- the contrast reading conditional on the (already known) old measured values ---
print('\nContrast reading conditional on the old measured values (they are fixed numbers now):')
rows = []
for q in (3, 5, 7):
    p = out['%d_1_10000' % q]; m = MEAS[(q, 1, 10 ** 4)][0]
    dlaw = p['law'] - pooled(OLD, 'old', q, 1, 10 ** 4)
    dtail = (p['law'] - pooled(NEW, 'new', q, 1, 10 ** 4, 'notail')) - (pooled(OLD, 'old', q, 1, 10 ** 4) - pooled(OLD, 'old', q, 1, 10 ** 4, 'notail'))
    sig = math.sqrt(p['sig24'] ** 2 + (0.4 * dtail) ** 2 + 0.0005 ** 2)
    rows.append((q, dlaw, sig, p['law'] - m, p['lead'] - m, p['MQ'] - m))
w = [d / s ** 2 for q, d, s, a, b, c in rows]; den = sum(wi * d for wi, (q, d, s, a, b, c) in zip(w, rows))
ET = sum(wi * a for wi, (q, d, s, a, b, c) in zip(w, rows)) / den
ETlead = sum(wi * b for wi, (q, d, s, a, b, c) in zip(w, rows)) / den
ETmq = sum(wi * c for wi, (q, d, s, a, b, c) in zip(w, rows)) / den
sT = math.sqrt(sum((wi * s) ** 2 for wi, (q, d, s, a, b, c) in zip(w, rows))) / abs(den)
for (q, d, s, a, b, c), wi in zip(rows, w):
    print('  q=%d  weight %.1f  sigma %.4f  law_new - old_meas %+.4f  lead_new - old_meas %+.4f  sieve_new - old_meas %+.4f' % (q, wi, s, a, b, c))
print('  T = sum w (W_new - W_old_meas) / %.4f;  E[T | law] = %.3f +- %.3f;  no height dependence: 0 (%.1f sigma);  leading order: %.2f (%.1f sigma);  sieve chain: %.2f (%.1f sigma)'
      % (den, ET, sT, ET / sT, ETlead, (ET - ETlead) / sT, ETmq, (ET - ETmq) / sT))
json.dump(dict(weights=w, den=den, ET=ET, sT=sT, ETlead=ETlead, ETmq=ETmq, rows=rows), open('data/reg_T.json', 'w'), indent=1)
