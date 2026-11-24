"""The within part of block transition tables of consecutive primes, from the Hardy-Littlewood conjecture.
Own code (imports only the own modules hlcore, chainlaw, triple, tailmod, components in this folder; nothing from the repo).
Stage 1 (components.py LAM ...) caches the HL components per height in data/comp_<lam>.json; this file assembles:
  (1) closed forms at leading order (from Lemke Oliver-Soundararajan's non-adjacent term):
        k = 1:  E[within]/dof = 1 - 1/((phi(q)-1) lam)
        k = 2:  E[within]/dof = 1 - 1/(3 (phi(q)-1) lam) + 2 B_q^2/(phi(q)-1)^2      (B_q = phi P(p_{n+1}=p_n) - 1)
  (2) the finite-lam law: W(H) = 1 + 2 sum_j (1 - j/N) tau_j + W_tail(H), N = H/lam, with
        tau_1 (k = 1) from the HL triple law, tau_j (j >= 2) from the fine-class chain, W_tail from the insertion responses
        (k = 2: tau_2 from the leading-order LOS term -1/(6 (C-1) lam), the rest as for k = 1);
  (3) the retrospective table against the nine measured windows (values copied from r135c.json / r135b_*.json),
  (4) the validation-stretch components, (5) the predictions for the registration windows.
usage: law.py  (after components.py has run for the six heights listed in LAMS_OLD and LAMS_NEW)"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, math, json, glob
HERE = os.path.dirname(os.path.abspath(__file__))
G = 0.5772156649015329
OLD = [2e11, 3e11, 4e11, 5.5e11, 6e11, 7e11, 8e11, 9e11, 1e12]
NEW = [6.0e9 + i * 1e8 for i in range(24)]
lam_of = lambda x0: math.log(x0 + 5e7)
HS = [10 ** 4, 10 ** 5, 10 ** 6]
# measured, nine windows (r135c.json 'all', win_q_k_H): mean, pooled block error, window-scatter error (std/3 of per_window)
MEAS = {(3, 1, 10 ** 4): (0.9549, 0.0045, 0.0062), (3, 1, 10 ** 5): (0.9562, 0.0144, 0.0087), (3, 1, 10 ** 6): (0.9611, 0.0442, 0.0352),
        (5, 1, 10 ** 4): (0.9844, 0.0015, 0.0019), (5, 1, 10 ** 5): (0.9728, 0.0048, 0.0038), (5, 1, 10 ** 6): (0.9486, 0.0152, 0.0070),
        (7, 1, 10 ** 4): (0.9859, 0.0009, 0.0016), (7, 1, 10 ** 5): (0.9798, 0.0029, 0.0027), (7, 1, 10 ** 6): (0.9735, 0.0089, 0.0055),
        (3, 2, 10 ** 4): (1.0047, 0.0047, 0.0043), (3, 2, 10 ** 5): (1.0174, 0.0157, 0.0114), (3, 2, 10 ** 6): (0.9555, 0.0465, 0.0552),
        (5, 2, 10 ** 4): (1.0060, 0.0016, 0.0016), (5, 2, 10 ** 5): (0.9980, 0.0049, 0.0048), (5, 2, 10 ** 6): (0.9938, 0.0160, 0.0112)}
# (5, 2) is the table mod 10 (the same table for odd primes)
def phi_q(q): return sum(1 for a in range(1, q) if math.gcd(a, q) == 1)
def los_c2_aa(q):
    ps = [p for p in (2, 3, 5, 7, 11, 13) if q % p == 0]; f = phi_q(q)
    return (f * math.log(q / (2 * math.pi)) + math.log(2 * math.pi)) / 2 - f / 2 * sum(math.log(p) / (p - 1) for p in ps)
def los_Bq(q, lam):
    f = phi_q(q); return (-(f - 1) / 2 * math.log(lam) + los_c2_aa(q)) / lam
def lead_k1(q, lam): return 1 - 1 / ((phi_q(q) - 1) * lam)
def lead_k2(q, lam): return 1 - 1 / (3 * (phi_q(q) - 1) * lam)
def lead_k2_wick(q, lam): return lead_k2(q, lam) + 2 * los_Bq(q, lam) ** 2 / (phi_q(q) - 1) ** 2

COMP = {}
for f in glob.glob(os.path.join(HERE, 'data', 'comp_*.json')):
    d = json.load(open(f)); COMP[d['lam']] = d
def comp_at(lam, region):
    """components interpolated (quadratic in 1/lam) from the cached heights of one region"""
    ls = sorted(l for l in COMP if (l > 24.5) == (region == 'old'))
    return ls
def W_parts(d, q, k, H, variant='law'):
    lam = d['lam']; N = H / lam; C = phi_q(q); c = d[str(q)]
    j = np.arange(61); wj = np.clip(1 - j / N, 0, None)
    if variant == 'MQ':
        tau = np.array(c['MQ_k%d' % k]); return 1 + 2 * float(np.sum(wj[1:] * tau[1:]))
    if k == 1:
        tau = np.array(c['chain_k1']).copy()
        t1 = d['tau1_triple_y13'][str(q)] if variant != 'y11' else d['tau1_triple_y11'][str(q)]
        tau[1] = t1
        tail = c['tail_k1'][str(H)] if variant not in ('y11', 'w3') else (c['tail_k1_y11'][str(H)] if variant == 'y11' else c['tail_k1_w3'][str(H)])
        if variant == 'notail': tail = 0.0
    else:
        tau = np.array(c['chain_k2']).copy(); tau[2] = -1 / (6 * (C - 1) * lam)
        tail = c['tail_k2'][str(H)] if variant != 'notail' else 0.0
    return 1 + 2 * float(np.sum(wj[1:] * tau[1:])) + tail
def interp(lams_have, vals, lam):
    x = 1 / np.array(lams_have); y = np.array(vals); p = np.polyfit(x, y, min(2, len(x) - 1)); return float(np.polyval(p, 1 / lam))
def pooled(windows, region, q, k, H, variant='law'):
    ls = comp_at(None, region)
    vals = [W_parts(COMP[l], q, k, H, variant) for l in ls]
    return float(np.mean([interp(ls, vals, lam_of(x0)) for x0 in windows]))
def law_err(windows, region, q, k, H):
    """the law's own uncertainty: truncation (y = 13 vs 11), the w2/w3 choice, and +-40 % of the tail (validation
    stretch: the model tail is 50-70 % of the measured tail bins for q = 5, 7)"""
    base = pooled(windows, region, q, k, H)
    e2 = 0.0
    if k == 1:
        e2 += (pooled(windows, region, q, k, H, 'y11') - base) ** 2 + (pooled(windows, region, q, k, H, 'w3') - base) ** 2
    tail = base - pooled(windows, region, q, k, H, 'notail')
    e2 += (0.4 * tail) ** 2
    if k == 2:      # tau_2 from the leading-order LOS term: validation shows +-0.001 per direction-sum -> 0.002 in W
        e2 += 0.002 ** 2
    return math.sqrt(e2)

if __name__ == '__main__':
    print('cached heights:', ', '.join('%.4f' % l for l in sorted(COMP)))
    lo = np.mean([lam_of(x) for x in OLD]); ln = np.mean([lam_of(x) for x in NEW])
    print('\n(1) closed forms, leading order (LOS non-adjacent term), at the mean heights lam = %.3f (old) and %.3f (new)' % (lo, ln))
    for q in (3, 5, 7):
        print('  q=%d  k=1: 1 - 1/((phi-1) lam) = %.4f / %.4f    k=2: 1 - 1/(3(phi-1)lam) = %.4f / %.4f,  + 2 B_q^2/(phi-1)^2 (LOS B_q = %+.4f / %+.4f) = %.4f / %.4f'
              % (q, lead_k1(q, lo), lead_k1(q, ln), lead_k2(q, lo), lead_k2(q, ln), los_Bq(q, lo), los_Bq(q, ln), lead_k2_wick(q, lo), lead_k2_wick(q, ln)))
    print('\n(2)-(3) RETROSPECTIVE (the measured values were known before the law was written): nine windows 2e11-1e12')
    print('  cell (q,k,H)     measured +- err(block; window scatter) | law +- law err | z      | lead.order  z  | no tail | sieve chain M_Q | 1')
    rows = []
    for (q, k, H), (m, eb, ew) in MEAS.items():
        L = pooled(OLD, 'old', q, k, H); le = law_err(OLD, 'old', q, k, H)
        lead = float(np.mean([lead_k1(q, lam_of(x)) if k == 1 else lead_k2_wick(q, lam_of(x)) for x in OLD]))
        nt = pooled(OLD, 'old', q, k, H, 'notail'); mq = pooled(OLD, 'old', q, k, H, 'MQ')
        e = max(eb, ew); z = (m - L) / math.sqrt(e ** 2 + le ** 2); zl = (m - lead) / e
        rows.append(dict(q=q, k=k, H=H, meas=m, err=e, law=L, law_err=le, z=z, lead=lead, z_lead=zl, notail=nt, MQ=mq))
        print('  (%d,%d,%.0e)   %.4f +- %.4f (%.4f; %.4f) | %.4f +- %.4f | %+5.1f | %.4f %+5.1f | %.4f | %.4f | 1 (z %+.1f)'
              % (q, k, H, m, e, eb, ew, L, le, z, lead, zl, nt, mq, (m - 1) / e))
    json.dump(rows, open(os.path.join(HERE, 'data', 'retro.json'), 'w'), indent=1)
    print('\n  slopes in H (law, old windows): ' + '  '.join('(%d,%d): %+.4f/%+.4f' % (q, k, pooled(OLD, 'old', q, k, 10 ** 5) - pooled(OLD, 'old', q, k, 10 ** 4), pooled(OLD, 'old', q, k, 10 ** 6) - pooled(OLD, 'old', q, k, 10 ** 5)) for q, k in ((3, 1), (5, 1), (7, 1), (3, 2), (5, 2))))
    print('  measured slopes:                ' + '  '.join('(%d,%d): %+.4f/%+.4f' % (q, k, MEAS[(q, k, 10 ** 5)][0] - MEAS[(q, k, 10 ** 4)][0], MEAS[(q, k, 10 ** 6)][0] - MEAS[(q, k, 10 ** 5)][0]) for q, k in ((3, 1), (5, 1), (7, 1), (3, 2), (5, 2))))
    # (4) validation stretch
    vm = json.load(open(os.path.join(HERE, 'data', 'val_measure.json'))); vl = min(COMP, key=lambda l: abs(l - vm['lam']))
    d = COMP[vl]
    print('\n(4) VALIDATION STRETCH [1.5e11, 1.5e11 + 1e8) (lam %.4f): components, model against primes' % vm['lam'])
    for q in (3, 5, 7):
        C = phi_q(q); tm = vm['tables']['%d_1' % q]['tau']; tm2 = vm['tables']['%d_2' % q]['tau']
        print('  q=%d k=1: tau_1 triple-law %+.5f (y=11: %+.5f; chain only %+.5f; LOS lead %+.5f) vs measured %+.5f +- %.5f;  tau_2 chain %+.5f vs %+.5f +- %.5f'
              % (q, d['tau1_triple_y13'][str(q)], d['tau1_triple_y11'][str(q)], d[str(q)]['chain_k1'][1], -1 / (2 * (C - 1) * vm['lam']),
                 tm['1'][0] / (C - 1) ** 2, tm['1'][1] / (C - 1) ** 2, d[str(q)]['chain_k1'][2], tm['2'][0] / (C - 1) ** 2, tm['2'][1] / (C - 1) ** 2))
        print('       k=2: tau_1 chain %+.5f (Wick B_q^2/(C-1)^2 = %+.5f) vs %+.5f +- %.5f;  tau_2 LOS %+.5f (chain %+.5f) vs %+.5f +- %.5f;  model B_q %+.5f'
              % (d[str(q)]['chain_k2'][1], d[str(q)]['Bq'] ** 2 / (C - 1) ** 2, tm2['1'][0] / (C - 1) ** 2, tm2['1'][1] / (C - 1) ** 2,
                 -1 / (6 * (C - 1) * vm['lam']), d[str(q)]['chain_k2'][2], tm2['2'][0] / (C - 1) ** 2, tm2['2'][1] / (C - 1) ** 2, d[str(q)]['Bq']))
    # (5) predictions at the registration windows
    print('\n(5) PREDICTIONS for the registration: 24 windows [x0, x0 + 1e8), x0 = 6.0e9, 6.1e9, ..., 8.3e9 (lam %.3f-%.3f)' % (lam_of(NEW[0]), lam_of(NEW[-1])))
    pred = {}
    for (q, k, H) in MEAS:
        L = pooled(NEW, 'new', q, k, H); le = law_err(NEW, 'new', q, k, H)
        lead = float(np.mean([lead_k1(q, lam_of(x)) if k == 1 else lead_k2_wick(q, lam_of(x)) for x in NEW]))
        mq = pooled(NEW, 'new', q, k, H, 'MQ'); old_law = pooled(OLD, 'old', q, k, H)
        pred['%d_%d_%d' % (q, k, H)] = dict(law=L, law_err=le, lead=lead, MQ=mq, old_meas=MEAS[(q, k, H)][0], old_law=old_law)
        print('  (%d,%d,%.0e): law %.4f +- %.4f   leading order %.4f   sieve chain %.4f   old-height law %.4f, old measured %.4f'
              % (q, k, H, L, le, lead, mq, old_law, MEAS[(q, k, H)][0]))
    json.dump(pred, open(os.path.join(HERE, 'data', 'predictions.json'), 'w'), indent=1)
