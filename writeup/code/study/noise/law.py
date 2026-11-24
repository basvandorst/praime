"""law.py: the class-count covariance law and the class-response chi^2 law (derivation.md), evaluated
(1) as closed forms, (2) retrospectively against the repo's reported anomalous statistics, (3) for the registration draft.
No repo imports. Reads (read-only) the repo's JSON results and this folder's simulation outputs:
  validate_cov.json (validate_cov.py), sl_*.json (simlaw.py), cr_1e9.json (classresp.py), check_T output (run separately).
Run: (ulimit -v 1500000; timeout 600 nice ../../../../.venv/bin/python -I law.py) > law_output.txt"""
import os, json, math
import numpy as np
from sympy import totient, factorint

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'data', 'study'))
EG = 0.5772156649015329
B = 1 - EG - math.log(2 * math.pi)                       # Montgomery-Soundararajan constant, -1.4148


def Sq(q):
    return sum(math.log(p) / (p - 1) for p in factorint(q))


def Lam(m):
    f = factorint(m)
    return math.log(list(f)[0]) if len(f) == 1 else 0.0


def v_class(q, H, lam):
    """Var(N_a)/E(N_a) for one reduced class a mod q, window length H, lam = ln x (eq. 4.3)."""
    return 1 - (math.log(H / q) + Sq(q) - B) / lam


def cov_ratio(q, c, lam):
    """Cov(N_a, N_{a+c})/E(N_a) for c != 0 mod q (eq. 4.4)."""
    m = q // math.gcd(c, q)
    return -Lam(m) / (int(totient(m)) * lam)


def contrast(q, H, lam):
    """E[chi^2]/dof of the equidistribution statistic over the phi(q) classes, total fitted (eq. 4.6)."""
    phi = int(totient(q))
    return v_class(q, H, lam) + (math.log(q) - Sq(q)) / ((phi - 1) * lam)


def cov_matrix(q, H, lam):
    red = [a for a in range(q) if math.gcd(a, q) == 1]
    n = len(red); C = np.zeros((n, n))
    for i, a in enumerate(red):
        for j, b in enumerate(red):
            C[i, j] = v_class(q, H, lam) if i == j else cov_ratio(q, b - a, lam)
    return red, C


def col_part(P, rho_cov):
    """E[column part] = tr(M^+ (I-P^T) Cov(delta) (I-P)), rows equal (uniform stationary law), per unit N/phi."""
    R = P.shape[0]
    M = sum(np.diag(P[a]) - np.outer(P[a], P[a]) for a in range(R))
    A = np.eye(R) - P.T
    return float(np.trace(np.linalg.pinv(M) @ A @ rho_cov @ A.T))


def load(name, base=HERE):
    return json.load(open(os.path.join(base, name)))


def section(t):
    print('\n' + '=' * 100 + '\n' + t + '\n' + '=' * 100)


def main():
    section('1. Closed forms: consistency identity sum_{b != a} Lambda(m)/phi(m) = ln q - S_q (needed for the total = MS)')
    for q in (3, 4, 5, 7, 9, 10, 12, 30, 210, 2310, 30030):
        red = [a for a in range(q) if math.gcd(a, q) == 1]
        s = sum(Lam(q // math.gcd(b - 1, q)) / int(totient(q // math.gcd(b - 1, q))) for b in red if b != 1)
        print('q=%6d  sum %.10f   ln q - S_q %.10f' % (q, s, math.log(q) - Sq(q)))

    section('2. Validation of the class-count law on primes, [5e11, 5e11+1e8) (validate_cov.py; not a test window)')
    V = load('validate_cov.json')
    print('%8s %4s | %-22s | %-22s | %-22s' % ('H', 'q', 'class Var/mean (law)', 'contrast chi2/dof (law)', 'total Var/mean (MS)'))
    for r in V['rows']:
        print('%8d %4d | %.3f +- %.3f (%.3f) | %.3f +- %.3f (%.3f) | %.3f +- %.3f (%.3f)' % (
            r['H'], r['q'], r['v_class_meas'], r['v_class_err'], r['v_class_law'], r['contrast_meas'], r['contrast_err'],
            r['contrast_law'], r['v_total_meas'], r['v_total_err'], r['v_total_law']))
    zc = [(r['contrast_meas'] - r['contrast_law']) / r['contrast_err'] for r in V['rows']]
    print('contrast z-scores:', ' '.join('%+.1f' % z for z in zc))

    section('3. RETROSPECTIVE (reported values known): spread of prime counts over the 48 classes mod 210 (r13f, round 36)')
    r13f = load('r13f.json', REPO); r36 = load('r36.json', REPO)
    wins = {'train': (1e7 + 4.5e7, 9e7), 'A': (1e10 + 2.5e7, 5e7), 'B': (1e12 + 2.5e7, 5e7), 'C': (1e15 + 1e7, 2e7)}
    print('window | measured | law (eq. 4.6) | round 36 (numerical HL) | sd (chi2_47) | z | sieve model 1-p | multinomial')
    for w, (xm, H) in wins.items():
        lam = math.log(xm); pred = contrast(210, H, lam); meas = r13f[w]['primes']; sd = pred * math.sqrt(2 / 47)
        print('%-6s | %.3f | %.3f | %.3f | %.3f | %+.1f | %.3f | 1' % (w, meas, pred, r36[w]['predicted_spread_HL'], sd,
                                                                   (meas - pred) / sd, 1 - (30030 / 5760) / lam))

    section('4. RETROSPECTIVE: last-digit tables, rounds 85-86 (Pearson chi^2 on 3 x 12 dof)')
    files = {'1e9': ('r85_1e9.json', 'sl_dig_9.json', 1e8), '1e12': ('r85_1e12.json', 'sl_dig_12.json', 2e9),
             '1e15': ('r86_1e15.json', 'sl_dig_15.json', 1e9)}
    print('scale k | chi2 meas | col part meas / within meas | law: chi2 (col, within) | sieve model M_Q | multinomial')
    for sc, (fr, fs, H) in files.items():
        r = load(fr, REPO); s = load(fs); lam = r['lambda']
        red, Cr = cov_matrix(10, H, lam)
        Cr = Cr - Cr.mean(1, keepdims=True) - Cr.mean(0, keepdims=True) + Cr.mean()   # contrasts only
        tot = {'meas': 0, 'law': 0, 'sieve': 0, 'col_meas': 0, 'col_law': 0}
        blocks_law = 0
        for k in ('1', '2', '3'):
            law = np.array(r['k'][k]['law']); dat = np.array(r['k'][k]['data']); N = np.array(r['k'][k]['row_counts'])
            Xc = dat * N[:, None]; E = law * N[:, None]
            chi = float(((Xc - E) ** 2 / E).sum())
            eps = Xc - E; c = eps.sum(0)
            M = sum(N[a] * (np.diag(law[a]) - np.outer(law[a], law[a])) for a in range(4))
            col = float(c @ np.linalg.pinv(M) @ c)
            st = s['stats']['trans:10:%s' % k]
            colL = col_part(law, Cr)                                     # law's column part (closed form, mod 10)
            print('%-5s %s | %6.2f | %5.2f / %5.2f | %5.2f (%4.2f, %5.2f) | %5.2f | 12' % (
                sc, k, chi, col, chi - col, st['chi2_law'], colL, st['chi2_law'] - colL, st['chi2_sieve_model']))
            tot['meas'] += chi; tot['law'] += st['chi2_law']; tot['sieve'] += st['chi2_sieve_model']
            tot['col_meas'] += col; tot['col_law'] += colL
            blocks_law = blocks_law + np.array(st['chi_blocks']) * st['chi2_law'] / st['chi2_sieve_model']
        sd = float(np.std(blocks_law, ddof=1))
        sds = float(np.std(sum(np.array(s['stats']['trans:10:%s' % k]['chi_blocks']) for k in ('1', '2', '3')), ddof=1))
        print('%-5s all | %6.2f | col %.2f | law %.2f +- %.2f (z %+.1f) | sieve %.2f +- %.2f (z %+.1f) | multinomial 36 (z %+.1f)' % (
            sc, tot['meas'], tot['col_meas'], tot['law'], sd, (tot['meas'] - tot['law']) / sd, tot['sieve'], sds,
            (tot['meas'] - tot['sieve']) / sds, (tot['meas'] - 36) / sds))
        print('        sum z^2 (48 binomial z): measured %.1f; law approx (4/3) x chi2_law = %.1f; multinomial 48' % (
            sum((np.array(r['k'][k]['z']) ** 2).sum() for k in ('1', '2', '3')), 4 / 3 * tot['law']))

    section('5. RETROSPECTIVE: block errors of the last-digit proportions (review 24: variance 0.75-0.87 of binomial)')
    for sc, (fr, fs, H) in files.items():
        r = load(fr, REPO); lam = r['lambda']; s = load(fs)
        for NB in (10, 20, 100):
            Hb = H / NB
            # per-table ratio at block length: scale FO_HL with v_Q(Hb)/v_Q(H) (approximation, cross terms small)
            vals = []
            for k in ('1', '2', '3'):
                st = s['stats']['trans:10:%s' % k]
                vQ_H = s['v_Q']; vQ_b = vQ_H + math.log(H / Hb) / lam
                vals.append((st['chi2_sieve_model'] - st['FO_sieve'] + st['FO_HL'] * vQ_b / vQ_H) / 12)
            print('%-5s NB=%3d (blocks of %.0e): predicted variance ratio per cell %.2f %.2f %.2f' % (sc, NB, Hb, *vals))

    section('6. RETROSPECTIVE: round 83, bands of D_k at [1e9, 1.1e9) (chi^2 on 19 dof; k=1 on 16)')
    r83 = load('r83.json', REPO); cr = load('sl_band_1e9.json')
    for k in ('1', '2', '3', '4'):
        d = r83['k'][k]; b = d['bands']
        lo = np.array([x['lo'] for x in b]); cnt = np.array([x['count'] for x in b], float); e = np.array([x['law'] for x in b])
        z = (cnt - e) / np.sqrt(e); P = e / e.sum(); cen = lo + 5.0
        u = np.sqrt(P) * (cen - (P * cen).sum()); u /= np.linalg.norm(u)
        st = cr['stats']['band:%s' % k]
        sdl = st['chi2_sd_per_block'] * st['chi2_law'] / st['chi2_sieve_model']
        print('k=%s measured %.1f (scale part %.2f) | law %.1f +- %.1f (z %+.1f) | sieve model %.1f +- %.1f, FO_M %.2f, FO_HL %.2f | nominal dof %d' % (
            k, d['law']['chi2'], (u @ z) ** 2, st['chi2_law'], sdl, (d['law']['chi2'] - st['chi2_law']) / sdl,
            st['chi2_sieve_model'], st['chi2_sd_per_block'], st['FO_sieve'], st['FO_HL'], d['law']['dof']))

    section('7. RETROSPECTIVE: round 13 residue x gap tables mod 210 (chi^2/dof, dof 401)')
    meas = {'train': 0.691, 'A': 0.857, 'B': 0.995, 'C': 0.894}
    r39 = load('r39.json', REPO)
    for w in ('train', 'A', 'B', 'C'):
        s = load('sl_r13_%s.json' % w)['stats']['resgap']
        print('%-5s measured %.3f | class-response law %.3f | round 39 (four-point HL) %.3f | sieve model %.3f | null 1.000 +- 0.058' % (
            w, meas[w], s['ratio_law'], r39[w]['implied_chi2_dof'], s['ratio_sieve']))

    section('8. REGISTRATION predictions (fresh windows; nothing measured here)')
    reg = {'W1': 3e11, 'W2': 6e11, 'W3': 9e11}
    Hb = 1e6; NB = 100
    print('S1: per-block equidistribution chi^2/(phi-1), blocks of 1e6, 100 per window')
    print('window q | law | sd (per window) | sieve model (1-p) | multinomial')
    for w, x0 in reg.items():
        lam = math.log(x0 + 5e7); p = (30030 / 5760) / lam
        for q in (3, 4, 10, 210):
            phi = int(totient(q)); L = contrast(q, Hb, lam)
            print('%s %4d | %.4f | %.4f | %.4f | 1' % (w, q, L, L * math.sqrt(2 / (phi - 1)) / math.sqrt(NB), 1 - p))
    runs = [load(f) for f in sorted(os.listdir(HERE)) if (f.startswith('sl_reg_W') or f.startswith('sl_regx_')) and f.endswith('.json')]
    specs = ['trans:3:1', 'trans:3:2', 'trans:10:1', 'trans:10:2']
    print('\nS2: block-homogeneity chi^2 of transition tables, 100 blocks of 1e6 per window; ratio = chi^2/((B-1) R(C-1))')
    print('short-range part E_M - FO_M pooled over %d simulation runs (%d blocks):' % (len(runs), sum(r['blocks_sim'] for r in runs)))
    base = {}
    for sp in specs:
        vals = np.array([r['stats'][sp]['chi2_sieve_model'] - r['stats'][sp]['FO_sieve'] for r in runs])
        errs = np.array([r['stats'][sp]['chi2_sieve_err'] for r in runs])
        w_ = 1 / errs ** 2; m_ = float((vals * w_).sum() / w_.sum()); e_ = float(1 / math.sqrt(w_.sum()))
        chi_runs = float((((vals - m_) / errs) ** 2).sum())
        base[sp] = (m_, e_)
        print('  %-10s E_M - FO_M = %.4f +- %.4f  (runs: %s; consistency chi2 %.1f on %d)' % (
            sp, m_, e_, ' '.join('%.3f' % v for v in vals), chi_runs, len(vals) - 1))
    tot = {sp: [] for sp in specs}
    for w in ('W1', 'W2', 'W3'):
        s_ = load('sl_reg_%s.json' % w); lam = s_['lam']; p = s_['p']
        print('\n%s  x0 = %.0e, lam = %.4f, p = %.4f, v_Q(1e6) incl. block term = %.4f' % (w, reg[w], lam, p, s_['v_Q']))
        print('  stat       | ratio: law (pred. err) [sd per window] | sieve M_Q | chain/multinomial | col/(C-1): law, M_Q, chain | within/((R-1)(C-1)): law')
        for sp in specs:
            st = s_['stats'][sp]; q = int(sp.split(':')[1]); R = int(totient(q)); dof = st['dof']
            chi_law = base[sp][0] + st['FO_HL']; ratio = chi_law / dof; perr = base[sp][1] / dof
            sd_blk = st['chi2_sd_per_block'] * chi_law / st['chi2_sieve_model']
            sdw = sd_blk * math.sqrt(NB) / ((NB - 1) * dof)
            P = np.array(st['P_pooled']).reshape(R, R)
            red, Cr = cov_matrix(q, Hb, lam)
            Cr = Cr + np.eye(R) * math.log(NB) / ((NB - 1) * lam)
            Cr = Cr - Cr.mean(1, keepdims=True) - Cr.mean(0, keepdims=True) + Cr.mean()
            colL = col_part(P, Cr)
            colS = col_part(P, (1 - p) * (np.eye(R) - np.ones((R, R)) / R))
            within = (chi_law - colL) / ((R - 1) ** 2)
            tot[sp].append((ratio, sdw, perr, colL / (R - 1), sd_blk))
            print('  %-10s | %.4f (+-%.4f) [%.4f] | %.4f | 1 | %.4f, %.4f (sim %.4f), 1 | %.4f' % (
                sp, ratio, perr, sdw, st['ratio_sieve'], colL / (R - 1), colS / (R - 1), st['col_part_sieve_per_block'] / (R - 1), within))
    print('\npooled over W1-W3 (mean of the three window ratios):')
    for sp in specs:
        a = np.array(tot[sp]); pr = a[:, 0].mean(); sd = math.sqrt((a[:, 1] ** 2).sum()) / 3
        print('  %-10s law %.4f  sd(measurement) %.4f  prediction err %.4f  | col/(C-1) law %.4f' % (sp, pr, sd, a[:, 2].mean(), a[:, 3].mean()))

if __name__ == '__main__':
    main()
