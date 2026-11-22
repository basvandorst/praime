"""Print the markdown tables of STUDY.md from the verdict files (no number in the tables is typed by hand)."""
import json, os, sys
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
J = lambda f: json.load(open(os.path.join(D, f)))

NAMES = {
    'S.gap2': 'gap 2', 'S.gap4': 'gap 4', 'S.gap6': 'gap 6', 'S.gap8': 'gap 8', 'S.gap10': 'gap 10',
    'S.gap12': 'gap 12', 'S.gap18': 'gap 18', 'S.gap30': 'gap 30', 'N.eqgap': 'two equal gaps in a row',
    'N.updown': 'P(next gap larger) − P(smaller)',
    'A.same3': 'next prime same class mod 3', 'A.same4': 'same class mod 4', 'A.same5': 'same class mod 5',
    'A.same7': 'same class mod 7', 'A.same8': 'same class mod 8',
    'S.var50': 'Var/mean, H = 50', 'S.var500': 'Var/mean, H = 500', 'S.var5000': 'Var/mean, H = 5000',
    'I.migap_lag1': 'mutual information of consecutive gaps (bits)', 'D.gain1_gap': 'Markov gain, order 1 (bits)',
    'G.chi2_BH': 'χ²/curve on Ulam diagonals', 'G.slope_BH': 'Bateman–Horn slope on Ulam diagonals',
    'T4.H_class': 'entropy of the gap class (bits)', 'T4.H_cond': 'conditional entropy H(next | gap) (bits)',
    'T3.cpap4': 'four primes in a row in arithmetic progression',
    'A.race4': 'race mod 4 (slope)', 'A.race3': 'race mod 3 (slope)', 'A.race8': 'race mod 8 (slope)',
}


def nm(s):
    if s.startswith('T1.tail'):
        return 'gaps ≥ %s' % s[7:]
    return NAMES.get(s, s)


def r1():
    v = J('verdicts_r1.json')
    print('| statistic | primes (train) | round-1 formula | rel. diff | z train | z A | z B |')
    print('|---|---|---|---|---|---|---|')
    for s, d in v['train'].items():
        if 'T' not in d or s.startswith('A.race') or s.startswith('G.corr'):
            continue
        z = [v[w][s].get('zT') for w in ('train', 'A', 'B')]
        rel = '%+.1f %%' % (100 * d['rel']) if d.get('rel') is not None else ''
        print('| %s | %.5g | %.5g | %s | %s |' % (nm(s), d['x'], d['T'], rel,
              ' | '.join('%.1f' % x if x is not None else '' for x in z)))


def r23():
    v2 = J('verdicts_theory_r2.json'); v3 = J('verdicts_r3.json'); b3 = J('bias_r3.json'); g = J('verdicts_theory_geo.json')
    stats = [s for s in v3['train'] if not s.startswith('T2') and not s.startswith('T1')]
    print('| statistic | primes (train) | cluster law | rel. diff train | z train | z A | z B | z C (10¹⁵) |')
    print('|---|---|---|---|---|---|---|---|')
    for s in stats:
        d = v3['train'][s]
        zs = []
        for w in ('train', 'A', 'B', 'C'):
            if s not in v3[w]:
                zs.append(''); continue
            z = v3[w][s]['zT']
            if s in b3[w]:
                z = b3[w][s]['z_corrected']
            zs.append('%.1f' % z)
        print('| %s | %.5g | %.5g | %+.2f %% | %s |' % (nm(s), d['x'], d['T'], 100 * d['rel'], ' | '.join(zs)))
    s = 'G.chi2_BH'
    print('| %s | %.4g | %.4g | %+.1f %% | %s | |' % (nm(s), g['train'][s]['x'], g['train'][s]['T'], 100 * g['train'][s]['rel'],
          ' | '.join('%.1f' % g[w][s]['zT'] for w in ('train', 'A', 'B'))))


def tails():
    v3 = J('verdicts_r3.json')
    print('| window | ln x | gaps ≥ 3 ln x: primes / law, z | ≥ 4 ln x | ≥ 5 ln x |')
    print('|---|---|---|---|---|')
    import math
    lm = {'train': math.log(5.5e7), 'A': math.log(1e10), 'B': math.log(1e12), 'C': math.log(1e15)}
    for w in ('train', 'A', 'B', 'C'):
        ts = [k for k in v3[w] if k.startswith('T1')]
        cells = ['h ≥ %s: %.4g / %.4g, %.1f' % (k[7:], v3[w][k]['x'], v3[w][k]['T'], v3[w][k]['zT']) for k in ts]
        print('| %s | %.1f | %s |' % (w, lm[w], ' | '.join(cells)))


def t2():
    v3 = J('verdicts_r3.json')
    print('| window | gap pairs in the cells | χ² | cells | χ²/cell | z |')
    print('|---|---|---|---|---|---|')
    for w in ('train', 'A', 'B', 'C'):
        d = v3[w]['T2.chi2_joint']
        print('| %s | %d | %.1f | %d | %.3f | %.1f |' % (w, d['sum_obs'], d['chi2'], d['dof'], d['chi2_per_dof'], d['z']))


def landau():
    v = J('verdicts_landau.json')
    lab = lambda s: {'L.ge60': 'n-gap ≥ 60', 'L.ge100': 'n-gap ≥ 100'}.get(s, ('n-gap %s' % s[5:]) if s.startswith('L.gap')
                                                                       else 'same n mod %s' % s[6:])
    print('| statistic | F1: data / law | t₉, z_bin | F2: data / law | t₉, z_bin | F3: data / law | t₉, z_bin |')
    print('|---|---|---|---|---|---|---|')
    for s in v['F1']:
        cells = []
        for w in ('F1', 'F2', 'F3'):
            d = v[w][s]
            cells.append('%.4g / %.4g | %.1f, %.1f' % (d['x'], d['T'], d['zT'], d['z_binom']))
        print('| %s | %s |' % (lab(s), ' | '.join(cells)))
    sm = v['_summary']
    print()
    print('%d tests; %d with |t₉| > 3 (expected by chance %.1f); largest |z_bin| = %.1f.' % (
        sm['tests'], sm['misses'], sm['expected_by_chance'],
        max(abs(v[w][s]['z_binom']) for w in ('F1', 'F2', 'F3') for s in v[w])))


def bits():
    b = J('bits_r6.json')
    keys = list(b['train'].keys())
    print('| model or measurement | ' + ' | '.join(['train', 'A', 'B', 'C (10¹⁵)']) + ' |')
    print('|---|---|---|---|---|')
    for k in keys:
        print('| %s | %s |' % (k, ' | '.join('%.4f' % b[w][k] for w in ('train', 'A', 'B', 'C'))))


def chain():
    r = J('r7.json')
    print('| window | measured H₁ / H₂ | H₁ − H₂ measured | mod 6: H₁ / H₂, H₁ − H₂ | mod 30 | mod 210 |')
    print('|---|---|---|---|---|---|')
    for w in ('train', 'A', 'B', 'C'):
        m = r['%s_Q6' % w]
        cells = []
        for Q in (6, 30, 210):
            d = r['%s_Q%d' % (w, Q)]
            cells.append('%.4f / %.4f, %.4f' % (d['H1_pred'], d['H2_pred'], d['H1_pred'] - d['H2_pred']))
        print('| %s | %.4f / %.4f | %.4f | %s |' % (w, m['H1_meas'], m['H2_meas'], m['H1_meas'] - m['H2_meas'], ' | '.join(cells)))


def gauss():
    v = J('verdicts_r11_dd.json'); v.update(J('verdicts_r11_rep_dd.json'))
    boxes = [('G1', 'a, b ∈ [10⁵, 10⁵ + 2000)'), ('G2', 'a, b ∈ [5·10⁶, 5·10⁶ + 1000)'),
             ('G3', 'a, b ∈ [2·10⁵, 2·10⁵ + 2000)'), ('G4', 'a ∈ [10⁵, +2000), b ∈ [3·10⁵, +2000)')]
    print('| statistic | ' + ' | '.join('%s: data / law, t₉' % b for b, _ in boxes) + ' |')
    print('|---|---|---|---|---|')
    for s in v['G1']:
        print('| %s | %s |' % (s.replace('d=', 'd = ').replace('d>=', 'd ≥ ').replace('mean_d', 'mean d'),
              ' | '.join('%.4g / %.4g, %.1f' % (v[b][s]['x'], v[b][s]['T'], v[b][s]['zT']) for b, _ in boxes)))


def r9():
    r = J('r9.json'); r.update(J('r9_U159.json'))
    names = {'U': 'Ulam diagonals 4n² + bn + c (all 198)', 'U159': 'Ulam diagonals, the 159 distinct curves',
             'E': "Euler's n² + n + c, c odd ≤ 199"}
    print('| family | curves | prime values | χ²/curve | law | z | Bateman–Horn: mean (obs − E)/E, z |')
    print('|---|---|---|---|---|---|---|')
    for k in ('U', 'U159', 'E'):
        d = r[k]
        print('| %s | %d | %s | %.3f | %.3f | %.1f | %+.5f, %.1f |' % (names[k], d['curves'], format(d['primes'], ','),
              d['chi2'], d['theory'], d['z'], d['mean_rel_dev_from_BH'], d['z_BH']))


def sg():
    v = J('verdicts_r14.json')
    lab = lambda s: ('gap %s' % s[6:]) if s.startswith('W1') else {'W2.le300': 'gap ≤ 300',
          'W3.same5': 'same class mod 5 (gaps ≤ 300)', 'W3.same7': 'same class mod 7 (gaps ≤ 300)'}[s]
    print('| statistic | SG primes (train) | law, order 1 / 3 | t₉ train | t₉ A | t₉ B | t₉ C |')
    print('|---|---|---|---|---|---|---|')
    for s, d in v['train'].items():
        print('| %s | %.5g | %.5g / %.5g | %s |' % (lab(s), d['x'], d['order1'], d['T'],
              ' | '.join('%.1f' % v[w][s]['zT'] for w in ('train', 'A', 'B', 'C'))))


def eis():
    v = J('verdicts_r16_dd.json'); v.update(J('verdicts_r16_rep_dd.json'))
    boxes = ['E1', 'E2', 'E3', 'E4']
    print('| statistic | ' + ' | '.join('%s: data / law, t₉' % b for b in boxes) + ' |')
    print('|---|---|---|---|---|')
    for s in v['E1']:
        print('| %s | %s |' % (s.replace('d=', 'd = ').replace('d>=', 'd ≥ ').replace('mean_d', 'mean d'),
              ' | '.join('%.4g / %.4g, %.1f' % (v[b][s]['x'], v[b][s]['T'], v[b][s]['zT']) for b in boxes)))


def r17():
    r = J('r17.json')
    print('| statistic | primes near 10¹⁸ | cluster law | rel. diff | t₉ | z (binomial) |')
    print('|---|---|---|---|---|---|')
    for s, d in r['stats'].items():
        print('| %s | %.5g | %.5g | %+.2f %% | %.1f | %s |' % (nm(s), d['x'], d['T'], 100 * d['rel'], d['t9'],
              '%.1f' % d['z_binom'] if d['z_binom'] is not None else ''))


def scoreboard():
    from scipy.stats import t as st, norm
    p9 = float(2 * st.sf(3, 9)); pn = float(2 * norm.sf(3))
    rows = []

    def summ(f, label, rnd):
        v = J(f)
        if '_summary' in v:
            sm = v['_summary']; rows.append((rnd, label, sm['tests'], sm['misses'], sm['tests'] * p9))
        else:
            zs = [d.get('zT') for w in v for d in v[w].values() if isinstance(d, dict) and d.get('zT') is not None]
            rows.append((rnd, label, len(zs), sum(abs(z) > 3 for z in zs), len(zs) * p9))
    summ('verdicts_r3_t9.json', 'gap law, classes, variance, information, tails, CPAP-4 (4 windows; registered test)', '2-3')
    summ('verdicts_twins.json', 'consecutive twin primes', '4')
    summ('verdicts_landau.json', 'consecutive n with n² + 1 prime', '5')
    summ('verdicts_r6.json', 'second differences', '6')
    g1 = J('verdicts_r11_dd.json'); g2 = J('verdicts_r11_rep_dd.json')
    zs = [d['zT'] for v in (g1, g2) for w in v if not w.startswith('_') for d in v[w].values()]
    rows.append(('11', 'Gaussian primes, nearest neighbours (4 boxes)', len(zs), sum(abs(z) > 3 for z in zs), len(zs) * p9))
    summ('verdicts_r14.json', 'consecutive Sophie Germain primes', '14')
    ve = J('verdicts_euler41.json')
    cnt = {json.loads(l)['window']: json.loads(l)['stats']['count'] for l in open(os.path.join(D, 'results_euler41.jsonl'))
           if json.loads(l)['block'] == -1}
    zs = [d['zT'] for w in ('F1', 'F2', 'F3') for d in ve[w].values() if d['T'] * (cnt[w] - 1) >= 5]
    rows.append(('18', 'consecutive n with n² + n + 41 prime (tests with expected count ≥ 5)', len(zs), sum(abs(z) > 3 for z in zs), len(zs) * p9))
    e1 = J('verdicts_r16_dd.json'); e2 = J('verdicts_r16_rep_dd.json')
    zs = [d['zT'] for v in (e1, e2) for w in v if not w.startswith('_') for d in v[w].values()]
    rows.append(('16', 'Eisenstein primes, nearest neighbours (4 boxes)', len(zs), sum(abs(z) > 3 for z in zs), len(zs) * p9))
    r = J('r17.json')
    rows.append(('17', 'the law at 10¹⁸', len(r['stats']), sum(abs(d['t9']) > 3 for d in r['stats'].values()), len(r['stats']) * p9))
    sk = J('r10_skew.json')
    zs = [sk[w][H]['t9'] for w in ('train', 'A', 'B', 'C') for H in ('50', '500')]
    rows.append(('10', 'third cumulant of interval counts', len(zs), sum(abs(z) > 3 for z in zs), len(zs) * p9))
    v1 = J('r10_var.json'); v2 = J('r10_var_3000000000.json')
    zs = [d[H]['t9'] for d in (v1, v2) for H in ('10000', '100000', '1000000')]
    rows.append(('10', 'variance in long intervals (2 windows)', len(zs), sum(abs(z) > 3 for z in zs), len(zs) * p9))
    r9 = J('r9.json')
    zs = [r9[k]['z'] for k in ('U', 'E')] + [r9[k]['z_BH'] for k in ('U', 'E')]
    rows.append(('9', 'variance and means along quadratics (fresh data)', len(zs), sum(abs(z) > 3 for z in zs), len(zs) * pn))
    r15 = J('r15.json')
    zs = [r15[w]['z'] for w in r15]
    rows.append(('15', 'joint law of three gaps, 729 cells (4 windows)', len(zs), sum(abs(z) > 3 for z in zs), len(zs) * pn))
    import numpy as np
    r20rows = {}
    r20rows.update(J('r20_H30.json'))
    for w in ('train', 'A', 'B', 'C'):
        r20rows.update(J('r20_H60_%s.json' % w))
    zs = [t for r in r20rows.values() for t, l in zip(r['t9'], r['law']) if l >= 1e-3]
    rows.append(('20', 'count distribution in short intervals (H = 30, 60; cells with P ≥ 10⁻³)', len(zs), sum(abs(z) > 3 for z in zs), len(zs) * p9))
    r21 = J('r21.json')
    zs = [t for w in ('train', 'A', 'B', 'C') for t, l in zip(r21[w]['t9'], r21[w]['law_exact']) if l >= 1e-3]
    r21b = J('r21b.json')
    zs += [t for w in ('train', 'A', 'B', 'C') for t, l in zip(r21b[w]['t9'], r21b[w]['law_exact']) if l >= 1e-3]
    rows.append(('21', 'the same with the exact Hardy–Littlewood sums (H = 30 and 60)', len(zs), sum(abs(z) > 3 for z in zs), len(zs) * p9))
    a70 = J('r22_HM70.json'); a52 = J('r22_HM52.json')
    zs = [d['t9'] for d in a70['stats'].values()] + [a52['stats']['T1.tail54']['t9']]
    for w in ('A', 'B', 'C'):
        zs += [d['t9'] for d in J('r22_HM70_%s.json' % w)['stats'].values()]
    rows.append(('22', 'exact gap law: gaps 2–70 and the tail ≥ 72 in four windows (and ≥ 54 at train)', len(zs), sum(abs(z) > 3 for z in zs), len(zs) * p9))
    r24 = J('r24.json')['U_summary']
    zs = [(r24['mean'] - 1) / r24['se'], 0.0 if r24['ks_p'] > 0.01 else 99.0]
    rows.append(('24', 'first occurrences of the gaps 180–454 (mean of U, KS test against Exp(1))', 2, sum(abs(z) > 3 for z in zs), 0.0027 + 0.01))
    r13 = J('r13.json')
    zs = [d['z'] for d in r13.values()]
    rows.append(('13', 'residue chain mod 30 and 210, cell by cell (χ² too low counts as a miss)', len(zs), sum(abs(z) > 3 for z in zs), len(zs) * pn))
    rows.sort(key=lambda r: int(r[0].split('-')[0]))
    print('| round | what | tests | misses | expected by chance |')
    print('|---|---|---|---|---|')
    T = M = E = 0
    for rnd, lab, n, m, e in rows:
        print('| %s | %s | %d | %d | %.1f |' % (rnd, lab, n, m, e)); T += n; M += m; E += e
    print('| | **all** | **%d** | **%d** | **%.1f** |' % (T, M, E))


def euler():
    v = J('verdicts_euler41.json')
    lab = lambda s: {'L.ge15': 'n-gap ≥ 15'}.get(s, ('n-gap %s' % s[5:]) if s.startswith('L.gap') else 'same n mod %s' % s[6:])
    print('| statistic | F1: data / law | t₉, z_bin | F2: data / law | t₉, z_bin | F3: data / law | t₉, z_bin |')
    print('|---|---|---|---|---|---|---|')
    for s in v['F1']:
        cells = []
        for w in ('F1', 'F2', 'F3'):
            d = v[w][s]
            t = '%.1f' % d['zT'] if abs(d['zT']) < 100 else '(degenerate)'
            cells.append('%.4g / %.4g | %s, %.1f' % (d['x'], d['T'], t, d['z_binom']))
        print('| %s | %s |' % (lab(s), ' | '.join(cells)))


def r20():
    import numpy as np
    rows = {}
    rows.update(J('r20_H30.json'))
    for w in ('train', 'A', 'B', 'C'):
        rows.update(J('r20_H60_%s.json' % w))
    print('| window | H | P(0): data / law / Poisson | P(2): data / law / Poisson | distance to the data: law / Poisson | block misses (cells with P ≥ 10⁻³) |')
    print('|---|---|---|---|---|---|')
    for H in (30, 60):
        for w in ('train', 'A', 'B', 'C'):
            r = rows['%s_H%d' % (w, H)]
            d, l, p, t = map(np.array, (r['data'], r['law'], r['poisson'], r['t9']))
            tvl = 0.5 * np.abs(d - l).sum(); tvp = 0.5 * np.abs(d - p).sum()
            sel = l >= 1e-3
            print('| %s | %d | %.4f / %.4f / %.4f | %.4f / %.4f / %.4f | %.4f / %.4f | %d of %d |' % (
                w, H, d[0], l[0], p[0], d[2], l[2], p[2], tvl, tvp, int(np.sum(np.abs(t[sel]) > 3)), int(sel.sum())))


def r21():
    print('| H | window | k = 0 | 1 | 2 | 3 | 4 | 5 | block misses |')
    print('|---|---|---|---|---|---|---|---|---|')
    for H, f in ((30, 'r21.json'), (60, 'r21b.json')):
        r = J(f)
        for w in ('train', 'A', 'B', 'C'):
            d = r[w]
            miss = sum(abs(t) > 3 for t, l in zip(d['t9'], d['law_exact']) if l >= 1e-3)
            n = sum(l >= 1e-3 for l in d['law_exact'])
            print('| %d | %s | %s | %d of %d |' % (H, w, ' | '.join('%.4f / %.4f' % (d['data'][k], d['law_exact'][k]) for k in range(6)), miss, n))


def r23():
    r = J('r23.json')
    print('| gap g | N(g) below 4·10¹⁸ (counted) | exact law | relative difference | (count − law)/√count |')
    print('|---|---|---|---|---|')
    for g, d in sorted(r.items(), key=lambda kv: int(kv[0])):
        print('| %s | %s | %s | %+.1e | %+.1f |' % (g, format(d['data'], ','), format(round(d['law']), ','), d['rel'], d['z_sqrtN']))


def _main():
    {'r1': r1, 'r23': r23, 't2': t2, 'tails': tails, 'twins': twins, 'landau': landau, 'bits': bits, 'chain': chain, 'gauss': gauss, 'r9': r9, 'sg': sg, 'eis': eis, 'r17': r17, 'scoreboard': scoreboard, 'euler': euler, 'r20': r20, 'r21': r21, 'r23': r23}[sys.argv[1]]()


def twins():
    v = J('verdicts_twins.json')
    print('| statistic | twins (train) | law, order 1 / 3 | z train | z A | z B | z C |')
    print('|---|---|---|---|---|---|---|')
    lab = lambda s: ('twin gap %s' % s[6:]) if s.startswith('W1') else {'W2.le300': 'twin gap ≤ 300',
          'W3.same5': 'same class mod 5 (gaps ≤ 300)', 'W3.same7': 'same class mod 7 (gaps ≤ 300)'}[s]
    for s, d in v['train'].items():
        print('| %s | %.5g | %.5g / %.5g | %s |' % (lab(s), d['x'], d['order1'], d['T'],
              ' | '.join('%.1f' % v[w][s]['zT'] for w in ('train', 'A', 'B', 'C'))))


if __name__ == '__main__':
    _main()
