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


def _main():
    {'r1': r1, 'r23': r23, 't2': t2, 'tails': tails, 'twins': twins, 'landau': landau}[sys.argv[1]]()


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
