"""Consistency check of STUDY.md against the data files.

1. Every table that tables.py prints must appear verbatim in STUDY.md.
2. Every number quoted in the prose (list PROSE below) must equal the value recomputed from the data files,
   rounded as quoted.
3. Small exact identities: S({0,2}) = 2 C2, S_a averages to S, S_f({0}) against the Bateman-Horn constant of n^2+1.
Exit status 1 on the first failure class with all failures printed.
"""
import io, json, os, sys, contextlib, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..', '..', '..')
D = os.path.join(HERE, '..', '..', 'data', 'study')
sys.path.insert(0, HERE)
import tables

study = open(os.path.join(ROOT, 'STUDY.md')).read()
J = lambda f: json.load(open(os.path.join(D, f)))
fail = []

# 1. tables
for name in ('r1', 'r3', 't2', 'tails', 'twins', 'landau', 'bits', 'chain', 'gauss', 'r9', 'sg', 'eis', 'r17', 'scoreboard', 'euler', 'r20', 'r21', 'r23'):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        getattr(tables, name)()
    for line in buf.getvalue().strip().splitlines():
        if line.strip() and line.replace('|', '\\|') not in study and line not in study:
            # the r23 table escapes the | inside H(next | gap)
            if line.replace('H(next | gap)', 'H(next \\| gap)') not in study:
                fail.append('table %s: line not in STUDY.md: %s' % (name, line[:90]))

# 2. prose numbers: (quoted text, value, format)
v1 = J('verdicts_r1.json'); v3 = J('verdicts_r3.json'); th3 = J('theory_r3.json'); g = J('verdicts_theory_geo.json')
tw = J('verdicts_twins.json'); r6 = J('verdicts_r6.json'); r7 = J('r7.json'); r8 = J('r8_tail.json'); b8 = J('r8_blocks.json')
rows = {json.loads(l)['window']: json.loads(l) for l in open(os.path.join(D, 'results_r3.jsonl')) if json.loads(l)['block'] == -1}
lan = {json.loads(l)['window']: json.loads(l)['stats'] for l in open(os.path.join(D, 'results_landau.jsonl')) if json.loads(l)['block'] == -1}
twr = {json.loads(l)['window']: json.loads(l)['stats'] for l in open(os.path.join(D, 'results_twins.jsonl')) if json.loads(l)['block'] == -1}
PROSE = [
    ('5,096,876 primes', rows['train']['members'], '{:,} primes'),
    ('2,171,134', rows['A']['members'], '{:,}'),
    ('1,808,833', rows['B']['members'], '{:,}'),
    ('578,837', rows['C']['members'], '{:,}'),
    ('0.70556 / 0.77451 / 0.81207', [th3[w]['-1']['S.var50'] for w in ('train', 'A', 'B')], '{:.5f} / {:.5f} / {:.5f}'),
    ('0.70526 / 0.77551 / 0.81187', [v3[w]['S.var50']['x'] for w in ('train', 'A', 'B')], '{:.5f} / {:.5f} / {:.5f}'),
    ('1.0010, 1.0005, 1.0003 and 1.0002', [th3[w]['-1']['norm_gap_law'] for w in ('train', 'A', 'B', 'C')],
     '{:.4f}, {:.4f}, {:.4f} and {:.4f}'),
    ('0.406 at train, 0.694 at A and 0.831 at B', [g[w]['G.chi2_BH']['T'] for w in ('train', 'A', 'B')],
     '{:.3f} at train, {:.3f} at A and {:.3f} at B'),
    ('0.386, 0.715 and 0.721', [g[w]['G.chi2_BH']['x'] for w in ('train', 'A', 'B')], '{:.3f}, {:.3f} and {:.3f}'),
    ('381,332 / 124,843 / 86,479 / 22,013', [twr[w]['n_twins'] for w in ('train', 'A', 'B', 'C')], '{:,} / {:,} / {:,} / {:,}'),
    ('0.72637 against 0.72634', [tw['train']['W2.le300']['x'], tw['train']['W2.le300']['T']], '{:.5f} against {:.5f}'),
    ('order 1 (independent interior pairs) gives 0.73687', tw['train']['W2.le300']['order1'], 'order 1 (independent interior pairs) gives {:.5f}'),
    ('0.3046 among twin gaps', tw['train']['W3.same5']['x'], '{:.4f} among twin gaps'),
    ('0.3043 predicted', tw['train']['W3.same5']['T'], '{:.4f} predicted'),
    ('94,545, 37,081 and 5,974', [lan[w]['count'] for w in ('F1', 'F2', 'F3')], '{:,}, {:,} and {:,}'),
    ('0.0683 and 0.0285 at train', [r6['train']['N.d2_+2']['T'], r6['train']['N.d2_+0']['x']], '{:.4f} and {:.4f} at train'),
    ('0.068348 and 0.068196 at train', [r6['train']['N.d2_-2']['x'], r6['train']['N.d2_+2']['x']], '{:.6f} and {:.6f} at train'),
    ('0.1320 against 0.1321', [r7['train_Q210']['H1_pred'] - r7['train_Q210']['H2_pred'],
                               r7['train_Q210']['H1_meas'] - r7['train_Q210']['H2_meas']], '{:.4f} against {:.4f}'),
    ('0.0971 against 0.0971 at 10¹⁰', [r7['A_Q210']['H1_pred'] - r7['A_Q210']['H2_pred'],
                                      r7['A_Q210']['H1_meas'] - r7['A_Q210']['H2_meas']], '{:.4f} against {:.4f} at 10¹⁰'),
    ('0.0746 against 0.0744 at 10¹²', [r7['B_Q210']['H1_pred'] - r7['B_Q210']['H2_pred'],
                                      r7['B_Q210']['H1_meas'] - r7['B_Q210']['H2_meas']], '{:.4f} against {:.4f} at 10¹²'),
    ('0.0511 against 0.0506 at 10¹⁵', [r7['C_Q210']['H1_pred'] - r7['C_Q210']['H2_pred'],
                                      r7['C_Q210']['H1_meas'] - r7['C_Q210']['H2_meas']], '{:.4f} against {:.4f} at 10¹⁵'),
    ('gaps ≥ 54: 0.035827 → 0.035722; primes 0.035494', [r8['round3']['tail54'], r8['round8']['tail54'], r8['primes']['tail54']],
     'gaps ≥ 54: {:.6f} → {:.6f}; primes {:.6f}'),
    ('from −4.09 to −2.80 (≥ 54), −2.30 to −1.19 (≥ 72) and −1.41 to −0.37 (≥ 90)',
     [b8['round3']['T1.tail54']['t9'], b8['round8']['T1.tail54']['t9'], b8['round3']['T1.tail72']['t9'],
      b8['round8']['T1.tail72']['t9'], b8['round3']['T1.tail90']['t9'], b8['round8']['T1.tail90']['t9']],
     'from {:.2f} to {:.2f} (≥ 54), {:.2f} to {:.2f} (≥ 72) and {:.2f} to {:.2f} (≥ 90)'),
]
# 2b. prose numbers of rounds 45-62 (added after the seventeenth review)
r45b = J('r45b.json'); r50 = J('r50.json'); r52 = J('r52.json'); r53g = J('r53g.json')['compare']; r51f = J('r51f.json')['summary']
r54 = J('r54.json')['bands']; r55 = J('r55.json')['bands']; r56 = J('r56.json'); r58t = J('r58_theory.json'); r59 = J('r59.json')['compare']
r59b = J('r59b.json'); r60 = J('r60.json'); r61 = J('r61.json'); r48 = J('r48.json')['bands']
PROSE += [
    ('= 0.7627474', r52['K_nats'], '= {:.7f}'),
    ('= 0.26965 bits', r56['I_inf_bits'], '= {:.5f} bits'),
    ('= **0.26965 bits**', r56['I_inf_bits'], '= **{:.5f} bits**'),
    ('= **0.15093 bits**', r58t['J_inf_bits'], '= **{:.5f} bits**'),
    ('0.11434 nats = 0.16496 bits', [r60['K_tw_nats'], r60['K_tw_bits']], '{:.5f} nats = {:.5f} bits'),
    ('χ² = 30.7 on 30 bands', [r59['chi2'], r59['dof']], 'χ² = {:.1f} on {} bands'),
    ('gives **2,294**', round(r59b['chi2_pair_factor_baseline']), 'gives **{:,}**'),
    ('(largest 0.063 % at g = 640)', 100 * r53g['640']['rel_o3'], '(largest {:.3f} % at g = 640)'),
    ('Mean U is **1.160**', r51f['order3_overlap']['mean'], 'Mean U is **{:.3f}**'),
    ('72,370 | 73,950 (z −5.8) | 72,432', [r48['120-178_o3']['count'], r48['120-178_o3']['law']], '{:,} | 73,950 (z −5.8) | {:,.0f}'),
    ('| **pooled** | 9.7–18.2 | 4,446 | 4,547 | 0.978 |', [r54['pooled']['count'], r54['pooled']['law'], r54['pooled']['ratio']], '| **pooled** | 9.7–18.2 | {:,} | {:,.0f} | {:.3f} |'),
    ('| **pooled** | 10.0–16.3 | 195 | 182.0 | 1.07 |', [r55['pooled']['count'], r55['pooled']['law'], r55['pooled']['ratio']], '| **pooled** | 10.0–16.3 | {} | {:.1f} | {:.2f} |'),
    ('agree to within 0.0016', max(abs(x['diff']) for x in r50['hit']), 'agree to within {:.4f}'),
    ('| 100 | 1,000 | 10,000 |', None, None),
]
for quoted, val, fmt in PROSE:
    if fmt is None:
        if quoted not in study:
            fail.append('prose: quoted text not found in STUDY.md: %s' % quoted)
        continue
    got = fmt.format(*val) if isinstance(val, list) else fmt.format(val)
    got = got.replace('-', '−')
    if quoted not in study:
        fail.append('prose: quoted text not found in STUDY.md: %s' % quoted)
    elif got != quoted:
        fail.append('prose: STUDY.md says "%s", data give "%s"' % (quoted, got))

# 3. identities
import theory2 as T2
if abs(T2.S1([0, 2]) - 2 * 0.66016181584686957393) > 1e-6:
    fail.append('S(0,2) != 2 C2')
for B in ([0, 2], [0, 6], [0, 2, 6], [0, 4, 6, 10]):
    avg = sum(T2.S_cond(a, 30)(np.array([B]))[0] for a in range(30)) / 30
    if abs(avg - T2.S1(B)) > 1e-9 * T2.S1(B):
        fail.append('S_a does not average to S for %s' % B)
import theory_landau as TL
if abs(TL.Sf_many(np.array([[0]]))[0] - 1.3728134628) > 5e-5:
    fail.append('S_f({0}) far from the Bateman-Horn constant of n^2+1')

print('\n'.join(fail) if fail else 'all checks pass: %d table generators, %d prose numbers, 3 identities' % (18, len(PROSE)))
sys.exit(1 if fail else 0)
