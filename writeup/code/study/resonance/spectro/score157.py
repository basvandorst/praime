# Round 157: score the blind run against the registered decisions (out_blind157.json from spectrometer_blind.py).
import json, sys, statistics, cypari2
pari = cypari2.Pari()
d = json.load(open(sys.argv[1])); ts = json.load(open(sys.argv[2]))
assert set(d) == set(ts) and len(d) == 40
rows = []
for name in sorted(d, key=lambda s: d[s]['N_true']):
    v = d[name]; E = pari.ellinit(ts[name]); N = int(pari.ellglobalred(E)[0]); assert N == v['N_true']
    fa = pari.factor(N); semist = all(int(e) == 1 for e in fa[1])
    hit = v['N_hat'] == v['N_true'] and v['eps_hat'] == v['eps_true']
    top3 = v['rank_true'] is not None and v['rank_true'] <= 2
    # sign at the true conductor: z_true is z(N, eps_true); z(N, -eps) = -z(N, eps)
    sign_ok = v['z_true'] is not None and v['z_true'] > 0
    rows.append((name, N, v['eps_true'], v['c_true'], semist, hit, top3, sign_ok, v['z_true'], v['rank_true'], v['N_hat'], v['eps_hat'], v['z_hat'], v['c_hat'], v['runners_up'][0]))
    print('%-10s N=%4d eps=%+d c=%.3f %-10s %s  z_true=%5.1f rank %s | best N=%4d eps=%+d z=%5.1f c_hat=%.3f | next %d%+d z=%.1f' % (
        name, N, v['eps_true'], v['c_true'], 'semistable' if semist else 'additive', 'HIT ' if hit else 'MISS', v['z_true'], v['rank_true'],
        v['N_hat'], v['eps_hat'], v['z_hat'], v['c_hat'], v['runners_up'][0]['N'], v['runners_up'][0]['eps'], v['runners_up'][0]['z']))
H = sum(r[5] for r in rows); T3 = sum(r[6] for r in rows); SG = sum(r[7] for r in rows)
dev = [abs(r[13] / r[3] - 1) for r in rows if r[5]]
print('D1 hits %d of 40 (pass >= 30, falsifier < 20): %s' % (H, 'PASS' if H >= 30 else 'FAIL'))
print('D2 truth among the three largest: %d of 40 (pass >= 34): %s' % (T3, 'PASS' if T3 >= 34 else 'FAIL'))
print('D3 root number at the true conductor: %d of 40 (pass >= 36): %s' % (SG, 'PASS' if SG >= 36 else 'FAIL'))
print('D4 median |c_hat/c - 1| over the hits: %.4f (pass <= 0.05): %s ; largest %.3f; within 0.10: %d of %d' % (
    statistics.median(dev), 'PASS' if statistics.median(dev) <= 0.05 else 'FAIL', max(dev), sum(x <= 0.10 for x in dev), len(dev)))
print('by stratum: N <= 500: %d of %d ; N > 500: %d of %d' % (sum(r[5] for r in rows if r[1] <= 500), sum(1 for r in rows if r[1] <= 500),
      sum(r[5] for r in rows if r[1] > 500), sum(1 for r in rows if r[1] > 500)))
print('by reduction: semistable %d of %d ; additive %d of %d' % (sum(r[5] for r in rows if r[4]), sum(1 for r in rows if r[4]),
      sum(r[5] for r in rows if not r[4]), sum(1 for r in rows if not r[4])))
print('z of the truth: smallest %.1f, median %.1f ; margin z_true - z_next for the hits: smallest %.1f' % (
    min(r[8] for r in rows), statistics.median([r[8] for r in rows]), min(r[12] - r[14]['z'] for r in rows if r[5])))
