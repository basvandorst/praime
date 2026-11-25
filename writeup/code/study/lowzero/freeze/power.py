"""power.py: expected errors of the registered readings for a test of n characters, from the development characters
only (dryrun_dev.json, made by freeze_test.py on the 32 development characters).
Method: draw n characters with replacement from the 32 development records (numpy default_rng(4242)), recompute every
reading (freeze_test.readings with the weights of the development run), 400 draws; the sd over draws is the expected
se at size n (the development design: 12, 11, 9, 10 of 32 with p = 2, 3, 5, 7 | D; the test design has 96, 81, 64,
56 of 240, the same proportions within 0.1).  Also: the n at which H_W and H_G (gamma 0 and 1) are 5 se apart, and
at which a halfway truth (gamma = 0.5) is 5 se from both.
usage: power.py dryrun_dev.json predictions_dev.json"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import freeze_test as ft


def main():
    dev = json.load(open(sys.argv[1])); pred = json.load(open(sys.argv[2]))
    recs = dev['chars']
    full, ctl = ft.core(recs, np.arange(len(recs)))
    rng0 = np.random.default_rng(1)
    boots = [ft.core(recs, rng0.integers(0, len(recs), len(recs))) for _ in range(400)]
    W = {}
    for s in ft.STATS:
        cb = np.array([b[0][s][0] for b in boots]); rb = np.array([b[0][s][1] for b in boots])
        W['c_' + s] = 1 / np.var(cb, axis=0); W['d_' + s] = 1 / np.var(cb - rb, axis=0)
    keys = ['G_M', 'G_lnZ_3', 'G_varM', 'G_lnZ_1.5', 'G_lnZ_2', 'A_crit_M', 'A_crit_lnZ_3', 'A_crit_slope32',
            'C1_A_EZ1', 'c_M_2', 'c_M_7', 'c_lnZ_3_2', 'c_varM_2']
    rng = np.random.default_rng(4242)
    res = {}
    for n in (32, 120, 240, 400):
        vals = {k: [] for k in keys}
        for _ in range(400):
            idx = rng.integers(0, len(recs), n)
            co, ct = ft.core(recs, idx)
            rd = ft.readings(co, ct, pred, W)
            for k in keys:
                vals[k].append(rd[k])
        res[n] = {k: float(np.std(vals[k])) for k in keys}
        print('n = %3d  ' % n + '  '.join('%s %.4f' % (k, res[n][k]) for k in keys), flush=True)
    se240 = res[240]
    print('\nat n = 240:')
    print('  gamma (relational, H_W 0 vs H_G 1): M %.3f, lnZ_3 %.3f, varM %.3f -> H_W/H_G apart by %.1f, %.1f, %.1f se'
          % (se240['G_M'], se240['G_lnZ_3'], se240['G_varM'], 1 / se240['G_M'], 1 / se240['G_lnZ_3'],
             1 / se240['G_varM']))
    for k, gap in (('A_crit_M', 0.165), ('A_crit_lnZ_3', 0.164)):
        print('  %s se %.4f: H_lead - H_crit = %.3f = %.1f se (inside the theory band only %.3f = %.1f se)' % (
            k, se240[k], gap, gap / se240[k], gap - 0.08, (gap - 0.08) / se240[k]))
    print('  A_crit_slope32 se %.4f: naive (3.8) and spread (3.3) are %.0f and %.0f se from frozen (1)' % (
        se240['A_crit_slope32'], 2.8 / se240['A_crit_slope32'], 2.3 / se240['A_crit_slope32']))
    for k in ('G_M', 'G_lnZ_3'):
        s32 = res[32][k]
        print('  %s: n for 5 se between 0 and 1: %d; n for a halfway truth 5 se from both: %d' % (
            k, int(np.ceil(32 * (s32 / 0.2) ** 2)), int(np.ceil(32 * (s32 / 0.1) ** 2))))
    json.dump({str(n): res[n] for n in res}, open('power_out.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
