#!/usr/bin/env python3
"""evaluate.py REPORT.json -- applies the registered readings of round_draft.txt to a map.py report of the fresh windows.
Self-contained (numpy + standard library).  Thresholds are fixed here and in round_draft.txt; nothing is fitted."""
import sys, json, math
import numpy as np

REG_WINDOWS = [int(3.0e10 + i * 1e8) for i in range(24)]
EXPECT = {   # exploratory values (seen data, combined over the two old heights), for display only
    ('5_1', 'rc'): -0.0089, ('8_1', 'rr'): -0.0071, ('9_1', 'conj'): -0.0090, ('11_1', 'eq'): -0.0036,
    ('11_1', 'rc'): -0.0046, ('11_1', 'cc'): -0.0025, ('13_1', 'cc'): -0.0019}
SECONDARY = [('8_1', 'rr'), ('9_1', 'conj'), ('11_1', 'eq'), ('11_1', 'rc'), ('11_1', 'cc'), ('13_1', 'cc')]
P1M_CENTRE = -0.0087


def jk(v):
    v = np.asarray(v, float); n = len(v); loo = (v.sum() - v) / (n - 1)
    return math.sqrt((n - 1) / n * ((loo - loo.mean()) ** 2).sum())


def main(path):
    R = json.load(open(path)); C = R['cells']
    ok_w = sorted(R['windows']) == REG_WINDOWS
    print('report %s: %d windows, lam %.4f; registered windows: %s' % (path, len(R['windows']), R['lam'], 'yes' if ok_w else 'NO'))
    if not ok_w: print('  WARNING: the windows differ from the registration; the readings below are not the registered test')
    a = C['5_1']['types']['rc']
    E, s = a['excess'], a['err_tot']
    print('\nPRIMARY P1  mod 5, k = 1, mixed directions (type rc, 4 directions)')
    print('  measured slope %+.5f +- %.5f, LR slope %+.5f +- %.5f, excess %+.5f +- %.5f, z = %+.2f'
          % (a['slope'][0], a['slope'][1], a['lr'], a['lr_err'], E, s, E / s))
    print('  P1  (sign):      excess <= -2 err            -> %s' % ('PASS' if E <= -2 * s else 'FAIL'))
    print('  P1+ (strong):    excess <= -3 err            -> %s' % ('yes' if E <= -3 * s else 'no'))
    print('  P1m (magnitude): |excess - (%+.4f)| <= 2.5 err -> %s' % (P1M_CENTRE, 'PASS' if abs(E - P1M_CENTRE) <= 2.5 * s else 'FAIL'))
    print('  FALSIFIER F:     |excess| < 1 err            -> %s' % ('FALSIFIED' if abs(E) < s else 'not triggered'))
    # contrast
    T = C['5_1']['types']; others = [t for t in T if t != 'rc']; nd = {t: T[t]['ndir'] for t in T}
    s_rc = np.array(T['rc']['slope_pw']); s_ot = sum(nd[t] * np.array(T[t]['slope_pw']) for t in others) / sum(nd[t] for t in others)
    l_rc = T['rc']['lr']; l_ot = sum(nd[t] * T[t]['lr'] for t in others) / sum(nd[t] for t in others)
    le = math.hypot(T['rc']['lr_err'], math.sqrt(sum((nd[t] * T[t]['lr_err']) ** 2 for t in others)) / sum(nd[t] for t in others))
    c = (s_rc - l_rc) - (s_ot - l_ot); ce = math.hypot(jk(c), le)
    print('\nSECONDARY S0  mod 5, k = 1: contrast = excess(rc) - excess(the other 5 directions) = %+.5f +- %.5f (z %+.2f) -> %s'
          % (c.mean(), ce, c.mean() / ce, 'PASS' if c.mean() <= -2 * ce else 'FAIL'))
    print('\nSECONDARY S1-S6 (sign kept at >= 2 errors: excess <= -2 err)')
    npass = 0
    for cell, t in SECONDARY:
        b = C[cell]['types'][t]; p = b['excess'] <= -2 * b['err_tot']; npass += p
        print('  q,k = %-5s %-5s excess %+.5f +- %.5f (z %+.2f; exploratory %+.4f) -> %s' % (cell, t, b['excess'], b['err_tot'], b['z'], EXPECT[(cell, t)], 'PASS' if p else 'FAIL'))
    print('  %d of 6 pass (expected about 4.8 if the exploratory values hold, 1.9 if half as large, 0.14 if noise)' % npass)
    print('\nREPORTED, NO READING: every other (q, k, type) cell, and the lag-50 controls')
    for cell in sorted(C, key=lambda x: tuple(map(int, x.split('_')))):
        for t, b in C[cell]['types'].items():
            if (cell, t) in SECONDARY or (cell, t) == ('5_1', 'rc') or 'excess' not in b: continue
            print('  %-5s %-5s %3d dirs  excess %+.5f +- %.5f (z %+.2f)' % (cell, t, b['ndir'], b['excess'], b['err_tot'], b['z']))
        if 'excess' in C[cell]['all']:
            b = C[cell]['all']; print('  %-5s all   %3d dirs  excess %+.5f +- %.5f (z %+.2f)  slope/LR %.2f' % (cell, b['ndir'], b['excess'], b['err_tot'], b['z'], b['slope'][0] / b['lr']))


if __name__ == '__main__':
    main(sys.argv[1])
