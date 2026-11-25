"""twins_discovery.py: the constellation statistics of twin_lib on the DISCOVERY range p < 3e7 (seen; not blind).
Window: Hann in v = sqrt(n) over [sqrt(6e5), sqrt(3e7)].  Readings = smooth sum / (predicted constant * MT_unit).
Errors: rms of the 8 control-frequency readings (same weights), and the random-model error sqrt(noise/2)/|MT|.
Output: out_twins_discovery.json, printed table."""
import sys, os, json, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import spectro_lib as S
import twin_lib as T

X0, X1 = 0, 30_000_000
v0, v1 = math.sqrt(6e5), math.sqrt(3e7)
for nm in ('11a1', '15a1'):
    z = np.load(os.path.join(S.HERE, 'aptab_%s_30000000.npz' % nm)); T.APTAB[nm] = (z['p'], z['a'])
t0 = time.time(); tot = None; part = {}
CH = 5_000_000
for L in range(X0, X1, CH):
    U = min(X1, L + CH)
    r = T.chunk_sums(L, U, v0, v1)
    tot = T.add_rec(tot, r); part[U] = json.loads(json.dumps(tot))
    print('chunk', L, U, '%.0fs' % (time.time() - t0), flush=True)

consts = {}
for lab in T.LABELS:
    c = S.label_info(lab)['c']
    consts[lab] = dict(single=c, twinlo=T.c_const(lab, 'twinlo'), twinhi=T.c_const(lab, 'twinhi'), sg=T.c_const(lab, 'sg_lo'),
                       shift2=(T.c_const(lab, 'shift2') if lab in T.LATTICE else None), shift6=(T.c_const(lab, 'shift6') if lab in T.LATTICE else None), factorised_twin=T.TWO_C2 * c)
    print(lab, consts[lab], flush=True)
res = {}
for lab in T.LABELS:
    MT = S.main_term_unit(lab, v0, v1)
    K1, kappa = S.K1_kappa(lab)
    res[lab] = {}
    for st, rec in tot[lab].items():
        pred = consts[lab][st]
        sm = [complex(a, b) for a, b in rec['smooth']]
        rd = sm[0] / MT / pred
        ctrl = [z / MT / pred for z in sm[1:]]
        rms = math.sqrt(np.mean([abs(z) ** 2 for z in ctrl]) / 2)
        rnd = math.sqrt(rec['noise'] / 2) / abs(MT) / pred
        # sharp running readings: S_sharp(x) / (pred K1 x^{3/4})
        runs = {}
        for U, pr in part.items():
            sh = complex(*pr[lab][st]['sharp']); runs[U] = sh / (pred * K1 * U ** 0.75)
        res[lab][st] = dict(pred_const=pred, reading=[rd.real, rd.imag], err_ctrl=rms, err_random=rnd, count=rec['count'],
                            ctrl=[[z.real, z.imag] for z in ctrl], sharp_running={str(k): [v.real, v.imag] for k, v in runs.items()})
        print('%-8s %-7s const %.4f  reading %+.3f %+.3fi   err(ctrl rms) %.3f  err(random) %.3f   n=%d   sharp@3e7 %+.3f %+.3fi'
              % (lab, st, pred, rd.real, rd.imag, rms, rnd, rec['count'], runs[X1].real, runs[X1].imag))
json.dump(dict(window=[v0, v1], consts=consts, res=res), open(os.path.join(S.HERE, 'out_twins_discovery.json'), 'w'), indent=1)
print('done %.0fs' % (time.time() - t0))
