"""(a) direct double-sum X (erf window) on 110 OO functions (5,995 pairs) and on 2,100 complex-design pairs vs the t-grid X;
(b) own readings R1-R7, C1-C7 with own regression / jackknife, compared with readings.json."""
import sys, math, time, json
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import numpy as np
from rlib import *

OUT = __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/'
RJ = json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/data/study/r143/readings.json'))
t0 = time.time()
P = np.load(FZ + 'pred_OO.npz'); Ds = [int(x) for x in P['Ds']]
C = np.load(FZ + 'pred_C.npz'); clabs = [str(x) for x in C['labels']]


def show(name, c, se, ref=None):
    s = ', '.join('%.4f%+.4fi' % (complex(x).real, complex(x).imag) if abs(complex(x).imag) > 1e-12 else '%.4f' % complex(x).real for x in c)
    e = ', '.join('%.4f' % x for x in se)
    r = ''
    if ref is not None:
        r = '   | proposer: ' + ', '.join('%.4f%+.4fi' % (v[0], v[1]) if abs(v[1]) > 1e-12 else '%.4f' % v[0] for v in ref['coef']) + ' +- ' + ', '.join('%.4f' % x for x in ref['se_re'])
    print('  %-44s %s +- %s%s' % (name, s, e, r))


# ---------------- (a) direct double sums
Z = load_group(R120F + ['r134_B.npz', 'r136_Dn.npz'])
rng = np.random.default_rng(143)
sub = np.sort(rng.choice(546, 110, replace=False))
for sig in (0.5, 1.0):
    Xg = np.load(OUT + 'Xown_OO_%.1f.npy' % sig)
    ph = {}
    for i in sub:
        D = Ds[i]; g = Z['D%d' % D]; g = g[(g > 90) & (g < 600)]
        ph[i] = (g, np.exp(2j * theta(g, abs(D), 0 if D > 0 else 1)))
    d = []
    for a_ in range(len(sub)):
        for b_ in range(a_ + 1, len(sub)):
            i, j = sub[a_], sub[b_]
            x = X_direct(ph[i][0], ph[i][1], ph[j][0], ph[j][1], sig, 100.0, 594.0)
            d.append((x, Xg[i, j]))
    d = np.array(d)
    print('OO sigma %.1f: direct double sum vs t-grid on %d pairs: max |diff| %.2e, rms X %.4f (%.0fs)' % (sig, len(d), np.abs(d[:, 0] - d[:, 1]).max(), np.sqrt(np.mean(np.abs(d[:, 0]) ** 2)), time.time() - t0), flush=True)
ZC = load_group(CPLXF + R120F)
for sig in (0.5,):
    Xg = np.load(OUT + 'Xown_C_%.1f.npy' % sig)
    sel = rng.choice(236, 66, replace=False)
    ph = {}
    for i in sel:
        lab = clabs[i]
        if lab.startswith('D'):
            D = int(lab[1:]); g = ZC[lab]; q, a = abs(D), 0 if D > 0 else 1
        else:
            N, m = int(lab.split('_')[0][1:]), int(lab.split('_')[1][1:]); z = ZC['N%d_m%d' % (N, m)]
            g = np.sort(-z[z < 0]) if lab.endswith('_conj') else np.sort(z[z > 0]); q = N; a = int(C['a'][i])
        g = g[(g > 96) & (g < 500)]; ph[i] = (g, np.exp(2j * theta(g, q, a)))
    d = []
    for a_ in range(len(sel)):
        for b_ in range(a_ + 1, len(sel)):
            i, j = sel[a_], sel[b_]
            d.append((X_direct(ph[i][0], ph[i][1], ph[j][0], ph[j][1], sig, 106.0, 494.0), Xg[i, j]))
    d = np.array(d)
    print('C sigma %.1f: direct double sum vs t-grid on %d pairs: max |diff| %.2e, rms X %.4f (%.0fs)' % (sig, len(d), np.abs(d[:, 0] - d[:, 1]).max(), np.sqrt(np.mean(np.abs(d[:, 0]) ** 2)), time.time() - t0), flush=True)

# ---------------- (b) readings OO / HI
for des in ('OO', 'HI'):
    I, J = P['I'], P['J']
    if des == 'HI':
        m = (I < 181) & (J < 181)
        excl = np.zeros(546, bool); excl[[Ds.index(-395), Ds.index(197)]] = True
    else:
        m = np.ones(len(I), bool); excl = np.zeros(546, bool)
    I0, J0 = I[m], J[m]; keep = ~(excl[I0] | excl[J0]); I1, J1 = I0[keep], J0[keep]
    q = np.abs(np.array(Ds)); a = (np.array(Ds) < 0).astype(int)
    print('%s: %d pairs' % (des, len(I1)))
    for sig in (0.5, 1.0):
        X = np.load(OUT + 'Xown_%s_%.1f.npy' % (des, sig))
        x = X[I1, J1]; p = P['pred_%.1f' % sig][m][keep]; p1 = P['pred1_%.1f' % sig][m][keep]; one = np.ones(len(x))
        tag = 's%.1f' % sig; R = RJ[des]
        c, se, _, _ = ols_jack([one, p], x.real, I1, J1); show(tag + ' R1', c, se, R.get(tag + ' R1 slope of Re X on the law'))
        c, se, _, reps = ols_jack([one, p1, p - p1], x.real, I1, J1); show(tag + ' R2', c, se, R.get(tag + ' R2 split [H=1 part, H-1 part]'))
        if des == 'OO' and sig == 0.5:
            cc = np.corrcoef(reps[:, 1], reps[:, 2])[0, 1]
            print('    jackknife correlation of [H=1 coef, H-1 coef]: %.3f; corr(pred1, pred-pred1) over pairs %.3f; corr(1, pred1) %.3f' % (cc, np.corrcoef(p1, p - p1)[0, 1], 0))
        dl = np.abs(np.log(q[J1] / q[I1]))
        for lo, hi in ((0, 0.1), (0.1, 0.3), (0.3, 0.6), (0.6, 1.0), (1.0, 9.0)):
            mm = (dl >= lo) & (dl < hi)
            if mm.sum() > 50:
                c, se, _, _ = ols_jack([one, p], x.real, I1, J1, mm); show(tag + ' R3 %.1f-%.1f' % (lo, hi), c, se, R.get(tag + ' R3 slope, %.1f <= |ln(q2/q1)| < %.1f' % (lo, hi)))
        c, se, _, _ = ols_jack([one, p], x.imag, I1, J1); show(tag + ' R4 Im', c, se, R.get(tag + ' R4 slope of Im X on the law'))
        par = a[I1] + a[J1]
        for nm, mm in (('even-even', par == 0), ('odd-odd', par == 2), ('even-odd', par == 1)):
            c, se, _, _ = ols_jack([one, p], x.real, I1, J1, mm); show(tag + ' R5 ' + nm, c, se, R.get(tag + ' R5 slope, ' + nm))
        c, se, _, _ = ols_jack([p], x.real, I1, J1, p < -0.01); show(tag + ' R6', c, se, R.get(tag + ' R6 pairs with law < -0.01: slope of Re X on the law, no intercept'))
        c, se, _, _ = ols_jack([p], x.real, I1, J1); show(tag + ' R7', c, se, R.get(tag + ' R7 slope of Re X on the law, no intercept'))
    print('  (%.0fs)' % (time.time() - t0), flush=True)

# ---------------- complex design
excl = np.array([l == 'N173_m2_conj' for l in clabs])
unit = {}; uid = np.array([unit.setdefault(l.replace('_conj', ''), len(unit)) for l in clabs])
kk = np.array([kappa(int(qq)) for qq in C['q']])
# own eps: Gauss sums of own Conrey tables (complex) or 1 (real); conj side conj
own_eps = []
for l in clabs:
    if l.startswith('D'): own_eps.append(1 + 0j); continue
    N, m_ = int(l.split('_')[0][1:]), int(l.split('_')[1][1:]); tab = conrey_table(N, m_); a_ = 0 if abs(tab[N - 1] - 1) < 1e-9 else 1
    e = gauss_eps(tab, a_); own_eps.append(np.conj(e) if l.endswith('_conj') else e)
own_eps = np.array(own_eps)
print('C: max |own eps - frozen eps| %.1e; max |own kappa - frozen| %.1e' % (np.abs(own_eps - C['eps']).max(), np.abs(kk - C['kappa']).max()))
for sig in (0.5, 1.0):
    X = np.load(OUT + 'Xown_C_%.1f.npy' % sig); tag = 's%.1f' % sig; R = RJ['C']
    I = np.r_[C['I_cc'], C['I_cr']]; J = np.r_[C['J_cc'], C['J_cr']]
    JJ = np.r_[C['J_cc_%.1f' % sig], C['J_cr_%.1f' % sig]]
    keep = ~(excl[I] | excl[J]); I, J, JJ = I[keep], J[keep], JJ[keep]
    iscc = np.r_[np.ones(len(C['I_cc']), bool), np.zeros(len(C['I_cr']), bool)][keep]
    base = kk[I] * kk[J] * own_eps[I] * np.conj(own_eps[J]) / (4 * math.pi ** 2)
    p = base * JJ                                    # law rebuilt from frozen J with own eps, kappa
    pf = np.r_[C['pred_cc_%.1f' % sig], C['pred_cr_%.1f' % sig]][keep]
    if sig == 0.5: print('  law rebuilt (own eps, kappa) vs frozen pred: max diff %.1e' % np.abs(p - pf).max())
    pn = kk[I] * kk[J] * JJ / (4 * math.pi ** 2); p1 = base * sig * math.sqrt(2 * math.pi)
    x = X[I, J]; one = np.ones(len(x), complex); ui, uj = uid[I], uid[J]
    mirror = ui == uj
    c, se, sei, _ = ols_jack([one, p], x, ui, uj); show(tag + ' C1', c, se, R.get(tag + ' C1 complex slope of X on the law (CC + CR)'))
    c, se, sei, _ = ols_jack([one, p, pn], x, ui, uj); show(tag + ' C2 [law, noroot]', c, se, R.get(tag + ' C2 two predictors [law, no root numbers]'))
    c, se, sei, _ = ols_jack([one, base * JJ.real, base * 1j * JJ.imag], x, ui, uj); show(tag + ' C3', c, se, R.get(tag + ' C3 [Re J part, i Im J part] (direction: law 1, conj J -1)'))
    c, se, sei, _ = ols_jack([one, p], x, ui, uj, mirror); show(tag + ' C4 mirror', c, se, R.get(tag + ' C4 (chi, conj chi) pairs: complex slope'))
    c, se, sei, _ = ols_jack([one, p], x, ui, uj, ~iscc); show(tag + ' C5 CR', c, se, R.get(tag + ' C5 complex x real pairs: complex slope'))
    c, se, sei, _ = ols_jack([one, p], x, ui, uj, iscc & ~mirror); show(tag + ' C6 CC', c, se, R.get(tag + ' C6 complex x complex, chi2 != conj chi1: complex slope'))
    c, se, sei, _ = ols_jack([one, p1, p - p1], x, ui, uj); show(tag + ' C7', c, se, R.get(tag + ' C7 split [H=1 part, H-1 part] (CC + CR)'))
print('done %.0fs' % (time.time() - t0))
