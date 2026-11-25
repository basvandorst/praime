# Reviewer's regressions: one-point (O1, O2, per line, per band), c_E dependence, X1/X5/X6 with my own
# delete-one-curve jackknife and a dyadic (node-clustered) sandwich, comparator from measured lines,
# height splits, influence, subgroups.
import os, numpy as np
SC = os.path.dirname(os.path.abspath(__file__))+'/'
Z = np.load(SC + 'rv_zeros.npz'); LW = np.load(SC + 'rv_law.npz'); CU = np.load(SC + 'rv_curves.npz')
PF = np.load(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/ecross/pred_fresh.npz')
ids = Z['ids']; nc = len(ids)
assert np.array_equal(ids, LW['ids']) and np.array_equal(ids, CU['ids'])
idx = {s: k for k, s in enumerate(ids)}
pf2my = np.array([idx[s] for s in PF['ids']])
# frozen pairs -> my curve indices
PI, PJ = pf2my[PF['I']], pf2my[PF['J']]
w, N, rank, c, amp = CU['w'], CU['N'], CU['rank'], LW['c'], LW['amp']
RAT = [tuple(x) for x in Z['rat']]
LINES13 = [tuple(x) for x in LW['lines']]
ri13 = [RAT.index(x) for x in LINES13]


def Xmat(key):
    M = np.zeros((nc, nc), complex); X = Z[key]; M[Z['I'], Z['J']] = X; M[Z['J'], Z['I']] = np.conj(X)
    return M


def fit(Y, Xd, units, nunits=nc):
    """OLS; delete-one-unit jackknife (units: list of index arrays per row, 1 or 2 units) and clustered sandwich."""
    Xd = np.atleast_2d(Xd.T).T if Xd.ndim == 1 else Xd
    A = Xd.T @ Xd; b = Xd.T @ Y
    th = np.linalg.solve(A, b)
    e = Y - Xd @ th
    S = Xd * e[:, None]
    Ak = np.zeros((nunits,) + A.shape); bk = np.zeros((nunits, len(b))); Sk = np.zeros((nunits, len(b)))
    for col in range(units.shape[1]):
        u = units[:, col]
        np.add.at(Ak, u, Xd[:, :, None] * Xd[:, None, :])
        np.add.at(bk, u, Xd * Y[:, None])
        np.add.at(Sk, u, S)
    thk = np.array([np.linalg.solve(A - Ak[k], b - bk[k]) for k in range(nunits)])
    jk = np.sqrt((nunits - 1) / nunits * ((thk - thk.mean(0)) ** 2).sum(0))
    meat = Sk.T @ Sk
    if units.shape[1] == 2:
        meat -= S.T @ S
    Ai = np.linalg.inv(A)
    sw = np.sqrt(np.diag(Ai @ meat @ Ai))
    return th, jk, sw, thk


def show(name, th, jk, sw):
    print('  %-62s %s' % (name, '  '.join('%9.4f +- %.4f (sw %.4f)' % (t, j, s) for t, j, s in zip(th, jk, sw))))


T1 = 380.0
units1 = np.arange(nc)[:, None]
print('=== one-point (my theta, my c_E, my rho) ===')
for bn, Tb in (('b20_400', 380.), ('b20_200', 180.), ('b200_400', 200.)):
    G = Z['G_' + bn]
    law = amp * Tb / (2 * np.pi)
    th, jk, sw, _ = fit(G[:, ri13[0]].real, np.c_[np.ones(nc), law[:, 0]], units1)
    show('O1 r=1 %s [intercept, slope]' % bn, th, jk, sw)
    # O2 pooled 12 lines, no intercept
    Y = G[:, ri13[1:]].real.ravel(); Xd = law[:, 1:].ravel()[:, None]
    un = np.repeat(np.arange(nc), 12)[:, None]
    th, jk, sw, _ = fit(Y, Xd, un)
    show('O2 lines r!=1 pooled %s' % bn, th, jk, sw)
    Y = G[:, ri13].real.ravel(); Xd = law.ravel()[:, None]
    un = np.repeat(np.arange(nc), 13)[:, None]
    th, jk, sw, _ = fit(Y, Xd, un)
    show('O12 all 13 lines pooled %s' % bn, th, jk, sw)
    print('  Im G mean (13 lines) %s: %.3f, rms Re(G-law) r=1: %.2f' % (bn, G[:, ri13].imag.mean(), np.std(G[:, ri13[0]].real - law[:, 0])))
# per line
G = Z['G_b20_400']; law = amp * T1 / (2 * np.pi)
print('  per line slope (no intercept):', ' '.join('%d/%d:%.3f' % (a, b, np.sum(G[:, ri13[k]].real * law[:, k]) / np.sum(law[:, k] ** 2)) for k, (a, b) in enumerate(LINES13)))

print('=== c_E dependence of the r = 1 line (band 20-400) ===')
y = G[:, ri13[0]].real
base = -w * T1 / (2 * np.pi)
LS1 = PF['LS1'][np.argsort(pf2my)]; B = PF['B'][np.argsort(pf2my)]
assert np.allclose(PF['c'][np.argsort(pf2my)], c, rtol=2e-3)
for nm, f in (('c_E = B/L(1,Sym2)', c), ('1 (constant)', np.ones(nc)), ('1/L(1,Sym2)', 1 / LS1), ('B alone', B), ('c_E^2', c ** 2), ('sqrt c_E', np.sqrt(c))):
    x = base * f
    th, jk, sw, _ = fit(y, x[:, None], units1)
    rss = np.sum((y - th[0] * x) ** 2)
    print('  amplitude ~ %-20s slope %.4f +- %.4f  RSS %.0f' % (nm, th[0], jk[0], rss))
th, jk, sw, _ = fit(y, np.c_[base * c, base * (1 / LS1 - c)], units1)
show('[c_E, 1/L(1,Sym2) - c_E]', th, jk, sw)
th, jk, sw, _ = fit(y, np.c_[base * c, base * (np.ones(nc) - c)], units1)
show('[c_E, 1 - c_E]', th, jk, sw)
gam = np.linspace(0, 2, 201); rss = []
for g in gam:
    x = base * c ** g; s = np.sum(x * y) / np.sum(x * x); rss.append(np.sum((y - s * x) ** 2))
rss = np.array(rss); s2 = rss.min() / (nc - 2)
ok = gam[(rss - rss.min()) / s2 <= 1.0]
print('  power law amplitude ~ c^gamma: best gamma %.2f, 1-sigma range %.2f..%.2f' % (gam[rss.argmin()], ok.min(), ok.max()))

print('=== pairs (my X, frozen pred) ===')
units2 = np.c_[PI, PJ]
for sig, key in ((0.5, 'X_0.5'), (1.0, 'X_1.0')):
    M = Xmat(key)
    ReX = M[PI, PJ].real; ImX = M[PI, PJ].imag
    pred = PF['pred_%s' % sig]; pred1 = PF['pred1_%s' % sig]
    th, jk, sw, thk = fit(ReX, np.c_[np.ones(len(ReX)), pred], units2)
    show('s%s X1 [a, slope]' % sig, th, jk, sw)
    if sig == 0.5:
        thk_X1 = thk; th_X1 = th
    Tw = 376.0
    Gw = Z['G_w20_396']
    kern = lambda a, b: sig * np.sqrt(2 * np.pi) * np.exp(-sig ** 2 * np.log(a / b) ** 2 / 2)
    P0 = kern(1, 1) * (Gw[PI, ri13[0]] * np.conj(Gw[PJ, ri13[0]])).real / Tw ** 2
    th, jk, sw, _ = fit(ReX - P0, np.c_[np.ones(len(ReX)), pred - pred1], units2)
    show('s%s X5 [a, H-1 coefficient]' % sig, th, jk, sw)
    # X6 with the 13 lines
    meas13 = sum(kern(a, b) * (Gw[PI, ri13[k]] * np.conj(Gw[PJ, ri13[k]])).real / Tw ** 2 for k, (a, b) in enumerate(LINES13))
    law13 = sum(kern(a, b) * amp[PI, k] * amp[PJ, k] / (4 * np.pi ** 2) for k, (a, b) in enumerate(LINES13))
    X6 = ReX - meas13 - (pred - law13)
    th, jk, sw, _ = fit(X6, np.ones((len(X6), 1)), units2)
    show('s%s X6 continuum mean (13 lines)' % sig, th, jk, sw)
    print('    corr(pred, law13) = %.3f ; mean(pred - law13) = %.4f' % (np.corrcoef(pred, law13)[0, 1], np.mean(pred - law13)))
    th, jk, sw, _ = fit(ReX, np.c_[np.ones(len(ReX)), meas13], units2)
    show('s%s Re X on measured 13-line products (no law)' % sig, th, jk, sw)
    th, jk, sw, _ = fit(ImX, np.ones((len(X6), 1)), units2)
    show('s%s Im X mean' % sig, th, jk, sw)
    if sig == 0.5:
        for tag in ('_lo', '_hi'):
            Mh = Xmat('X_0.5' + tag)
            th, jk, sw, _ = fit(Mh[PI, PJ].real, np.c_[np.ones(len(ReX)), pred], units2)
            show('s0.5 X1 window %s' % ('20-208' if tag == '_lo' else '208-396'), th, jk, sw)
        # influence
        print('  X1 slope without the k curves of largest jackknife influence:')
        infl = thk_X1[:, 1] - th_X1[1]
        order = np.argsort(-np.abs(infl))
        print('   top-5 curve influences on slope:', ' '.join('%s(c=%.2f):%+.4f' % (ids[k], c[k], infl[k]) for k in order[:5]))
        for kdrop in (5, 10):
            keep = ~(np.isin(PI, order[:kdrop]) | np.isin(PJ, order[:kdrop]))
            th, jk, sw, _ = fit(ReX[keep], np.c_[np.ones(keep.sum()), pred[keep]], units2[keep])
            show('   drop %d most influential curves (%d pairs)' % (kdrop, keep.sum()), th, jk, sw)
        top = np.argsort(-np.abs(pred))
        for kdrop in (10, 100):
            keep = np.ones(len(pred), bool); keep[top[:kdrop]] = False
            th, jk, sw, _ = fit(ReX[keep], np.c_[np.ones(keep.sum()), pred[keep]], units2[keep])
            show('   drop %d pairs of largest |pred|' % kdrop, th, jk, sw)
        # share of the slope's leverage in top pairs
        lev = (pred - pred.mean()) ** 2; lev = np.sort(lev)[::-1]
        print('   leverage share of top 10 / 100 / 1000 pairs: %.3f %.3f %.3f' % (lev[:10].sum() / lev.sum(), lev[:100].sum() / lev.sum(), lev[:1000].sum() / lev.sum()))
        # subgroups
        for nm, sel in (('both rank 0', (rank[PI] == 0) & (rank[PJ] == 0)), ('both rank >= 1', (rank[PI] >= 1) & (rank[PJ] >= 1)),
                        ('rank 2 involved', (rank[PI] == 2) | (rank[PJ] == 2)),
                        ('both N < 1500', (N[PI] < 1500) & (N[PJ] < 1500)), ('both N >= 1500', (N[PI] >= 1500) & (N[PJ] >= 1500)),
                        ('both c > 1.5', (c[PI] > 1.5) & (c[PJ] > 1.5)), ('both c < 1', (c[PI] < 1) & (c[PJ] < 1)),
                        ('equal conductor', N[PI] == N[PJ])):
            if sel.sum() < 30:
                print('   %s: %d pairs: ReX %s pred %s' % (nm, sel.sum(), np.round(ReX[sel], 3), np.round(pred[sel], 3)))
                continue
            th, jk, sw, _ = fit(ReX[sel], pred[sel][:, None], units2[sel])
            show('   subgroup %s (%d pairs), slope through 0' % (nm, sel.sum()), th, jk, sw)
        np.savez(SC + 'rv_pairs.npz', ReX=ReX, pred=pred, pred1=pred1, P0=P0, meas13=meas13, law13=law13, PI=PI, PJ=PJ)

print('=== one-point subgroups (r = 1 and 13 lines pooled, band 20-400) ===')
for nm, sel in (('rank 0', rank == 0), ('rank 1', rank == 1), ('rank 2', rank == 2), ('w=+1', w == 1), ('w=-1', w == -1),
                ('N<1500', N < 1500), ('N>=1500', N >= 1500), ('c>1.5', c > 1.5), ('c<1', c < 1)):
    k = np.where(sel)[0]
    Y = G[k][:, ri13].real.ravel(); Xd = law[k].ravel()[:, None]
    un = np.repeat(np.arange(len(k)), 13)[:, None]
    th, jk, sw, _ = fit(Y, Xd, un, nunits=len(k))
    th1, jk1, sw1, _ = fit(G[k, ri13[0]].real, law[k, 0][:, None], np.arange(len(k))[:, None], nunits=len(k))
    print('   %-8s (%3d curves): r=1 slope %.4f +- %.4f ; 13 lines %.4f +- %.4f' % (nm, len(k), th1[0], jk1[0], th[0], jk[0]))
