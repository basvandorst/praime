# finite-X model: local factors at l <= L0 from the observed types of the same fields; marginal (independent) vs joint (l <= 7)
import sys, os, time, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pred4 as P, cypari2
from collections import Counter
H = os.path.dirname(os.path.abspath(__file__))
pari = cypari2.Pari(); pari.allocatemem(200_000_000); pari('ltype(nf, l) = vecsort(apply(q -> [q.e, q.f], idealprimedec(nf, l)))')
import qcounts as QC
ue = np.array([0.0] + [(2 * j + 1) / 24 for j in range(24)]); ye = ue**3
L1 = json.load(open(os.path.join(H, 'look1.json')))
base = P.quartic_model(T=1000.0); w = base.c + 1j * base.t
LS = [int(l) for l in P.PR[P.PR <= 47]]
for r2 in (1, 2, 0):
    z = np.load(os.path.join(H, 'data', 'disc_r%d.npz' % r2)); n = len(z['D'])
    rows = [l.split() for l in open(os.path.join(H, 'data', 'disc_r%d_fields.txt' % r2)) if not l.startswith('#')]
    rows = sorted(set((int(D), f) for D, f in rows), key=lambda r: (abs(r[0]), r[1]))
    T = []
    for D, f in rows:
        nf = pari.nfinit(pari(f)); T.append([QC.type_name([(int(e), int(ff)) for e, ff in pari('ltype')(nf, l)]) for l in LS])
    T = np.array(T); assert (T[:, :4] == z['types']).all()
    m = np.array(L1[str(r2)]['m']); s = np.array(L1[str(r2)]['s']); mus = P.SIG_MU[r2]
    def report(name, mod):
        pr = P.bin_means(mod, mus, ye); out = []
        for lo, hi in ((4, 13), (13, 24), (4, 24)):
            sl = slice(lo, hi); ww = 1 / s[sl]**2; a = np.sum(pr[sl] * m[sl] * ww) / np.sum(pr[sl]**2 * ww); ae = 1 / np.sqrt(np.sum(pr[sl]**2 * ww))
            out.append('%.3f+-%.3f chi2 %.1f/%d' % (a, ae, np.sum((m[sl] - pr[sl])**2 * ww), hi - lo))
        print('  r2=%d %-14s A(1)x%.3f y^-1/3 %.2f | %s' % (r2, name, mod.D[0].real / base.D[0].real, P.tail_coeffs(mod, mus)['0.3333'], ' | '.join(out)), flush=True)
    for L0 in (7, 13, 23, 47):
        ov = {}
        for k, l in enumerate(LS):
            if l > L0: break
            fr = {nm: c / n for nm, c in Counter(T[:, k]).items()}; ov[l] = (lambda X, fr=fr, l=l: P.E_local(l, X, fr))
        report('marginal l<=%d' % L0, P.quartic_model(T=1000.0, local_override=ov))
    # joint over l <= 7
    k4 = [LS.index(l) for l in (2, 3, 5, 7)]; tup = Counter(tuple(r[k4]) for r in T)
    Xs = {l: np.exp(-w * np.log(l)) for l in (2, 3, 5, 7)}
    num = sum(c / n * np.prod([P.F_type(nm, l, Xs[l]) for nm, l in zip(tp, (2, 3, 5, 7))], axis=0) for tp, c in tup.items())
    den = np.prod([P.E_local(l, Xs[l]) for l in (2, 3, 5, 7)], axis=0)
    j3 = sum(c / n * np.prod([P.F_type(nm, l, l**(-1/3)) for nm, l in zip(tp, (2, 3, 5, 7))]) for tp, c in tup.items()) / np.prod([P.E_local(l, l**(-1/3)) for l in (2, 3, 5, 7)])
    mj = P.Model(base.c, base.t, base.D * num / den, [base.poles[0], (1 / 3, base.poles[1][1] * j3)])
    report('joint l<=7', mj)
