"""Predictions of the alternatives for every sample (no zeros used here except for the noise model):
   BK      the Bogomolny-Keating law (gamma = 1 by construction)
   AMP     the bad-prime off-diagonal change scaled by 0.932 (the "plain amplitude")
   LP3     the bad-prime off-diagonal change evaluated at L + 0.3 (diagonal at L)
For each sample: BK full/diag per L-function (cached), the readings gamma, between, within, pooled that each alternative
would give if it were exactly true (noise-free), and the expected statistical error of gamma from a noise model:
residual rows of the rounds 117/120 sample (100-600) scaled by sqrt(zeros in reference / zeros in the new band).
Samples: B (400<|D|<=700, 100-600), R120 (100<=|D|<=400, 100-600), C (R120's D, 600-1000), Bp (B's D, 600-1000),
Dn (700<|D|<=1000, 100-600).  usage: alt.py SAMPLE  (caches alt_SAMPLE.npz) | alt.py report"""
import os, sys, json, math, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'bk'))
from common import *
import gridcache

def fund_discs(lo, hi):
    out = []
    for D in list(range(-hi, -lo + 1)) + list(range(lo, hi + 1)):
        if D in (0, 1): continue
        if D % 4 == 1 and all(D % (p * p) for p in range(2, int(abs(D) ** .5) + 1)): out.append(D)
        elif D % 4 == 0:
            m = D // 4
            if m % 4 in (2, 3) and all(m % (p * p) for p in range(2, int(abs(m) ** .5) + 2)): out.append(D)
    return sorted(out, key=lambda x: (abs(x), x))

P = json.load(open(os.path.join(LOW, 'bk', 'predictions_per_lfunction.json')))
DB = P['B']['D']; DC = P['C']['D']
SAMPLES = {'B': (DB, 100, 600), 'R120': (DC, 100, 600), 'C': (DC, 600, 1000), 'Bp': (DB, 600, 1000),
           'Dn': (fund_discs(701, 1000), 100, 600)}

def nzeros(D, lo, hi):
    k = abs(D); N = lambda t: t / (2 * math.pi) * math.log(k * t / (2 * math.pi * math.e))
    return N(hi) - N(lo)

def readings(S, F, Dg, X):
    g = gamma(S, F, Dg, X)
    H = X @ np.linalg.pinv(X.T @ X) @ X.T; x = F - Dg; y = S - Dg
    xb = H @ x; xw = x - xb; yb = H @ y; yw = y - yb
    xbc = xb - xb.mean(0); ybc = yb - yb.mean(0)
    return g, (xbc * ybc).sum() / (xbc ** 2).sum(), (xw * yw).sum() / (xw ** 2).sum(), pooled(S, F, Dg)

if sys.argv[1] != "report":
    nm = sys.argv[1]; Ds, lo, hi = SAMPLES[nm]; t0 = time.time()
    G = gridcache.get()
    F, Dg = predict(G, Ds, lo, hi)
    Fl, _ = predict(G, Ds, lo, hi, off_dL=0.3)
    np.savez(os.path.join(HERE, 'alt_%s.npz' % nm), F=F, Dg=Dg, Flp3=Fl, D=np.array(Ds))
    print(nm, len(Ds), 'done %.0fs' % (time.time() - t0))
    if nm in ('B', 'C'):
        idx = [P['l'].index(l) for l in LS]
        print('   check against the frozen predictions: max |full - frozen| %.2e, max |diag - frozen| %.2e' % (
            np.abs(F - np.array(P[nm]['full'])[:, idx]).max(), np.abs(Dg - np.array(P[nm]['diag'])[:, idx]).max()))
    sys.exit()

# ---------------- report ----------------
z = np.load(os.path.join(HERE, 'r120_100_600.npz')); Dref = list(z['D'])
Xr = design(Dref); Sr = z['S']; Fr = z['F']
Ar = np.linalg.pinv(Xr.T @ Xr) @ Xr.T
resid_ref = (Sr - Fr) - Xr @ (Ar @ (Sr - Fr))                       # residual rows, rounds 117/120 at 100-600
nref = np.mean([nzeros(D, 100, 600) for D in Dref])
zb = load_B(); XB = design(zb[0]); AB = np.linalg.pinv(XB.T @ XB) @ XB.T
resid_B = (zb[1] - zb[2]) - XB @ (AB @ (zb[1] - zb[2])); nB = np.mean([nzeros(D, 100, 600) for D in zb[0]])
print('noise model: residual rows of rounds 117/120 at 100-600 (mean %.0f zeros per set) and of B (%.0f), scaled by sqrt(n_ref/n_new)' % (nref, nB))
print('%-5s %4s %6s %6s | %-34s | %-34s | %s' % ('', 'n', 'meanL', 'zeros', 'AMP 0.932: gamma between within pooled', 'L+0.3: gamma between within pooled', 'expected se of gamma (ref R120 / ref B)'))
rows = {}
for nm, (Ds, lo, hi) in SAMPLES.items():
    f = os.path.join(HERE, 'alt_%s.npz' % nm)
    if not os.path.exists(f): print(nm, 'not computed'); continue
    a = np.load(f); F, Dg, Fl = a['F'], a['Dg'], a['Flp3']; X = design(Ds); n = len(Ds)
    Lm = np.mean([float(np.dot(*L_nodes(abs(D), lo, hi, 8))) for D in Ds]); nz = np.mean([nzeros(D, lo, hi) for D in Ds])
    amp = readings(Dg + 0.932 * (F - Dg), F, Dg, X); lp3 = readings(Fl, F, Dg, X)
    ses = []
    for res, nr in ((resid_ref, nref), (resid_B, nB)):
        rng = np.random.default_rng(5); sc = math.sqrt(nr / nz); gs = []
        for _ in range(400):
            gs.append(gamma(F + sc * res[rng.integers(0, len(res), n)], F, Dg, X))
        ses.append(float(np.std(gs)))
    rows[nm] = dict(n=n, L=Lm, zeros=nz, amp=amp, lp3=lp3, se=ses)
    print('%-5s %4d %6.3f %6.0f | %.3f %.3f %.3f %.3f            | %.3f %.3f %.3f %.3f            | %.3f / %.3f' % ((nm, n, Lm, nz) + tuple(amp) + tuple(lp3) + tuple(ses)))
json.dump({k: {kk: (list(map(float, vv)) if isinstance(vv, (tuple, list)) else float(vv)) for kk, vv in v.items()} for k, v in rows.items()}, open(os.path.join(HERE, 'alt_report.json'), 'w'), indent=1)
