"""Step 4 of the review: twins split by the residue class b of p mod q (q = 3, 4, 5, 7, 105), holdout (2.0e9, 2.0e10].
Class-level form of the claim (local densities in classes, the heuristic of section 2.1 applied to one class b):
   E_b(x) = R_b(x) - S x / A_q  ~  (S phi(q) / A_q) * [ Delta(x; q, b) + Delta(x; q, b + 2) ],
   Delta(x; q, c) = theta(x; q, c) - x / phi(q),  A_q = number of admissible classes b (b, b+2 units mod q).
Summed over b this gives the character formula; per class, odd characters do NOT cancel.
Fit: per class b, WLS of the cell increments of E_b/sqrt(x) on the increments of Delta(x; q, c)/sqrt(x) for every unit c mod q
(free coefficient per c) + 4 nuisance columns.  Prediction: S phi/A on c = b and c = b + 2, 0 on the other c.
Alternatives: (T) class twins follow the total prime error only: coefficient 2S/A_q on every c;
              (H) 'each member hears its own class with the full weight 2S phi/A of one member' is not a consistent model; not used.
usage: rc_class.py rc_full.npz
"""
import sys, math
import numpy as np

D = np.load(sys.argv[1])
X = D['X'].astype(np.float64)
K = np.arange(2, 4717)
Xk = X[K]; dx = X[K] - X[K - 1]; u = np.log(Xk); sq = np.sqrt(Xk); w = 1 / u ** 2
thc = D['thc'][K]; twc = D['twc'][K]
M = 420
S = 2 * 0.6601618158468695739
uc = (u - u.mean()) / (u.max() - u.min())
NUI = np.column_stack([np.ones_like(u), uc, uc ** 2, np.exp(-(u - u[0]) / 2)])

def wls(y, Z, blocks=40):
    sw = np.sqrt(w); A = Z * sw[:, None]; b = y * sw
    G = np.linalg.inv(A.T @ A); beta = G @ (A.T @ b); res = b - A @ beta
    s2 = res @ res / (len(b) - Z.shape[1])
    n = len(y); edges = np.linspace(0, n, blocks + 1).astype(int); bs = []
    for g in range(blocks):
        m = np.ones(n, bool); m[edges[g]:edges[g + 1]] = False
        bs.append(np.linalg.solve(A[m].T @ A[m], A[m].T @ b[m]))
    bs = np.array(bs); cj = (blocks - 1) / blocks * ((bs - bs.mean(0)).T @ (bs - bs.mean(0)))
    return beta, s2 * G, cj

def classsum(arr, q, c):
    return arr[:, [k for k in range(M) if k % q == c]].sum(1)

for q in (3, 4, 5, 7, 105):
    unitsq = [c for c in range(q) if math.gcd(c, q) == 1]
    phi = len(unitsq)
    adm = [b for b in range(q) if math.gcd(b * (b + 2), q) == 1]
    A = len(adm)
    Dl = np.column_stack([(classsum(thc, q, c) - dx / phi) / sq for c in unitsq])
    pred_own = S * phi / A
    print('\n=== q = %d: %d admissible classes b, phi = %d; predicted coefficient %.4f on Delta(b) and Delta(b+2), 0 elsewhere; '
          'alternative T: %.4f on every class' % (q, A, phi, pred_own, 2 * S / A))
    own_est, oth_est, chi2_tot, dof = [], [], 0.0, 0
    for b in adm:
        y = (classsum(twc, q, b) - S / A * dx) / sq
        if q == 105:
            # collapse: own = Delta(b) + Delta(b+2), other = sum of the other 46 classes
            iown = [unitsq.index(b), unitsq.index((b + 2) % q)]
            own = Dl[:, iown].sum(1); oth = Dl.sum(1) - own
            Z = np.column_stack([own, oth, NUI])
            be, cv, cj = wls(y, Z)
            own_est.append((be[0], cv[0, 0], cj[0, 0])); oth_est.append((be[1], cv[1, 1], cj[1, 1]))
            d = np.array([be[0] - pred_own, be[1]]); chi2_tot += d @ np.linalg.solve(cv[:2, :2], d); dof += 2
            continue
        Z = np.column_stack([Dl, NUI])
        be, cv, cj = wls(y, Z)
        pred = np.array([pred_own if c in (b, (b + 2) % q) else 0.0 for c in unitsq])
        d = be[:phi] - pred
        chi2 = d @ np.linalg.solve(cv[:phi, :phi], d)
        dT = be[:phi] - 2 * S / A
        chi2T = dT @ np.linalg.solve(cv[:phi, :phi], dT)
        chi2_tot += chi2; dof += phi
        print('  b=%3d: ' % b + '  '.join('c=%d %.3f+-%.3f' % (c, be[i], math.sqrt(cv[i, i])) for i, c in enumerate(unitsq)) +
              '   chi2 vs claim %.1f / %d, vs T %.1f' % (chi2, phi, chi2T))
        for i, c in enumerate(unitsq):
            (own_est if c in (b, (b + 2) % q) else oth_est).append((be[i], cv[i, i], cj[i, i]))
    def pool(lst):
        v = np.array(lst); wt = 1 / v[:, 1]
        m = (v[:, 0] * wt).sum() / wt.sum()
        return m, math.sqrt(1 / wt.sum())
    mo, so = pool(own_est); mt, st = pool(oth_est) if oth_est else (float('nan'), float('nan'))
    print('  pooled (inverse-variance, classes treated as independent): own %.4f +- %.4f (pred %.4f), other %.4f +- %.4f (pred 0); '
          'chi2 of all coefficients vs claim %.1f on %d' % (mo, so, pred_own, mt, st, chi2_tot, dof))

# Odd characters per class mod 5: the claim predicts a nonzero odd-character coefficient in each class, cancelling in the sum.
print('\n=== odd characters mod 5 class by class (twins with p = b mod 5); regressors theta(x, chi)/sqrt(x) for chi mod 5')
g = {1: 0, 2: 1, 4: 2, 3: 3}                       # discrete log base 2 mod 5
chis = {j: np.array([0] + [np.exp(2j * np.pi * j * g[c] / 4) for c in (1, 2, 3, 4)]) for j in range(4)}
TH5 = {j: sum(chis[j][c] * classsum(thc, 5, c) for c in (1, 2, 3, 4)) / sq for j in range(4)}
r = (TH5[0] - dx) / sq * sq / sq                    # theta(x) - x  (all p > 5)
r = (classsum(thc, 5, 1) + classsum(thc, 5, 2) + classsum(thc, 5, 3) + classsum(thc, 5, 4) - dx) / sq
Z = np.column_stack([r, TH5[2].real, 2 * TH5[1].real, -2 * TH5[1].imag, NUI])   # zeta, chi5 (j=2, real even), odd pair j=1 (Re, Im)
tot = np.zeros(Z.shape[1]); totp = np.zeros(4)
for b in (1, 2, 4):
    y = (classsum(twc, 5, b) - S / 3 * dx) / sq
    be, cv, cj = wls(y, Z)
    # prediction: (S/3) [ (theta - x) * 2 + sum_{chi != 1} (conj chi(b) + conj chi(b+2)) theta(x, chi) ]  (characters mod 5)
    pz = 2 * S / 3
    p5 = S / 3 * (np.conj(chis[2][b]) + np.conj(chis[2][(b + 2) % 5])).real
    p1 = S / 3 * (np.conj(chis[1][b]) + np.conj(chis[1][(b + 2) % 5]))
    pred = [pz, p5, p1.real, p1.imag]
    print('  b=%d: zeta %.3f+-%.3f (pred %.3f) | chi5 %.3f+-%.3f (pred %.3f) | odd pair Re %.3f+-%.3f (pred %.3f) Im %.3f+-%.3f (pred %.3f)'
          % (b, be[0], math.sqrt(cv[0, 0]), pred[0], be[1], math.sqrt(cv[1, 1]), pred[1], be[2], math.sqrt(cv[2, 2]), pred[2], be[3], math.sqrt(cv[3, 3]), pred[3]))
    tot += be; totp += np.array(pred)
print('  sum over the three classes: zeta %.3f (pred %.3f), chi5 %.3f (pred %.3f), odd pair Re %.3f (pred %.3f), Im %.3f (pred %.3f)'
      % (tot[0], totp[0], tot[1], totp[1], tot[2], totp[2], tot[3], totp[3]))
