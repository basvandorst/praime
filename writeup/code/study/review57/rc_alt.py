"""Step 3 of the review: pseudo-data for alternatives, run through my estimator (rc_fit.py model: zeta + 95 characters mod 420 + nuisance),
and a shifted-regressor null for the errors.  Inputs: rc_full.npz (my recount) and rc_fit_out.npz (regressors, residuals).
Pseudo-data are built from the REAL prime regressors (explicit formula = real primes) plus noise:
  noise N1 = real twin residuals, circularly shifted by a random offset (keeps the real, sub-Poisson spectrum);
  noise N2 = Poisson counts at cell level (plain noise, no rigidity).
Alternatives (coefficient on theta(x) - x, and on theta(x, chi)):
  A0 independence        0                   no characters
  A1 one member (thin)   S                   S * c1(chi)  (only the member n; odd chi do NOT cancel: chi_-3 gets -S)
  A2 claim               2S                  S (c1 + c2)
  A3 no L-lines          2S                  0
  A4 original guess      2S                  claim + chi_-3 at -2S
  A5 coefficient 4 (mean of 2S(h) over even h; 'no singular-series dependence')   4, characters as claim
usage: rc_fit.py must have been run first; rc_alt.py rc_full.npz rc_fit_out.npz
"""
import sys, math, itertools
import numpy as np
rng = np.random.default_rng(155)
D = np.load(sys.argv[1]); F = np.load(sys.argv[2])
r, CH, NUI, w, res2 = F['r'], F['CH'], F['NUI'], F['w'], F['res2']
X = D['X'].astype(np.float64); K = np.arange(2, 4717); Xk = X[K]; dx = X[K] - X[K - 1]; sq = np.sqrt(Xk); u = np.log(Xk)
M = 420; S = 2 * 0.6601618158468695739

# --- rebuild the character bookkeeping (same construction as rc_fit.py) ---
def dlog_table(p, g):
    t = {}; x = 1
    for k in range(p - 1):
        t[x] = k; x = x * g % p
    return t
tabs = {3: dlog_table(3, 2), 5: dlog_table(5, 2), 7: dlog_table(7, 3)}
def comp_val(p, j, c):
    if p == 4:
        if c % 2 == 0: return 0
        return 1 if j == 0 else (1 if c % 4 == 1 else -1)
    if c % p == 0: return 0
    return np.exp(2j * np.pi * j * tabs[p][c % p] / (p - 1))
chars = []
for j4, j3, j5, j7 in itertools.product(range(2), range(2), range(4), range(6)):
    vals = np.array([comp_val(4, j4, c) * comp_val(3, j3, c) * comp_val(5, j5, c) * comp_val(7, j7, c) for c in range(M)], dtype=complex)
    cond = (4 if j4 else 1) * (3 if j3 else 1) * (5 if j5 else 1) * (7 if j7 else 1)
    chars.append(dict(j=(j4, j3, j5, j7), v=vals, q=cond))
nonpr = [i for i, ch in enumerate(chars) if ch['q'] > 1]
colinfo = []; seen = set()
for i in nonpr:
    if i in seen: continue
    jj = [k for k in nonpr if np.allclose(chars[k]['v'], np.conj(chars[i]['v']))][0]
    seen.add(i); seen.add(jj)
    if jj == i: colinfo.append((i, 're'))
    else: colinfo.append((i, 're')); colinfo.append((i, 'im'))
assert len(colinfo) == CH.shape[1]
adm = [n for n in range(M) if math.gcd(n * (n + 2), M) == 1]
def cmean(ch, shift):
    return np.mean([np.conj(ch['v'][(n + shift) % M]) for n in adm])
def vec(f):
    out = []
    for i, part in colinfo:
        b = f(chars[i]); out.append(b.real if part == 're' else b.imag)
    return np.array(out)
claim = vec(lambda ch: S * (cmean(ch, 0) + cmean(ch, 2)))
one = vec(lambda ch: S * cmean(ch, 0))
k3 = [k for k, (i, p) in enumerate(colinfo) if chars[i]['j'] == (0, 1, 0, 0)][0]
k4 = [k for k, (i, p) in enumerate(colinfo) if chars[i]['j'] == (1, 0, 0, 0)][0]
k5 = [k for k, (i, p) in enumerate(colinfo) if chars[i]['j'] == (0, 0, 2, 0)][0]
guess = claim.copy(); guess[k3] = -2 * S
alts = {'A0 independence': (0.0, 0 * claim), 'A1 one member': (S, one), 'A2 claim': (2 * S, claim),
        'A3 no L-lines': (2 * S, 0 * claim), 'A4 orig. guess': (2 * S, guess), 'A5 coef 4': (4.0, claim)}
print('one-member model: chi_-3 coefficient %.4f, chi_-4 %.4f, chi_5 %.4f' % (one[k3], one[k4], one[k5]))

Z = np.column_stack([r, CH, NUI]); sw = np.sqrt(w); A = Z * sw[:, None]; G = np.linalg.inv(A.T @ A)
def fit(y):
    b = G @ (A.T @ (y * sw)); res = y * sw - A @ b
    s2 = res @ res / (len(y) - Z.shape[1]); return b, np.sqrt(s2 * np.diag(G))
def loud(b):
    pres = np.abs(claim) > 1e-12
    return (b[1:1 + CH.shape[1]][pres] @ claim[pres]) / (claim[pres] @ claim[pres])

print('\n%-18s %-8s %10s %10s %10s %10s | separation from the claim in units of the measured error (zeta / chi_-3 / loudness)' % ('alternative', 'noise', 'zeta', 'chi_-3', 'chi_5', 'loudness'))
rows = {}
for name, (cz, cch) in alts.items():
    for noise in ('N1', 'N2'):
        zs = []
        for t in range(20):
            if noise == 'N1':
                nz = np.roll(res2, rng.integers(300, len(res2) - 300))
            else:
                lam = S * dx / u ** 2               # expected number of pairs per cell, each of weight ~ (ln x)^2
                nz = (rng.poisson(lam) - lam) * u ** 2 / sq
            y = cz * r + CH @ cch + nz
            b, se = fit(y)
            zs.append((b[0], b[1 + k3], b[1 + k5], loud(b), se[0], se[1 + k3]))
        zs = np.array(zs); m = zs.mean(0); sd = zs.std(0)
        rows[(name, noise)] = m
        print('%-18s %-8s %6.3f+-%.3f %6.3f+-%.3f %6.3f+-%.3f %6.3f+-%.3f' % (name, noise, m[0], sd[0], m[1], sd[1], m[2], sd[2], m[3], sd[3]))
real_b, real_se = fit((D['P2'][K] - S * dx) / sq)
print('\nREAL twins: zeta %.4f +- %.4f, chi_-3 %.4f +- %.4f, chi_5 %.4f, loudness(q|105) %.4f' % (real_b[0], real_se[0], real_b[1 + k3], real_se[1 + k3], real_b[1 + k5], loud(real_b)))
for name in alts:
    m = rows[(name, 'N1')]
    print('  distance of REAL from %-18s: zeta %+6.1f sigma, chi_-3 %+6.1f sigma' % (name, (real_b[0] - m[0]) / real_se[0], (real_b[1 + k3] - m[1]) / real_se[1 + k3]))

# shifted-regressor null: real twins against the prime regressors circularly shifted by s cells
print('\nshifted-regressor null (real twin series, all prime regressors rolled by s cells):')
y = (D['P2'][K] - S * dx) / sq
vals = []
for s in range(200, 4600, 200):
    Zs = np.column_stack([np.roll(r, s), np.roll(CH, s, axis=0), NUI])
    As = Zs * sw[:, None]; bs = np.linalg.solve(As.T @ As, As.T @ (y * sw))
    res = y * sw - As @ bs; s2 = res @ res / (len(y) - Zs.shape[1]); se = np.sqrt(s2 * np.diag(np.linalg.inv(As.T @ As)))
    vals.append((bs[0], se[0], bs[1 + k3], se[1 + k3], bs[1 + k5], se[1 + k5]))
vals = np.array(vals)
print('  %d shifts: zeta mean %.3f sd %.3f (mean reported error %.3f); chi_-3 sd %.3f (err %.3f); chi_5 sd %.3f (err %.3f)'
      % (len(vals), vals[:, 0].mean(), vals[:, 0].std(), vals[:, 1].mean(), vals[:, 2].std(), vals[:, 3].mean(), vals[:, 4].std(), vals[:, 5].mean()))
