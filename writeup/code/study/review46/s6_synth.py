"""(A1) Which part of the readings is the explicit formula below a cutoff NC?
dn_syn_i(t) = sum_{n = p^m <= NC} c_n(w) chi_i(n) cos((t + w/2) ln n),  c_n = -(2/pi) sin(w ln n/2)/(m sqrt(n))
(the explicit formula for S(t + w) - S(t), truncated; no zeros used). Pair covariances of the synthetic series on the
registered grid, of the real series, and of the remainder dn - dn_syn. Registered readings on each.
Also: the share of the predictor variance over pairs carried by primes above 300 / above 1.6e4 (Heisenberg at t = 100)."""
import sys, math, time
import numpy as np
sys.path.insert(0, __import__('os').path.dirname(__file__))
from mylib import *
from myreg import fit_all, jackknife

NC = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
t0 = time.time()
Z = load_zeros(); S1 = np.load(OUT + '/s1.npz'); S2 = np.load(OUT + '/s2_pairs.npz'); Ds = list(S1['Ds'])
Xk = np.load(OUT + '/kron.npy').astype(float); P = np.load(OUT + '/primes.npy')
I, J = S2['I'], S2['J']
PW = []
for j, p in enumerate(P[P <= NC]):
    m = 1
    while p ** m <= NC: PW.append((int(p) ** m, m, j)); m += 1
lnn = np.array([math.log(a) for a, _, _ in PW]); mm = np.array([b for _, b, _ in PW]); jj = np.array([c for _, _, c in PW])
CH = Xk[:, jj] ** mm[None, :]
print('NC %d: %d prime powers' % (NC, len(PW)), flush=True)
Csyn = np.zeros((16, len(Ds), len(Ds))); Crem = np.zeros_like(Csyn); Cx = np.zeros_like(Csyn)
for k, w in enumerate(W):
    tg = grid(w); real = np.array([dn(Z[D], D, w, tg) for D in Ds])
    cn = -(2 / math.pi) / (mm * np.exp(lnn / 2)) * np.sin(w * lnn / 2)
    coef = CH * cn[None, :]
    syn = np.zeros_like(real)
    for s in range(0, len(lnn), 400):
        B = np.cos(np.outer(tg + w / 2, lnn[s:s + 400]))
        syn += coef[:, s:s + 400] @ B.T
    rem = real - syn
    Csyn[k] = syn @ syn.T / len(tg); Crem[k] = rem @ rem.T / len(tg); Cx[k] = (syn @ rem.T + rem @ syn.T) / len(tg) / 2
    if k == 0 or k == 15:
        print('w=%.1f: var real %.3f, var syn %.3f, var rem %.3f, corr(real, syn) %.3f  [%.0fs]' % (w, real.var(1).mean(), syn.var(1).mean(), rem.var(1).mean(),
              np.mean([np.corrcoef(a, b)[0, 1] for a, b in zip(real, syn)]), time.time() - t0), flush=True)
args = (S2['PSI'], S2['REST'], S2['FULL'], S2['SW'], S1['ATAB'], S1['BTAB'], I, J)
for name, C in (('real', S1['C']), ('synthetic (explicit formula, n <= NC)', Csyn), ('remainder x remainder', Crem), ('2 x cross syn-rem', 2 * Cx)):
    Cp = C[:, I, J]
    full, se, coefs, cse, cr = jackknife(Cp, *args)
    print('%-40s ' % name + ' '.join('%s %+.4f(%.4f)' % (k.replace('gamma_', 'g'), full[k], se[k]) for k in ('gamma_2', 'gamma_3', 'gamma_5', 'gamma_7', 'gammaB_2', 'gamma_rest', 'gamma_full', 'small_w', 'gamma_sw')))
    print('    icpt mean %+.2e; residual sd per pair (w=0.3, 4.8): %.4f %.4f' % (full['icpt_mean'], *[np.std(Cp[k] - np.column_stack([np.ones(len(I))] + [S2['PSI'][:, a] for a in range(4)] + [S2['PSI'][:, a] ** 2 for a in range(4)] + [S2['REST'][k], S2['SW'][k]]) @ coefs[k]) for k in (0, 15)]))
np.savez(OUT + '/s6_synth_%d.npz' % NC, Csyn=Csyn, Crem=Crem, Cx=Cx)
# share of predictor variance from primes above thresholds
X = Xk; A = np.zeros((16, len(P))); B = np.zeros((16, len(P))); lp = np.log(P.astype(float))
for m in range(1, 40):
    c = 1.0 / (m * m * P.astype(float) ** m)
    t = c[None, :] * np.sin(m * W[:, None] * lp[None, :] / 2) ** 2
    if m % 2: A += t
    else: B += t
A *= 2 / math.pi ** 2; B *= 2 / math.pi ** 2
rng = np.random.default_rng(1); sel = rng.choice(len(I), 4000, replace=False)
for lo, hi, lab in ((11, 300, 'REST from 7 < p <= 300'), (300, 16000, '300 < p <= 1.6e4'), (16000, 200001, 'p > 1.6e4 (Heisenberg range at t >= 100)')):
    m = (P >= lo) & (P < hi); m2 = (P > 7)
    def pr(mask, k):
        Xa = X[I[sel]][:, mask]; Xb = X[J[sel]][:, mask]
        return ((Xa * Xb) * A[k, mask]).sum(1) + ((Xa ** 2 * Xb ** 2) * B[k, mask]).sum(1)
    sh = [pr(m, k).var() / pr(m2, k).var() for k in (0, 1, 8, 15)]
    shf = [pr(m, k).var() / pr(P > 0, k).var() for k in (0, 1, 8, 15)]
    print('variance share of %-42s in REST: w=0.3 %.4f, 0.6 %.4f, 2.7 %.4f, 4.8 %.4f | in FULL: %.4f %.4f %.4f %.4f' % (lab, *sh, *shf))
print('done %.0fs' % (time.time() - t0))
