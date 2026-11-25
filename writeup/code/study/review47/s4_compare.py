"""Compare, pair by pair, X from zeros and X from primes (no zeros), both with the Hann product window on [100, 594];
and both against the frozen law. Also: the xi-lines (major arcs) of the cross-spectrum, and the random-sign noise surrogate."""
import sys, math, glob, os
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import numpy as np
from rlib import *

D_ = __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/'
W1, W2 = 100.0, 594.0; NW = 3 / 8 * (W2 - W1)
P = np.load(FZ + 'pred_OO.npz'); Ds = [int(x) for x in P['Ds']]
pos = {(int(i), int(j)): k for k, (i, j) in enumerate(zip(P['I'], P['J']))}
F = [np.load(D_ + 'gr/GR_%d.npz' % k) for k in range(25) if os.path.exists(D_ + 'gr/GR_%d.npz' % k)]
nf = len(F); XI = F[0]['xi']; dxi = XI[1] - XI[0]
G = np.array([f['G'] for f in F]); R = np.array([f['R'] for f in F]); idx = [int(f['i']) for f in F]
print('%d functions; per-function rms|G - R| / rms|G| over xi in [-7.5, 7.5]: %s' % (nf, np.round(np.sqrt(np.mean(np.abs(G - R) ** 2, 1) / np.mean(np.abs(G) ** 2, 1)), 4)))


def xmat(A, sig, lo=-99, hi=99):
    w = np.exp(-sig ** 2 * XI ** 2 / 2) * dxi * ((XI >= lo) & (XI <= hi))
    w[0] *= 0.5; w[-1] *= 0.5
    return (A * w) @ A.conj().T * sig / math.sqrt(2 * math.pi) / NW


II, JJ = np.triu_indices(nf, 1)
law = {s: np.array([P['pred_%.1f' % s][pos[(min(idx[a]), max(idx[a]))]] if False else P['pred_%.1f' % s][pos[tuple(sorted((idx[a], idx[b])))]] for a, b in zip(II, JJ)]) for s in (0.5, 1.0)}
law1 = {s: np.array([P['pred1_%.1f' % s][pos[tuple(sorted((idx[a], idx[b])))]] for a, b in zip(II, JJ)]) for s in (0.5, 1.0)}
# sign convention: pred is symmetric for real pairs; X_ab with a<b in our order corresponds to X_ij with i = idx[a]
Z = load_group(R120F)
for sig in (0.5, 1.0):
    XZ = xmat(G, sig); XP = xmat(R, sig)
    xz, xp, L = XZ[II, JJ], XP[II, JJ], law[sig]
    # identity check: direct double sum with Hann product weights for 5 pairs
    chk = []
    for a, b in list(zip(II, JJ))[:5]:
        ga = Z['D%d' % Ds[idx[a]]]; gb = Z['D%d' % Ds[idx[b]]]
        ga = ga[(ga > W1) & (ga < W2)]; gb = gb[(gb > W1) & (gb < W2)]
        ha = np.sin(math.pi * (ga - W1) / (W2 - W1)) ** 2 * np.exp(2j * theta(ga, abs(Ds[idx[a]]), int(Ds[idx[a]] < 0)))
        hb = np.sin(math.pi * (gb - W1) / (W2 - W1)) ** 2 * np.exp(2j * theta(gb, abs(Ds[idx[b]]), int(Ds[idx[b]] < 0)))
        u = gb[None, :] - ga[:, None]
        d = np.sum(ha[:, None] * np.conj(hb)[None, :] * np.exp(-u ** 2 / (2 * sig ** 2))) / NW
        chk.append(abs(d - XZ[a, b]))
    print('\nsigma %.1f (%d pairs): identity check (xi-integral vs direct double sum, zeros): max diff %.1e' % (sig, len(xz), max(chk)))
    print('  rms X(zeros) %.4f; rms law %.4f; rms [X(zeros) - X(primes)] %.4f (Re %.4f); rms [Re X(zeros) - law] %.4f; rms Im X(zeros) %.4f' % (
        np.sqrt(np.mean(np.abs(xz) ** 2)), np.sqrt(np.mean(L ** 2)), np.sqrt(np.mean(np.abs(xz - xp) ** 2)), np.sqrt(np.mean((xz - xp).real ** 2)),
        np.sqrt(np.mean((xz.real - L) ** 2)), np.sqrt(np.mean(xz.imag ** 2))))
    print('  corr(Re X(zeros), Re X(primes)) = %.5f; corr(Re X(zeros), law) = %.4f; corr(Re X(primes), law) = %.4f' % (
        np.corrcoef(xz.real, xp.real)[0, 1], np.corrcoef(xz.real, L)[0, 1], np.corrcoef(xp.real, L)[0, 1]))
    one = np.ones(len(L))
    for nm, y in (('X(zeros)', xz.real), ('X(primes)', xp.real)):
        c, se, _, _ = ols_jack([one, L], y, II, JJ)
        c2, se2, _, _ = ols_jack([one, law1[sig], L - law1[sig]], y, II, JJ)
        print('  %-9s on law: slope %.4f +- %.4f (icpt %.4f); split [H=1, H-1]: %.4f +- %.4f, %.4f +- %.4f' % (nm, c[1], se[1], c[0], c2[1], se2[1], c2[2], se2[2]))
    c, se, _, _ = ols_jack([one, xp.real], xz.real, II, JJ)
    print('  Re X(zeros) on Re X(primes): slope %.5f +- %.5f, icpt %.5f' % (c[1], se[1], c[0]))
    # sharp-window X of the registered statistic for the same pairs
    Xs = np.load(D_ + 'Xown_OO_%.1f.npy' % sig)
    xs = np.array([Xs[idx[a], idx[b]] if idx[a] < idx[b] else np.conj(Xs[idx[b], idx[a]]) for a, b in zip(II, JJ)])
    print('  registered (sharp-window) X on these pairs: corr with Hann X(zeros) %.3f; rms [Re X_sharp - law] %.4f' % (np.corrcoef(xs.real, xz.real)[0, 1], np.sqrt(np.mean((xs.real - L) ** 2))))
    # contribution of the xi = 0 line (|xi| < 0.05) vs the H = 1 law, and of the rest vs H - 1
    X0z = xmat(G, sig, -0.05, 0.05)[II, JJ].real; X0p = xmat(R, sig, -0.05, 0.05)[II, JJ].real
    c, se, _, _ = ols_jack([one, law1[sig]], X0z, II, JJ)
    c1, se1, _, _ = ols_jack([law1[sig]], X0z, II, JJ)
    print('  xi-window |xi| < 0.05 (the xi = 0 line): Re X0 on the H = 1 law, no icpt: %.4f +- %.4f; zeros vs primes rms %.1e' % (c1[0], se1[0], np.sqrt(np.mean((X0z - X0p) ** 2))))
    # product of the two one-point sums predicts the xi=0 line pair by pair
    G0 = np.array([np.sum(g[np.abs(XI) < 1e-9]) for g in G])
    rest = xz.real - X0z
    c, se, _, _ = ols_jack([one, L - law1[sig]], rest, II, JJ)
    print('  rest (|xi| >= 0.05) on the H - 1 law: slope %.4f +- %.4f' % (c[1], se[1]))

# per-function one-point (xi = 0) check: G(0)/(law) with Hann: -kappa * NWlin/2pi, NWlin = int w = T/2
k0 = np.argmin(np.abs(XI))
kap = np.array([kappa(abs(Ds[i])) for i in idx])
print('\none point (Hann): Re G(0) / (-kappa T/4pi): mean %.4f sd %.4f; |G(0) - R(0)| max %.3f' % (np.mean(G[:, k0].real / (-kap * (W2 - W1) / 2 / (2 * math.pi))), np.std(G[:, k0].real / (-kap * (W2 - W1) / 2 / (2 * math.pi))), np.abs(G[:, k0] - R[:, k0]).max()))

# lines at xi = -ln 2 and +ln 2: integrate G_a conj G_b over |xi - x0| < 0.05, regress on psi(2) = chi_a(2) chi_b(2);
# predicted weights (H's Dirichlet coefficients): 2 psi(2) at xi = -ln 2, conj psi(2)/2 at xi = +ln 2 (times the H=1 scale)
chi2 = np.array([kron(Ds[i], 2) for i in idx]); psi2 = chi2[II] * chi2[JJ]
for x0, pred_w in ((-math.log(2), 2.0), (math.log(2), 0.5)):
    for nm, A in (('zeros', G), ('primes', R)):
        sig = 0.5
        Xl = xmat(A, sig, x0 - 0.05, x0 + 0.05)[II, JJ].real
        base = law1[sig] * math.exp(-sig ** 2 * x0 ** 2 / 2)       # H = 1 scale times K-hat at x0
        m = psi2 != 0
        c, se, _, _ = ols_jack([psi2[m] * base[m]], Xl[m], II[m], JJ[m])
        print('line at xi = %+.3f (%s): coefficient of psi(2) %.3f +- %.3f (H predicts %.2f); pairs with psi(2)=0: mean %.4f' % (x0, nm, c[0], se[0], pred_w, Xl[~m].mean() if (~m).any() else float('nan')))

# random-sign surrogates
RS = [np.load(f) for f in sorted(glob.glob(D_ + 'gr/GR_*_rand.npz'))]
if RS:
    Rr = np.array([f['R'] for f in RS]); kr = [int(f['i']) for f in RS]
    for sig in (0.5,):
        A = np.vstack([Rr, R])
        XX = xmat(A, sig); nr = len(Rr)
        ss = XX[:nr, :nr][np.triu_indices(nr, 1)]; sr = XX[:nr, nr:].ravel()
        print('\nsurrogates (random multiplicative signs), sigma %.1f: rms X(rand x rand) %.4f (%d pairs), rms X(rand x real) %.4f; mean %.4f / %.4f' % (
            sig, np.sqrt(np.mean(np.abs(ss) ** 2)), len(ss), np.sqrt(np.mean(np.abs(sr) ** 2)), ss.real.mean(), sr.real.mean()))
        print('  predicted per-pair noise sd for Re X: sqrt(2 Var(rand x real) - Var(rand x rand))/sqrt2-ish = %.4f' % math.sqrt(max(0, 2 * np.mean(sr.real ** 2) - np.mean(ss.real ** 2))))
