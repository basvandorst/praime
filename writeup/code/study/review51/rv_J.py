# Reviewer's own pair integral J12 = int H12(iu) exp(-u^2/2 s^2) du, H12 = prod_p h_p(iu),
# h_p(iu) = sum_v rho1(p^v) rho2(p^v) p^-|v| p^{iuv}; exact local factors (v = -2..V) for all p <= 5e5.
import os, time, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/')
from rv_law import lam, bad, pr, primes, local_rho, ids, w

SC = os.path.dirname(os.path.abspath(__file__))+'/'
u = np.arange(0, 8.0001, 0.02)
lnp = np.log(pr)


def logH(i, j, Pmax=500000, V=40):
    m = pr <= Pmax
    out = np.zeros(len(u), complex)
    # small primes exact with long v range, large primes v in -2..3
    for (lo, hi, VV) in ((0, 3000, V), (3000, Pmax + 1, 3)):
        mm = m & (pr > lo) & (pr <= hi)
        p = pr[mm]; lp = lnp[mm]
        R1, _ = local_rho(lam[i, mm], bad[i, mm], p, VV)
        R2, _ = local_rho(lam[j, mm], bad[j, mm], p, VV)
        vs = np.arange(-2, VV + 1)
        coef = R1 * R2 * p[None, :] ** (-np.abs(vs)[:, None].astype(float))     # (nv, np)
        for b0 in range(0, len(p), 4000):
            sl = slice(b0, b0 + 4000)
            ph = np.exp(1j * np.outer(u, lp[sl]))          # (nu, nb)
            h = np.zeros((len(u), ph.shape[1]), complex)
            for a, v in enumerate(vs):
                h += coef[a, sl][None, :] * ph ** v
            out += np.log(h).sum(1)
    return out


def Jint(H, sig):
    K = np.exp(-u ** 2 / (2 * sig * sig))
    f = (H.real * K)
    return 2 * np.trapezoid(f, u)


if __name__ == '__main__':
    PF = np.load(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/ecross/pred_fresh.npz')
    idx = {s: k for k, s in enumerate(ids)}
    rng = np.random.default_rng(148)
    sel = rng.choice(len(PF['I']), 24, replace=False)
    # add the pairs with the largest |J| and largest |pred|
    sel = np.unique(np.concatenate([sel, np.argsort(-np.abs(PF['J_0.5']))[:3], np.argsort(-np.abs(PF['pred_0.5']))[:3]]))
    t0 = time.time()
    print(' pair             N1    N2    J0.5 own   frozen    J1.0 own   frozen   pred0.5 own  frozen   H(0) own  frozen')
    rows = []
    for q in sel:
        a, b = PF['ids'][PF['I'][q]], PF['ids'][PF['J'][q]]
        i, j = idx[a], idx[b]
        L = logH(i, j)
        H = np.exp(L)
        J5, J1 = Jint(H, 0.5), Jint(H, 1.0)
        from rv_law import cE
        p5 = w[i] * w[j] * cE(i) * cE(j) * J5 / (4 * np.pi ** 2)
        rows.append((J5, PF['J_0.5'][q], J1, PF['J_1.0'][q], p5, PF['pred_0.5'][q]))
        print('%5d %22s %5s %5s %10.4f %10.4f %10.4f %10.4f %10.5f %10.5f %8.4f %8.4f' % (q, '', a.split('_')[0], b.split('_')[0], J5, PF['J_0.5'][q], J1, PF['J_1.0'][q], p5, PF['pred_0.5'][q], H[0].real, PF['H0'][q]))
    rows = np.array(rows)
    print('max |dJ| / |J| (0.5): %.2e  (1.0): %.2e ; max |dpred| %.2e' % (np.max(np.abs(rows[:, 0] - rows[:, 1]) / np.abs(rows[:, 1])),
          np.max(np.abs(rows[:, 2] - rows[:, 3]) / np.abs(rows[:, 3])), np.max(np.abs(rows[:, 4] - rows[:, 5]))))
    print('time', time.time() - t0)
