# Pool the prime-side resonance: (a) the 8 curves to 1e7 at t_x = 400 and at x = 1e7;
# (b) all 138 curves at x = 5e5 (a_p from rv_curves), slope of the projected sum on the law, jackknife over curves.
import os, glob
import numpy as np
from scipy.special import beta
SC = os.path.dirname(os.path.abspath(__file__))+'/'
rows = []
for f in sorted(glob.glob(SC + 'rv_prime_N*.npy')):
    d = np.load(f, allow_pickle=True).item()
    res = d['res']; N, c = d['N'], d['c']
    law = lambda x: (2 * np.sqrt(2) / 3) * c * (x ** 3 / N) ** 0.25
    i400 = np.argmin(np.abs(res[:, 1] - 400)); il = len(res) - 1
    rows.append([res[i400, 1], res[i400, 2], res[i400, 3], law(res[i400, 0]), res[il, 2], res[il, 3], law(res[il, 0]),
                 d['Gs'].real / d['lawG'], d['Rs'].real / d['lawG'], d['Gh'].real / (d['lawG'] / 2), d['Rh'].real / (d['lawG'] / 2)])
rows = np.array(rows)
print('8 curves: t_x nearest 400 per curve:', np.round(rows[:, 0], 0))


def pooled(r, L):
    b = np.sum(r * L * L) / np.sum(L * L)
    n = len(r)
    bk = np.array([(np.sum(r * L * L) - r[k] * L[k] ** 2) / (np.sum(L * L) - L[k] ** 2) for k in range(n)])
    return b, np.sqrt((n - 1) / n * np.sum((bk - bk.mean()) ** 2))


print('  pooled Re(S0/law) at t_x~400: %.3f +- %.3f ; Riesz k=1: %.3f +- %.3f' % (pooled(rows[:, 1], rows[:, 3]) + pooled(rows[:, 2], rows[:, 3])))
print('  pooled Re(S0/law) at x = 1e7: %.3f +- %.3f ; Riesz k=1: %.3f +- %.3f' % (pooled(rows[:, 4], rows[:, 6]) + pooled(rows[:, 5], rows[:, 6])))
print('  same 8 curves, zero side: Re G/law sharp %.3f +- %.3f, primes band sharp %.3f +- %.3f, zeros Hann %.3f +- %.3f, primes Hann %.3f +- %.3f' % (
    pooled(rows[:, 7], np.ones(8)) + pooled(rows[:, 8], np.ones(8)) + pooled(rows[:, 9], np.ones(8)) + pooled(rows[:, 10], np.ones(8))))

# (b) all 138 curves at x = 5e5
CU = np.load(SC + 'rv_curves.npz'); LW = np.load(SC + 'rv_law.npz')
primes = CU['primes'].astype(float); AP = CU['ap'].astype(float); N = CU['N']; w = CU['w']; c = LW['c']
X = 5e5
proj, lawv, proj1 = [], [], []
for k in range(len(N)):
    lam = AP[k] / np.sqrt(primes); bad = (N[k] % CU['primes'] == 0)
    th = np.arccos(np.clip(lam / 2, -1, 1))
    ns, Ls = [], []
    for kk in range(1, 20):
        m = primes ** kk <= X
        if not m.any():
            break
        coef = lam[m] if kk == 1 else np.where(bad[m], lam[m] ** kk, 2 * np.cos(kk * th[m]))
        ns.append(primes[m] ** kk); Ls.append(coef * np.log(primes[m]))
    n = np.concatenate(ns); L = np.concatenate(Ls)
    ph = np.exp(-4j * np.pi * np.sqrt(n / N[k]))
    S0 = np.sum(L * ph); S1 = np.sum(L * ph * (1 - n / X))
    u = np.exp(-3j * np.pi / 4) * w[k]                     # direction of the law
    lw = (2 * np.sqrt(2) / 3) * c[k] * (X ** 3 / N[k]) ** 0.25
    proj.append((S0 * np.conj(u)).real); proj1.append((S1 * np.conj(u)).real / (0.75 * beta(0.75, 2))); lawv.append(lw)
proj, proj1, lawv = map(np.array, (proj, proj1, lawv))
for nm, y in (('raw S0', proj), ('Riesz k=1', proj1)):
    b = np.sum(y * lawv) / np.sum(lawv ** 2)
    n = len(y)
    bk = np.array([(np.sum(y * lawv) - y[k] * lawv[k]) / (np.sum(lawv ** 2) - lawv[k] ** 2) for k in range(n)])
    print('138 curves, primes only, x = 5e5 (t_x = %.0f..%.0f): %s slope on the law %.3f +- %.3f' % (
        2 * np.pi * np.sqrt(X / N.max()), 2 * np.pi * np.sqrt(X / N.min()), nm, b, np.sqrt((n - 1) / n * np.sum((bk - bk.mean()) ** 2))))
    # c_E dependence from primes only
    A = np.c_[lawv, lawv / c]
    th_ = np.linalg.lstsq(A, y, rcond=None)[0]
    print('     [law with c_E, law with c_E -> 1] coefficients: %.3f, %.3f' % tuple(th_))
