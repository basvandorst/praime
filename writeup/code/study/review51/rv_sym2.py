# L(1, Sym^2 E) (analytic normalisation) by PARI lfunsympow at s = 2, against the direct partial Euler product
# (good p: [(1-a^2/p)(1-1/p)(1-b^2/p)]^-1 with a b = 1, a + b = lam_p; multiplicative p: (1 - 1/p^2)^-1), and B = c L(1,Sym2).
import os, numpy as np, cypari2, time
SC = os.path.dirname(os.path.abspath(__file__))+'/'
CU = np.load(SC + 'rv_curves.npz'); LW = np.load(SC + 'rv_law.npz')
PF = np.load(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/ecross/pred_fresh.npz')
fz = {s: i for i, s in enumerate(PF['ids'])}
pari = cypari2.Pari(); pari.allocatemem(800 * 10**6)
pr = CU['primes'].astype(float)
for k in np.linspace(0, len(CU['ids']) - 1, 10).astype(int):
    t0 = time.time()
    s = CU['ids'][k]; N = int(CU['N'][k])
    a = [int(x.replace('m', '-')) for x in s.split('_')[1:]]
    Lp = float(pari('lfun(lfunsympow(ellinit(%s), 2), 2)' % str(a)))
    lam = CU['ap'][k] / np.sqrt(pr); bad = (N % CU['primes'] == 0)
    sym = lam * lam - 1        # a^2 + 1 + b^2
    loc = np.where(bad, 1 / (1 - 1 / pr ** 2), 1 / ((1 - 1 / pr) * (1 - (sym - 1) / pr + 1 / pr ** 2)))
    # (1-a^2 X)(1-b^2 X) = 1 - (a^2+b^2) X + X^2, a^2+b^2 = lam^2 - 2 = sym - 1
    direct = float(np.exp(np.sum(np.log(loc))))
    j = fz[s]
    print('%-22s N=%4d  PARI L(1,Sym2) %.5f  direct(p<=5e5) %.5f  frozen LS1 %.5f | c own %.5f frozen %.5f | B own %.5f frozen %.5f  (%.0fs)' % (
        s, N, Lp, direct, PF['LS1'][j], LW['c'][k], PF['c'][j], LW['c'][k] * Lp, PF['B'][j], time.time() - t0))
