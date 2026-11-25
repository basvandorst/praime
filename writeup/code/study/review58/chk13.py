import json, numpy as np, cypari2
z = np.load('rep/R_N1403_13.npz'); info = json.loads(str(z['info'])); g = z['g']
print({k: info[k] for k in ('eps', 'm0', 'A0', 'c_exp', 'n_first', 'added', 'added_extra', 'duplicates', 'n', 'drift', 'r_mean', 'status', 'retry_eps_fixed', 'c', 'c_end')})
pari = cypari2.Pari(); pari.allocatemem(500 * 10 ** 6)
L = pari('lfungenus2([(-1)*x^2+(-1)*x^3+(1)*x^4+(1)*x^5, x^0])')   # NO eps fix, independent init on [1, 63]
Li = pari.lfuninit(L, [1, 63]); h = np.sort(np.array([float(x) for x in pari.lfunzeros(Li, [0.01, 60])]))
gg = g[(g > 0) & (g <= 60)]
print('zeros <= 60: retry %d, independent %d, max |diff| %.2e' % (len(gg), len(h), np.max(np.abs(gg - h)) if len(gg) == len(h) else -1))
print('first zeros', np.round(gg[:5], 6), np.round(h[:5], 6))
