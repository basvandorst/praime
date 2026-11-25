import numpy as np
A = np.load('r24_leak_ext.npz'); B = np.load('r25_rep_leak.npz')
assert np.allclose(A['Ts'] if 'Ts' in A.files else B['Ts'], B['Ts'])
def win(O): return np.diff(O, axis=1)          # windows (10,20], ..., (140,150]: 14
rng = np.random.default_rng(152)
def chi2(W): n = len(W); m = W.mean(0); s = W.std(0, ddof=1) / np.sqrt(n); return np.sum((m / s) ** 2)
def null_chi2(W, K=4000):
    v = np.array([chi2(W * rng.choice([-1, 1], size=(len(W), 1))) for _ in range(K)]); return v
for lab, (Oo, Lo, Or, Lr) in (('raw', (A['OB'], 0 * A['LE'], B['OB'], 0 * B['LE'])), ('minus leak', (A['OB'], A['LE'], B['OB'], B['LE']))):
    Wo = win(Oo - Lo); Wr = win(Or - Lr)
    for part, f in (('Re', np.real), ('Im', np.imag)):
        wo, wr = f(Wo), f(Wr)
        r = np.corrcoef(wo.mean(0), wr.mean(0))[0, 1]
        nr = np.array([np.corrcoef(wo.mean(0), (wr * rng.choice([-1, 1], size=(len(wr), 1))).mean(0))[0, 1] for _ in range(4000)])
        c_obs = chi2(wr); nc = null_chi2(wr)
        print('%-10s %s  corr(old, rep) %+.3f  null p(corr >= obs) %.4f | rep window chi2 %.1f, sign-flip null p %.4f | old chi2 %.1f p %.4f' % (
            lab, part, r, (nr >= r).mean(), c_obs, (nc >= c_obs).mean(), chi2(wo), (null_chi2(wo) >= chi2(wo)).mean()))
print('rep window means Re raw', np.round(win(B['OB']).real.mean(0), 2))
print('old window means Re raw', np.round(win(A['OB']).real.mean(0), 2))
print('rep window means Im raw', np.round(win(B['OB']).imag.mean(0), 2))
print('old window means Im raw', np.round(win(A['OB']).imag.mean(0), 2))
print('rep leak mean Im      ', np.round(win(B['LE']).imag.mean(0), 2))
