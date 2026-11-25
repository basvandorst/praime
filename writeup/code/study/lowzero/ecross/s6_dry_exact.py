"""Exact expectations of the 'lines' synthetic (modulated Poisson, lines r = 2^a 3^b with the law's amplitudes, intensity
clipped at 0): E G at the registered lines and E X_ij (independent processes: E[G_i conj G_j] = E G_i conj E G_j), by
Gaussian smoothing of the intensity on a fine grid (FFT).  Compared with (a) the law restricted to p = 2, 3
(w w c c J23/4pi^2, frozen) and (b) the sampled synthetic X of ecross_test.py (run_<design>_lines).
usage: s6_dry_exact.py DESIGN"""
import sys, math, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ecross_test as E
des = sys.argv[1]
F = E.load_frozen(os.path.dirname(os.path.abspath(__file__)), des)
n = len(F['ids']); h = 0.002; top = 402.0
t = np.arange(0.5, top, h) + h / 2
W1, W2 = E.WIN; T = W2 - W1; tg = np.arange(W1, W2 + E.DT / 2, E.DT)
EG = {s: np.zeros((n, len(tg)), complex) for s in E.SIGMAS}; EL = np.zeros((n, len(F['lines'])), complex); clip = []
law_lines = [tuple(x) for x in F['lines']]
for k in range(n):
    N, w, c = int(F['N'][k]), int(F['w'][k]), float(F['c'][k])
    Lt = math.log(N) - 2 * math.log(2 * math.pi) + 2 * np.log(np.abs(1 + 1j * t)); dens = Lt / (2 * math.pi)
    M = np.ones(len(t), complex)
    for p in (2, 3):
        M *= E._msum(t, F['M%d_v' % p], F['M%d' % p][k], p)
    th = E.theta_E(N, t)
    dens = dens + (1 / math.pi) * np.real(-w * c * M * np.exp(-2j * th)); clip.append(np.mean(dens < 0)); dens = np.maximum(dens, 0)
    f = dens * h * np.exp(2j * th)
    b = (t > E.BAND[0]) & (t <= E.BAND[1])
    for m, (a, bb) in enumerate(law_lines):
        EL[k, m] = np.sum(f[b] * np.exp(-1j * math.log(a / bb) * t[b]))
    for s_ in E.SIGMAS:
        s = s_ / math.sqrt(2)
        kern_t = np.arange(-int(8 * s / h), int(8 * s / h) + 1) * h
        g = np.convolve(f, np.exp(-kern_t ** 2 / (2 * s * s)), mode='same')
        EG[s_][k] = np.interp(tg, t, g.real) + 1j * np.interp(tg, t, g.imag)
print('%s: clipped fraction of t: mean %.4f max %.4f' % (des, np.mean(clip), np.max(clip)))
I, J = F['I'], F['J']
lawL = F['lineamp'] * (E.BAND[1] - E.BAND[0]) / (2 * math.pi)
is23 = np.array([all(q in (2, 3) for q in _f(a) + _f(b)) if True else False for (a, b) in law_lines]) if False else None
# lines r = 2^a 3^b of the registered list: 1, 2, 1/2, 3, 1/3, 4, 1/4, 3/2, 2/3
sel = [m for m, (a, b) in enumerate(law_lines) if all(x in (1, 2, 3, 4) for x in (a, b))]
d = (EL[:, sel].real - lawL[:, sel]) / np.abs(lawL[:, sel]).mean()
print('E G (synthetic, exact) against the law at lines r = 2^a 3^b: max |diff|/mean|law| %.4f; slope %.4f' % (np.abs(d).max(), np.sum(EL[:, sel].real * lawL[:, sel]) / np.sum(lawL[:, sel] ** 2)))
for s in E.SIGMAS:
    G = EG[s]; X = (G @ G.conj().T) * (E.DT / (math.sqrt(math.pi) * (s / math.sqrt(2)) * T))
    ex = X[I, J].real; law23 = F['w'][I] * F['w'][J] * F['c'][I] * F['c'][J] * F['J23_%.1f' % s] / (4 * math.pi ** 2)
    print('sigma %.1f: E X (exact, synthetic) vs law23: slope %.4f, max |diff| %.4f (law23 sd %.4f), corr %.5f; E Im X max %.4f' % (
        s, np.dot(ex, law23) / np.dot(law23, law23), np.abs(ex - law23).max(), law23.std(), np.corrcoef(ex, law23)[0, 1], np.abs(X[I, J].imag).max()))
    fx = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'run_%s_lines' % des, 'X_%.1f.npy' % s)
    if os.path.exists(fx):
        Xs = np.load(fx)[I, J].real
        A = np.column_stack([np.ones(len(ex)), ex]); co = np.linalg.lstsq(A, Xs, rcond=None)[0]
        A2 = np.column_stack([np.ones(len(ex)), law23]); co2 = np.linalg.lstsq(A2, Xs, rcond=None)[0]
        print('   sampled synthetic X: slope on E X %.4f (icpt %.4f); slope on law23 %.4f; residual sd %.4f' % (co[1], co[0], co2[1], np.std(Xs - A @ co)))
        np.save(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'run_%s_lines' % des, 'EX_%.1f.npy' % s), ex)
