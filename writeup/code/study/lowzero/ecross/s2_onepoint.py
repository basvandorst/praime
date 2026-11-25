"""Step 2: the one-point sums G_E = sum_{20 < gamma <= 400} e^{2 i theta_E(gamma)} (sharp and Hann sin^2 taper) on all
158 zero sets (20 development + 138 fresh), against the recipe law -w c_E (T2 - T1)/2pi (Hann: half of it) and the
alternatives: zero (H0), no root number (-c_E T/2pi), constant amplitude (-w T/2pi).  Jackknife over curves.
Output: onepoint_all.npz, s2_out.txt."""
import json, math, sys, glob, os
import numpy as np
sys.path.insert(0, '.')
import elib as L

T1, T2 = 20.0, 400.0; T = T2 - T1
C = np.load('curves.npz')
ids = [str(s) for s in C['ids']]; sets = [str(s) for s in C['set']]
DEVZ = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'r148', 'dev_zeros.npz'))
FZ = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'r144', 'fresh_zeros.npz'))
def zeros(i):
    return np.sort(DEVZ[ids[i]]) if sets[i] == 'dev' else np.sort(np.array(FZ[ids[i]], float))
G = np.zeros((len(ids), 2), complex); nz = np.zeros(len(ids), int)
for i in range(len(ids)):
    g = zeros(i); b = g[(g > T1) & (g <= T2)]; nz[i] = len(b)
    ph = np.exp(2j * L.theta_E(int(C['N'][i]), b)); x = (b - T1) / T
    G[i] = ph.sum(), (ph * np.sin(math.pi * x) ** 2).sum()
w = C['w'].astype(float); c = C['c']
law = -w * c * T / (2 * math.pi)
np.savez('onepoint_all.npz', ids=np.array(ids), set=np.array(sets), G=G, law=law, nz=nz)

out = open('s2_out.txt', 'w')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')

def fit(A, y, units=None):
    A = np.atleast_2d(A.T).T if A.ndim == 1 else A
    coef = np.linalg.lstsq(A, y, rcond=None)[0]
    n = len(y); reps = []
    for k in range(n):
        m = np.ones(n, bool); m[k] = False
        reps.append(np.linalg.lstsq(A[m], y[m], rcond=None)[0])
    reps = np.array(reps); se = np.sqrt((n - 1) / n * ((reps - reps.mean(0)) ** 2).sum(0))
    r = y - A @ coef
    return coef, se, float(np.sqrt(np.mean(r ** 2)))

def show(name, res):
    c_, s_, rs = res
    P('  %-58s %s  (resid sd %.2f)' % (name, ', '.join('%.4f +- %.4f' % (a, b) for a, b in zip(c_, s_)), rs))

P('one-point sums, band (20, 400], T = 380; law -w c_E T/2pi (c_E = B_E/L(1,Sym^2 E)); Hann law = law/2')
P('zeros per curve in the band: %d .. %d' % (nz.min(), nz.max()))
for nm, m in (('all 158', np.ones(len(ids), bool)), ('development 20', np.array([s == 'dev' for s in sets])), ('fresh 138', np.array([s == 'fresh' for s in sets]))):
    P(nm + ':')
    y = G[m, 0].real; yh = G[m, 1].real; l = law[m]; one = np.ones(m.sum())
    show('slope of Re G (sharp) on the law, no intercept', fit(l[:, None], y))
    show('[intercept, slope] Re G (sharp)', fit(np.column_stack([one, l]), y))
    show('slope of Re G (Hann) on law/2', fit((l / 2)[:, None], yh))
    show('Im G (sharp): mean', fit(one[:, None], G[m, 0].imag))
    show('Im G (Hann): mean', fit(one[:, None], G[m, 1].imag))
    show('no root number: slope on -c T/2pi', fit((-c[m] * T / (2 * math.pi))[:, None], y))
    show('[law, no-root predictor] jointly', fit(np.column_stack([l, -c[m] * T / (2 * math.pi)]), y))
    show('[-w T/2pi (constant amplitude), -w (c-1) T/2pi]', fit(np.column_stack([-w[m] * T / (2 * math.pi), -w[m] * (c[m] - 1) * T / (2 * math.pi)]), y))
    r = y - l; P('  Re G - law: rms %.2f (sharp); Hann: rms %.2f; |law| median %.1f; sqrt(zeros) %.1f' % (
        np.sqrt(np.mean(r ** 2)), np.sqrt(np.mean((yh - l / 2) ** 2)), np.median(np.abs(l)), np.sqrt(np.median(nz[m]))))
    P('  H0 (zero mean): chi2-like sum (Re G)^2 / (resid var) = %.1f on %d curves' % (np.sum(y ** 2) / np.mean(r ** 2), m.sum()))
    P('  per-root-number means of Re G: w=+1: %.2f (law %.2f), w=-1: %.2f (law %.2f)' % (
        y[w[m] == 1].mean(), l[w[m] == 1].mean(), y[w[m] == -1].mean(), l[w[m] == -1].mean()))
out.close()
