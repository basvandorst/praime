# Round 160: the registered readings for the prime family (the formulas of summarize_data.py, which needs a squarefree chunk too).
import numpy as np, sys
z = np.load(sys.argv[1], allow_pickle=True)
T, dim, T2, cnt = z['T'], z['dim'], z['T2'], z['cnt']
mean = T / dim; err = np.sqrt(np.maximum(T2 / dim - mean ** 2, 0) / np.maximum(cnt, 1))
mp, mz = z['prime'] / dim, z['zub'] / dim
sel = err > 0
def chi(m, k): return float(np.sum(((mean - m) / err)[k] ** 2))
cp, cz = chi(mp, sel), chi(mz, sel); rms = float(np.sqrt(np.mean((mean - mp) ** 2))); rz = float(np.sqrt(np.mean((mean - mz) ** 2)))
nb = int(sel.sum()); h = np.arange(len(mean)) < 24
print('levels %d, bins %d' % (int(z['nlev']), nb))
print('chi2 against M_prime %.1f, against M_zub %.1f (ratio %.1f); rms(data - M_prime) %.3f, rms(data - M_zub) %.3f' % (cp, cz, cz / cp, rms, rz))
print('first 24 bins (y <= 6): %.1f / %.1f ; last 24 bins: %.1f / %.1f' % (chi(mp, sel & h), chi(mz, sel & h), chi(mp, sel & ~h), chi(mz, sel & ~h)))
print('D1 chi2(M_prime) <= 84.0: %s ; D2 chi2(M_zub) >= 5 chi2(M_prime): %s ; D3 rms <= 0.10: %s' % (
    'PASS' if cp <= 84.0 else 'FAIL', 'PASS' if cz >= 5 * cp else 'FAIL', 'PASS' if rms <= 0.10 else 'FAIL'))
