"""Step 1: curve data for the 20 development and 138 fresh curves (PARI): N, w, bad primes, a_p (p <= 1e5), L(1,Sym^2 E);
the one-point constant c_E = B_E/L(1,Sym^2 E) and its direct partial-product check; the identity h_table = h_lines
(four-fold sum against the lines form) on random prime/value combinations.  Output: curves.npz, s1_out.txt."""
import os, json, math, time, sys
import numpy as np, cypari2
sys.path.insert(0, '.')
import elib as L

D2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'deg2') + '/'
pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6)
dev = json.load(open(D2 + 'dev_set.json'))['curves']; fresh = json.load(open(D2 + 'fresh_set.json'))['curves']
rows = []; t0 = time.time()
for setname, cs in (('dev', dev), ('fresh', fresh)):
    for c in cs:
        I = L.curve_info(pari, c['ainv'], 10 ** 5)
        assert I['N'] == c['N'] and I['w'] == c['w'] and sorted(I['bad']) == sorted(c['bad']), c['id']
        for p, a in c['ap'].items():
            k = int(np.searchsorted(I['primes'], int(p))); assert I['ap'][k] == a
        cE, B, direct = L.c_E(I)
        cE3, _, direct3 = L.c_E(I, 10 ** 3)
        rows.append(dict(set=setname, id=c['id'], N=I['N'], w=I['w'], bad=I['bad'], LS1=I['LS1'], B=B, c=cE,
                         direct=direct, direct1e3=direct3, ap=I['ap'], primes=I['primes']))
print('curves done %.0fs' % (time.time() - t0))
pr = rows[0]['primes']
np.savez('curves.npz', ids=np.array([r['id'] for r in rows]), set=np.array([r['set'] for r in rows]),
         N=np.array([r['N'] for r in rows]), w=np.array([r['w'] for r in rows]), LS1=np.array([r['LS1'] for r in rows]),
         B=np.array([r['B'] for r in rows]), c=np.array([r['c'] for r in rows]),
         direct=np.array([r['direct'] for r in rows]), primes=pr, ap=np.array([r['ap'] for r in rows], dtype=np.int16),
         bad=np.array([json.dumps(r['bad']) for r in rows]))
out = open('s1_out.txt', 'w')
def pr_(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
c = np.array([r['c'] for r in rows]); d = np.array([r['direct'] for r in rows]); LS = np.array([r['LS1'] for r in rows])
B = np.array([r['B'] for r in rows]); d3 = np.array([r['direct1e3'] for r in rows])
pr_('158 curves (20 dev + 138 fresh): root numbers, bad primes and a_p (p <= 13) agree with the set files')
pr_('L(1,Sym^2 E): min %.4f median %.4f max %.4f' % (LS.min(), np.median(LS), LS.max()))
pr_('B_E: min %.4f median %.4f max %.4f' % (B.min(), np.median(B), B.max()))
pr_('c_E = B/L(1,Sym2): min %.4f median %.4f mean %.4f max %.4f' % (c.min(), np.median(c), c.mean(), c.max()))
pr_('direct partial product to 1e5 / c_E: median %.4f, range [%.4f, %.4f], rms log %.4f' % (np.median(d / c), (d / c).min(), (d / c).max(), np.sqrt(np.mean(np.log(d / c) ** 2))))
pr_('direct partial product to 1e3 / c_E: median %.4f, rms log %.4f' % (np.median(d3 / c), np.sqrt(np.mean(np.log(d3 / c) ** 2))))
w = np.array([r['w'] for r in rows]); pr_('root numbers: %d of +1, %d of -1' % ((w == 1).sum(), (w == -1).sum()))
# identity h_table (four-fold sum) = h_lines (lines form)
rng = np.random.default_rng(1); u = np.linspace(0, 6, 7); x = 1j * u; md = 0
for p in (2, 3, 5, 7, 11, 13, 101, 997):
    amax = int(2 * math.sqrt(p))
    a = np.array(list(range(-amax, amax + 1)) + [1, -1]); b = np.array([False] * (2 * amax + 1) + [True, True])
    T = L.h_table(p, a, b, a, b, x); H = L.h_lines(p, a, b, a, b, x)
    md = max(md, float(np.abs(T - H).max()))
pr_('h_table (four-fold sum) vs h_lines (sum over v of S1 S2 p^{xv}): max |diff| = %.1e over p <= 997, all a_p values, bad and good' % md)
# first-order check of h_p for large p
p = 2003; a = np.arange(-89, 90); b = np.zeros(len(a), bool); H = L.h_lines(p, a, b, a, b, x)
lam = a / math.sqrt(p); first = 1 + np.outer(lam, lam)[:, :, None] * (2 * np.cos(u * math.log(p)) / p)[None, None, :]
pr_('p = 2003: max |h_p - (1 + 2 lam1 lam2 cos(u ln p)/p)| = %.2e (O(p^-2) = %.1e)' % (np.abs(H - first).max(), 1 / p ** 2))
out.close()
