"""Replication on the 45 genus-2 L-functions of the same box that were NOT drawn in round 152 (genus2.json 'primitive' minus
testset4 g2; includes the three timing/proxy curves 1549, 1989, 2981).  For each: zeros to T = 150 with the registered routine
(bz4_lib.zeros, imported read-only from the repo), c by the registered direct rule (exact local factors p <= 300 via lfuneuler,
first order (1 - l^2/p)/(1 - 1/p) to 5e4, log-Cesaro over [P/10, P] with bias 2/(sqrt P ln P)), eps, small-prime data.
Resumable: one npz per curve in rep/.  usage: r9_replicate.py budget_seconds"""
import sys, os, json, math, time, cmath
REPO = '/home/bas/praime/writeup/code/study/lowzero/universal/blind4'
sys.path.insert(0, REPO); sys.path.insert(0, os.path.dirname(REPO))
import numpy as np, cypari2
from bz4_lib import zeros, params
from ulib import to_c, local_poly, local_data, primes_upto
t0 = time.time(); budget = float(sys.argv[1]); par = int(sys.argv[2]) if len(sys.argv) > 2 else -1; stack = int(sys.argv[3]) if len(sys.argv) > 3 else 1100
only = int(sys.argv[4]) if len(sys.argv) > 4 else -1   # retry one index in a fresh process (added for curve 13)
if only >= 0:
    # Smaller demand for the retry: lfungenus2 leaves the root number unset (ldata[6] = 0), so lfuninit solves for it at the
    # precision of the full domain; for curve 13 (N 1403) that step overflows a 600 MB stack. Here the root number is set
    # to lfunrootres(L) (the value the routine reads as eps anyway) before lfuninit. Same L-function, same routine otherwise.
    import bz4_lib
    _build = bz4_lib.build
    def _build_eps(pari_, spec):
        L_ = _build(pari_, spec)
        if spec[0] != 'g2': return L_
        return pari_('my(L = lfungenus2([%s,%s])); L[6] = lfunrootres(L)[3]; L' % (spec[1], spec[2]))
    bz4_lib.build = _build_eps
pari = cypari2.Pari(); pari.allocatemem(stack * 10 ** 6)
G = json.load(open(os.path.join(REPO, 'genus2.json'))); TS = json.load(open(os.path.join(REPO, 'testset4.json')))
drawn = set((g[0], g[1], g[2]) for g in TS['g2'])
todo = [r for r in G['primitive'] if (r[0], r[1], r[2]) not in drawn]
os.makedirs('rep', exist_ok=True)
def c_direct(L, k, P0=300, P=50000, nu=2):
    an = pari.lfunan(L, P); lam = np.array([to_c(an[n - 1]).real for n in range(1, P + 1)]) / np.sqrt(np.arange(1, P + 1))
    pr = primes_upto(P); logs = []
    for p in pr[pr <= P0]:
        p = int(p); F, S, _ = local_data(local_poly(pari, L, p, k), p); logs.append(cmath.log(F) - math.log(1 - 1 / p))
    big = pr[pr > P0]; lb = np.log(1 - lam[big - 1] ** 2 / big) - np.log(1 - 1 / big)
    cum = np.concatenate([np.cumsum(logs), sum(logs) + np.cumsum(lb)])
    sel = pr >= P // 10; w = 1.0 / pr[sel]; corr = nu / (np.sqrt(pr[sel]) * np.log(pr[sel]))
    return complex(np.exp(np.sum((cum[sel] + corr) * w) / w.sum())).real, complex(np.exp(cum[-1])).real, lam
for i, (N, f, h, m2, feq) in enumerate(todo):
    key = 'R_N%d_%d' % (N, i); fn = os.path.join('rep', key + '.npz')
    if os.path.exists(fn) or (par >= 0 and i % 2 != par): continue
    if only >= 0 and i != only: continue
    if only < 0 and stack < 1000 and os.path.exists(fn[:-4] + '.big'): continue
    if time.time() - t0 > budget - 75: print('budget reached', flush=True); break
    spec = ['g2', f, h]
    try:
        g, info = zeros(pari, spec, 150)
    except cypari2.PariError as e:
        open(fn[:-4] + '.big', 'w').write(str(e)); print(key, 'stack overflow at', stack, 'MB: deferred', flush=True)
        pari.allocatemem(stack * 10 ** 6); continue
    L = pari('lfungenus2([%s,%s])' % (f, h)); Q, k, mu = params(pari, L)
    c, c_end, lam = c_direct(L, k)
    info.update(c=c, c_end=c_end, lam_small={str(p): float(lam[p - 1]) for p in (2, 3, 5, 7, 11, 13)}, key=key,
                retry_eps_fixed=(only >= 0))
    np.savez_compressed(fn, g=g, info=json.dumps(info))
    gg = g[(g > 0) & (g <= 150)]
    from scipy.special import loggamma
    th = gg / 2 * np.log(Q / np.pi ** 4) + sum(np.imag(loggamma((0.5 + m + 1j * gg) / 2)) for m in mu)
    Gs = np.exp(2j * th).sum(); eps = info['eps'][0]; pred = -eps * c * 150 / (2 * math.pi)
    print('%-12s N %4d eps %+d m0 %d %s n %d c %+.3f (end %+.3f) G %+7.2f %+7.2fi pred %+7.2f res %+6.2f  [%.0f s]' % (
        key, N, eps, info['m0'], info['status'], len(gg), c, c_end, Gs.real, Gs.imag, pred, Gs.real - pred, time.time() - t0), flush=True)
left = [1 for i, r in enumerate(todo) if not os.path.exists(os.path.join('rep', 'R_N%d_%d.npz' % (r[0], i)))]
print('left', len(left))
