"""Part 1: ascending scan with trial factoring to depth K (factors 2kp+1, k <= K), LL tests per exponent, cost.
Reads kmin_10000_10000000.npy (p < 10^4, exact to k = 10^7) and kmin_100000_1000000.npy (p < 10^5, k <= 10^6).
Cost model (computed on this machine):
  LL(p)   = wall time of one gmpy2 LL test, power-law fit to timings measured below;
  TF(p,K) = rho * min(K, kmin_p) candidates (rho = candidates per k after the mod-8 rule and the sieve <= 3000,
            measured from the runs) times the measured numpy time per candidate (log2 p float-assisted mulmods).
Cost of a scan in LL equivalents = #LL + sum_p TF(p,K)/LL(p)  (every p's TF cost in units of its own LL test)."""
import time, math, numpy as np, gmpy2, sys
sys.path.insert(0, '.')
from llind import ll
MEXP = [2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243]
E = set(MEXP)

# --- LL cost: gmpy2 time per squaring measured with llind.py time (quiet machine), log-log interpolation
TSQ = {1009: 0.35e-6, 4423: 1.68e-6, 9689: 4.81e-6, 21701: 15.13e-6, 44497: 40.19e-6, 86243: 101.81e-6, 99991: 125.30e-6}
_x = np.log(list(TSQ)); _y = np.log(list(TSQ.values()))
tsq = lambda p: math.exp(np.interp(math.log(max(p, 1009)), _x, _y)) * (min(p, 1009) / 1009 if p < 1009 else 1)
LLt = lambda p: max(p - 2, 1) * tsq(p)

# --- TF time per candidate (numpy pow2_mod from tfdeep.py)
src = open('tfdeep.py').read()
ns = {}
exec(src[src.index('def pow2_mod'):src.index('work =')], {'np': np}, ns)
tfc = {}
for p in [5003, 50021]:
    q = (np.uint64(2 * p) * np.arange(10 ** 6, 2 * 10 ** 6, dtype=np.uint64) + np.uint64(1))
    t = time.time(); ns['pow2_mod'](p, q); tfc[p] = (time.time() - t) / len(q)
print('TF seconds per candidate', {p: f'{v:.2e}' for p, v in tfc.items()})
tfcand = lambda p: tfc[5003] + (tfc[50021] - tfc[5003]) * (math.log2(p) - math.log2(5003)) / (math.log2(50021) - math.log2(5003))

def table(X, fname, Kfile, Ks, rho):
    P, kmin = np.load(fname)
    P = [int(p) for p in P]; kmin = [int(k) for k in kmin]
    # validation: every recorded factor divides M_p and is a proper divisor
    for p, k in zip(P, kmin):
        if k: q = 2 * k * p + 1; assert pow(2, p, q) == 1 and q < (1 << p) - 1, (p, k)
    nexp = sum(1 for m in MEXP if m < X)
    print(f'\nprimes < {X}: {len(P) + 1} (incl. 2), exponents {nexp}; table: K, LL tests (= survivors incl. p = 2), rate, LL per exponent, CPU s: TF + LL = total (model)')
    rows = []
    for K in Ks:
        surv = [2] + [p for p, k in zip(P, kmin) if k == 0 or k > K]
        assert all(m in surv for m in MEXP if m < X)
        tf_eq = sum(rho * min(K, k if k else K) * tfcand(p) / LLt(p) for p, k in zip(P, kmin))
        tf_s = sum(rho * min(K, k if k else K) * tfcand(p) for p, k in zip(P, kmin))
        ll_s = sum(LLt(p) for p in surv)
        hits = sum(1 for s in surv if s in E)
        rows.append((K, len(surv), hits))
        print(f'K={K:>9}  LL {len(surv):>5}  hits {hits}  rate {hits / len(surv):6.2%}  LL/exp {len(surv) / hits:7.1f}  CPU {tf_s:9.1f} + {ll_s:8.1f} = {tf_s + ll_s:9.1f} s')
    return rows

# rho: candidates per unit k that reach the powering step, from the survivors over the whole extended range
for X, f0, f1, K0, K1 in [(10000, 'kmin_10000_1000000.npy', 'tfwork_10000_1000000.npy', 10 ** 5, 10 ** 6),
                          (100000, 'kmin_100000_1000000.npy', 'tfwork_100000_1000000.npy', 10 ** 4, 10 ** 6)]:
    P, km = np.load(f0); w = np.load(f1)
    m = km == 0
    rho = float(w[m].sum() / (m.sum() * (K1 - K0)))
    print(f'X={X}: rho = {rho:.4f} candidates per k')
    if X == 10000:
        table(X, 'kmin_10000_10000000.npy', 10 ** 7, [0, 1, 10, 100, 1000, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7], rho)
    else:
        table(X, 'kmin_100000_1000000.npy', 10 ** 6, [0, 1, 10, 100, 1000, 10 ** 4, 10 ** 5, 10 ** 6], rho)
