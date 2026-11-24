"""Validation of the class-count covariance law (derivation.md section 4) on primes in a window that is none of the
registered or proposed test windows: [5e11, 5e11 + 1e8). Windows of H = 1e5 (1000 of them) and H = 1e6 (100).
For each q: per-class Var/mean across windows (law: v_q), total Var/mean (law: v_tot, Montgomery-Soundararajan),
and the equidistribution statistic X = sum_a (N_a - N/phi)^2/(N/phi) per window, mean X/(phi-1) (law: v_q + K_q/((phi-1) lam)).
Errors: from the spread of the per-window values (windows treated as independent; adjacent windows are weakly
anti-correlated, which makes the quoted errors slightly conservative for means of positive quantities).
Also a check of the sieve against sympy on a short range."""
import math, json, sys, time
import numpy as np
from sympy import isprime, totient, factorint
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from sieve import primes_in

EG = 0.5772156649015329
B = 1 - EG - math.log(2 * math.pi)
X0 = int(5e11); SPAN = 10 ** 8
QS = (3, 4, 5, 7, 10, 12, 30, 210)


def law(q, H, lam):
    phi = int(totient(q)); Sq = sum(math.log(p) / (p - 1) for p in factorint(q))
    v = 1 - (math.log(H / q) + Sq - B) / lam
    K = math.log(q) - Sq
    vtot = 1 - (math.log(H) - B) / lam
    return {'v_class': v, 'contrast_per_dof': v + K / ((phi - 1) * lam) if phi > 1 else float('nan'), 'v_total': vtot}


def main():
    # sieve check
    lo = X0; hi = X0 + 20000
    P = np.concatenate(list(primes_in(lo, hi)))
    ref = [n for n in range(lo, hi) if isprime(n)]
    assert list(P) == ref, 'sieve mismatch'
    print('sieve check ok: %d primes in [%d, %d)' % (len(ref), lo, hi))
    t0 = time.time()
    Hs = (10 ** 5, 10 ** 6)
    cnt = {H: {q: np.zeros((SPAN // H, q), np.int64) for q in QS} for H in Hs}
    for P in primes_in(X0, X0 + SPAN):
        for H in Hs:
            w = (P - X0) // H
            for q in QS:
                np.add.at(cnt[H][q], (w, P % q), 1)
    print('sieved %.0f s' % (time.time() - t0), flush=True)
    lam = math.log(X0 + SPAN / 2)
    out = {'X0': X0, 'span': SPAN, 'lambda': lam, 'B': B, 'rows': []}
    for H in Hs:
        for q in QS:
            C = cnt[H][q]
            red = [a for a in range(q) if math.gcd(a, q) == 1]
            C = C[:, red].astype(float); phi = len(red); M = C.shape[0]
            tot = C.sum(1)
            vclass = C.var(0, ddof=1) / C.mean(0)                   # per class
            vtot = tot.var(ddof=1) / tot.mean()
            X = ((C - tot[:, None] / phi) ** 2 / (tot[:, None] / phi)).sum(1) / (phi - 1)
            L = law(q, H, lam)
            row = {'H': H, 'q': q, 'windows': M,
                   'v_class_meas': float(vclass.mean()), 'v_class_law': L['v_class'],
                   'v_class_err': float(vclass.mean() * math.sqrt(2 / (M - 1)) / math.sqrt(phi)),
                   'contrast_meas': float(X.mean()), 'contrast_law': L['contrast_per_dof'],
                   'contrast_err': float(X.std(ddof=1) / math.sqrt(M)),
                   'v_total_meas': float(vtot), 'v_total_law': L['v_total'],
                   'v_total_err': float(vtot * math.sqrt(2 / (M - 1)))}
            out['rows'].append(row)
            print('H=%7d q=%3d  class Var/mean %.3f (law %.3f, +-%.3f)  contrast chi2/dof %.3f (law %.3f, +-%.3f)  total %.3f (law %.3f, +-%.3f)'
                  % (H, q, row['v_class_meas'], row['v_class_law'], row['v_class_err'], row['contrast_meas'], row['contrast_law'],
                     row['contrast_err'], row['v_total_meas'], row['v_total_law'], row['v_total_err']), flush=True)
    json.dump(out, open('validate_cov.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
