"""Round 73 check, own code: reveal time T of p mod 30030 from past gaps, per-q T_q, coupon predictions,
P(gap = 0 mod q). usage: check73.py p1e9.npy"""
import sys, math, json
import numpy as np
import gmpy2
from fractions import Fraction

P = np.load(sys.argv[1])
QS = (3, 5, 7, 11, 13)
KMAX = 200

# spot-check the prime list with gmpy2 on a stretch
x = int(P[1000]); for_chk = [x]
for _ in range(2000):
    x = int(gmpy2.next_prime(x)); for_chk.append(x)
assert np.array_equal(P[1000:1000 + 2001], np.array(for_chk)), 'sieve mismatch'

g = np.diff(P)
out = {'n_primes': int(len(P)), 'n_gaps': int(len(g))}
out['P_gap0'] = {q: float(np.mean(g % q == 0)) for q in QS}
out['P_gap0_se'] = {q: float(math.sqrt(out['P_gap0'][q] * (1 - out['P_gap0'][q]) / len(g))) for q in QS}


def reveal(tg):
    Tq = {}
    for q in QS:
        T = np.full(len(tg), -1)
        seen = np.zeros((len(tg), q), bool); seen[:, 0] = True
        for k in range(1, KMAX + 1):
            seen[np.arange(len(tg)), (P[tg] - P[tg - k]) % q] = True
            T[(seen.sum(1) == q - 1) & (T < 0)] = k
        Tq[q] = T
    return Tq


def brute_T(t):
    """T by direct class-set intersection over 30030 (independent of the CRT argument)."""
    p = int(P[t]); A = [a for a in range(30030) if math.gcd(a, 30030) == 1]
    k = 0
    while len(A) > 1:
        k += 1; D = p - int(P[t - k])
        A = [a for a in A if math.gcd(a - D, 30030) == 1]
    assert A[0] == p % 30030
    return k


rng = np.random.default_rng(20261122)
tg = rng.choice(np.arange(KMAX + 1, len(P) - 1), 5000, replace=False)
Tq = reveal(tg); T = np.max(np.stack([Tq[q] for q in QS]), 0)
assert np.all(T > 0)
for i in range(30):
    assert brute_T(int(tg[i])) == T[i]
out['own_sample'] = {'n': len(tg), 'mean_T': float(T.mean()), 'se_T': float(T.std() / math.sqrt(len(T))), 'sd_T': float(T.std()),
                     **{'mean_T_q%d' % q: float(Tq[q].mean()) for q in QS}}
# the study's targets (seed 72, same choice expression)
rng = np.random.default_rng(72)
tg2 = rng.choice(np.arange(61, len(P) - 1), 20000, replace=False); tg2 = tg2[tg2 > 150]
Tq2 = reveal(tg2); T2 = np.max(np.stack([Tq2[q] for q in QS]), 0)
out['study_targets'] = {'n': len(tg2), 'mean_T': float(T2.mean()), 'sd_T': float(T2.std()), **{'mean_T_q%d' % q: float(Tq2[q].mean()) for q in QS},
                        'cdf': [float(np.mean(T2 <= k)) for k in range(151)]}
# per-q repeat structure for the study's targets: chance the step adds a new residue, by m
# coupon predictions: exact per-q mean (q-1) H_{q-2}, and mean of the max over q from the product of CDFs


def coupon_cdf(q, kmax):
    M = q - 1; dist = np.zeros(M + 1); dist[1] = 1.0; m = np.arange(M + 1)
    up = np.where(m >= 1, (q - 1 - m) / (q - 1), 0.0); F = [dist[M]]
    for k in range(kmax):
        new = dist * (1 - up); new[1:] += (dist * up)[:-1]; dist = new; F.append(dist[M])
    return np.array(F)


F = {q: coupon_cdf(q, 3000) for q in QS}
prod = np.prod(np.stack([F[q] for q in QS]), 0)
out['coupon_mean_T'] = float(np.sum(1 - prod))
out['coupon_mean_T_q'] = {q: float(Fraction(q - 1) * sum(Fraction(1, j) for j in range(1, q - 1))) for q in QS}
out['coupon_cdf_at'] = {k: float(prod[k]) for k in (20, 30, 34, 40, 50, 60, 80)}
out['data_cdf_at'] = {k: out['study_targets']['cdf'][k] for k in (20, 30, 34, 40, 50, 60, 80)}
out['max_cdf_diff_vs_coupon'] = float(np.max(np.abs(np.array(out['study_targets']['cdf']) - prod[:151])))
# independence across q on the study targets
joint = np.array([np.mean(T2 <= k) for k in range(151)]); pr = np.array([np.prod([np.mean(Tq2[q] <= k) for q in QS]) for k in range(151)])
out['max_joint_minus_product'] = float(np.max(np.abs(joint - pr)))
# lag-2 and lag-3 repeat probabilities in the data (for round 75)
for lag in (1, 2, 3):
    d = P[lag:] - P[:-lag]
    out['P_D%d_0' % lag] = {q: float(np.mean(d % q == 0)) for q in QS}
print(json.dumps({k: v for k, v in out.items() if k != 'study_targets'}, indent=1))
print(json.dumps({k: v for k, v in out['study_targets'].items() if k != 'cdf'}))
json.dump(out, open('check73.json', 'w'), indent=1)
