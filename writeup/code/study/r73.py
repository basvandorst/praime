"""Round 73 (preregister/round73.txt): the number T of past gaps needed to determine p mod 30030 from the gaps alone,
against the product of coupon-collector completion probabilities over q = 3, 5, 7, 11, 13. usage: r73.py [diag]
(diag, added after the run: per-q completion times on the same targets, to separate uniformity per q from independence across q)"""
import os, sys, json, math
import numpy as np
import common as C
from r72 import CLS, COP, Q

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
KMAX = 150


def predicted_cdf(kmax):
    cdf = np.ones(kmax + 1)
    for q in (3, 5, 7, 11, 13):
        M = q - 1; dist = np.zeros(M + 1); dist[1] = 1.0
        m = np.arange(M + 1); up = np.where(m >= 1, (q - 1 - m) / (q - 1), 0.0)
        done = np.zeros(kmax + 1); done[0] = dist[M]
        for k in range(1, kmax + 1):
            new = dist * (1 - up); new[1:] += (dist * up)[:-1]; dist = new; done[k] = dist[M]
        cdf *= done
    return cdf


def diag(P, tg):
    out = {}
    Tq = {}
    for q in (3, 5, 7, 11, 13):
        T = np.full(len(tg), -1)
        for b0 in range(0, len(tg), 2000):
            t = tg[b0:b0 + 2000]; seen = np.zeros((len(t), q), bool); seen[:, 0] = True; f = np.full(len(t), -1)
            for k in range(1, KMAX + 1):
                seen[np.arange(len(t)), (P[t] - P[t - k]) % q] = True
                f[(seen.sum(1) == q - 1) & (f < 0)] = k
            T[b0:b0 + len(t)] = f
        Tq[q] = T
        out['mean_T_q%d' % q] = float(np.mean(T))
    prod = np.ones(KMAX + 1); joint = np.zeros(KMAX + 1)
    Tall = np.max(np.stack([Tq[q] for q in Tq]), 0)
    for k in range(KMAX + 1):
        prod[k] = np.prod([np.mean(Tq[q] <= k) for q in Tq]); joint[k] = np.mean(Tall <= k)
    out['max_diff_joint_vs_product_of_empirical'] = float(np.max(np.abs(joint - prod)))
    pred = predicted_cdf(KMAX)
    for q in (3, 5, 7, 11, 13):
        M = q - 1; dist = np.zeros(M + 1); dist[1] = 1.0; m = np.arange(M + 1); up = np.where(m >= 1, (q - 1 - m) / (q - 1), 0.0); mean = 0.0
        for k in range(0, 2000):
            mean += 1 - dist[M]; new = dist * (1 - up); new[1:] += (dist * up)[:-1]; dist = new
        out['mean_T_q%d_coupon' % q] = mean
    return out


if __name__ == '__main__':
    P = np.nonzero(C.primes_mask(10 ** 9, 10 ** 8))[0].astype(np.int64) + 10 ** 9
    rng = np.random.default_rng(72)
    tg = rng.choice(np.arange(61, len(P) - 1), 20000, replace=False)          # the same targets as round 72
    tg = tg[tg > KMAX]
    if len(sys.argv) > 1 and sys.argv[1] == 'diag':
        res = json.load(open(os.path.join(OUT, 'r73.json'))); res['diag'] = diag(P, tg); print(res['diag'])
        json.dump(res, open(os.path.join(OUT, 'r73.json'), 'w'), indent=1); sys.exit()
    T = np.full(len(tg), -1)
    for b0 in range(0, len(tg), 1000):
        t = tg[b0:b0 + 1000]
        mask = np.ones((len(t), len(CLS)), bool); found = np.full(len(t), -1)
        for k in range(1, KMAX + 1):
            D = P[t] - P[t - k]
            mask &= COP[(CLS[None, :] - D[:, None]) % Q]
            one = (mask.sum(1) == 1) & (found < 0)
            found[one] = k
        T[b0:b0 + len(t)] = found
        # sanity: the surviving class is the true one
        assert np.all(CLS[mask.argmax(1)] == P[t] % Q)
    pred = predicted_cdf(KMAX)
    emp = np.array([np.mean((T >= 0) & (T <= k)) for k in range(KMAX + 1)])
    dmax = float(np.max(np.abs(emp - pred))); kat = int(np.argmax(np.abs(emp - pred)))
    mean_pred = float(np.sum(1 - pred[:-1])); mean_emp = float(np.mean(T[T >= 0]))
    out = {'targets': len(tg), 'undetermined_by_150': int(np.sum(T < 0)), 'max_cdf_diff': dmax, 'at_k': kat,
           'mean_T_data': mean_emp, 'mean_T_predicted': mean_pred, 'median_T_data': float(np.median(T[T >= 0])),
           'cdf_data': emp.tolist(), 'cdf_predicted': pred.tolist()}
    for k in (10, 20, 30, 40, 50, 60, 80, 100):
        print('k %3d  data %.4f  predicted %.4f' % (k, emp[k], pred[k]))
    print({a: (round(b, 4) if isinstance(b, float) else b) for a, b in out.items() if not a.startswith('cdf')})
    json.dump(out, open(os.path.join(OUT, 'r73.json'), 'w'), indent=1)
