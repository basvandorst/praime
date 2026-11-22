"""Round 61 (theory, no data): the information that the last k gaps carry about the next one, in the limit of large x,
from a hidden-residue model that is exact at leading order. At each odd prime q the local factor of the singular series
counts the residues left free by the prime positions, which is the law of: a hidden residue a (uniform, a != 0 mod q),
and prime positions i.i.d. uniform on the q - 1 residues other than a. Given the current prime and k before it, a is
uniform on the residues not seen; with m distinct residues seen, the next position has chance 1/(q-1) on each seen
residue and (q-m-1)/((q-1)(q-m)) on each unseen one. m grows by one with chance (q-1-m)/(q-1). So
  I_k = sum_{q >= 3} [H_q(1) - E H_q(m_k)],   I_inf(k -> inf) = sum_q ((q-2)/(q-1)) ln((q-1)/(q-2)) = infinity.
Checks: I_1 must equal round 56's I_inf = 0.1869077 nats and I_2 - I_1 round 58's J_inf = 0.1046134 nats.
usage: r61.py"""
import os, json, math
import numpy as np
import common as C

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')


def Hq(q, m):
    m = np.asarray(m, float)
    a = 1 / (q - 1); c = (q - m - 1) / ((q - 1) * (q - m))
    with np.errstate(divide='ignore', invalid='ignore'):
        t = np.where(c > 0, (q - m) * c * np.log(c), 0.0)
    return -(m * a * math.log(a) + t)


def I_k_q(q, ks):
    """information carried about the next gap by the last k gaps, local at q, for each k in ks (sorted)."""
    M = q - 1                                   # m runs over 1..q-1
    dist = np.zeros(M + 1); dist[1] = 1.0
    H = Hq(q, np.arange(M + 1)); H1 = H[1]
    out = {}; k = 0
    m = np.arange(M + 1); up = np.where(m >= 1, (q - 1 - m) / (q - 1), 0.0)
    for kk in ks:
        while k < kk:
            new = dist * (1 - up); new[1:] += (dist * up)[:-1]; dist = new; k += 1
        out[kk] = float(H1 - np.sum(dist[1:] * H[1:]))
    return out


if __name__ == '__main__':
    ks = [1, 2, 3, 4, 5, 10, 30, 100, 1000, 10000]
    qs = [int(q) for q in C.small_primes(200000)[1:]]
    tot = {k: 0.0 for k in ks}
    for q in qs:
        if q < 3:
            continue
        if q > 20000:
            # m_k is then almost always k+1 (all distinct): use H(1) - H(k+1) with a small-collision correction ignored
            for k in ks:
                if k + 1 < q:
                    tot[k] += float(Hq(q, 1) - Hq(q, k + 1))
            continue
        r = I_k_q(q, ks)
        for k in ks:
            tot[k] += r[k]
    full = {Q: float(sum(((q - 2) / (q - 1)) * math.log((q - 1) / (q - 2)) for q in qs if 3 <= q <= Q)) for Q in (10, 100, 1000, 10 ** 4, 10 ** 5)}
    out = {'I_k_nats': {str(k): v for k, v in tot.items()}, 'I_k_bits': {str(k): v / math.log(2) for k, v in tot.items()},
           'increments_bits': {str(k): (tot[k] - tot[ks[i - 1]]) / math.log(2) if i else tot[k] / math.log(2) for i, k in enumerate(ks)},
           'check_I1_vs_round56': tot[1] - 0.1869077, 'check_J_vs_round58': (tot[2] - tot[1]) - 0.1046134,
           'full_information_up_to_Q_nats': {str(Q): v for Q, v in full.items()},
           'note': 'primes up to 2e5; for q > 2e4 the all-distinct approximation; the sum over q of the k -> inf limit diverges like ln ln Q'}
    for k in ks:
        print('k %6d  I_k %.5f bits' % (k, tot[k] / math.log(2)))
    print('checks:', out['check_I1_vs_round56'], out['check_J_vs_round58'])
    print('full info (k -> inf) up to Q:', {Q: round(v / math.log(2), 4) for Q, v in full.items()}, 'bits')
    json.dump(out, open(os.path.join(OUT, 'r61.json'), 'w'), indent=1)
