"""Round 60 (preregister/round60.txt): the entropy of the gap between consecutive twin primes against
log2(e/(6 rho_tw)) - K_tw/ln 2 and against round 59's conditioned law. usage: r60.py"""
import os, json, math
import numpy as np
import common as C

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
PR = C.small_primes(20_000_000).astype(np.float64)
P3 = PR[PR > 2]; C2 = math.exp(float(np.sum(np.log1p(-1 / (P3 - 1) ** 2))))


def K_tw_with_tail():
    ps = PR[PR >= 5]
    f2 = 1 - 2 / ps; f3 = 1 - 3 / ps; f4 = 1 - 4 / ps
    mean = (f2 + 2 * f3 + (ps - 3) * f4) / ps
    w2, w3, w4 = f2 / mean, f3 / mean, f4 / mean
    terms = (w2 * np.log(w2) + 2 * w3 * np.log(w3) + (ps - 3) * w4 * np.log(w4)) / ps
    return float(terms.sum()), float(terms[-1] * ps[-1] ** 2)   # terms fall like about 3/p^3; beyond 2e7 negligible


if __name__ == '__main__':
    K, c = K_tw_with_tail()
    r = json.load(open(os.path.join(OUT, 'r59.json')))
    data = {int(k): v for k, v in r['data'].items()}; N = sum(data.values())
    cnt = np.array(list(data.values()), float); pp = cnt / N
    H_data = float(-np.sum(pp * np.log2(pp)) + (len(cnt) - 1) / (2 * N * math.log(2)))
    lam = r['compare']['mean_lambda']; rho = 2 * C2 / lam ** 2
    H_asym = (math.log(math.e / (6 * rho)) - K) / math.log(2)
    # round 59's conditioned law, averaged over the window as in r59 compare, up to d = 1200, plus an exponential tail
    LAMS = [math.log(1e9), math.log(1.5e9), math.log(2e9)]
    u = np.linspace(LAMS[0], LAMS[-1], 401); wt = np.exp(u) / u ** 2
    law = np.array([float(np.sum(np.interp(u, LAMS, r['law'][str(d)]) * wt) / np.sum(wt)) for d in range(6, 1201, 6)])
    tail_mass = 1 - law.sum(); q = 1 - 6 * rho
    tail = tail_mass * (1 - q) * q ** np.arange(0, 2000)
    Pl = np.concatenate([law, tail]); Pl = Pl[Pl > 0]
    H_law = float(-np.sum(Pl * np.log2(Pl)))
    # discrete leading-order law
    d = np.arange(6, 20001, 6); S = np.ones(len(d))
    for p in C.small_primes(20001):
        p = int(p)
        if p < 5:
            continue
        nu = np.where(d % p == 0, 2, np.where((d % p == 2) | (d % p == p - 2), 3, 4))
        S *= (1 - nu / p) / (1 - 2 / p) / (1 - 1 / p) ** 2
    Pd = S * np.exp(-rho * d); Pd /= Pd.sum()
    H_disc = float(-np.sum(Pd * np.log2(Pd)))
    out = {'K_tw_nats': K, 'K_tw_bits': K / math.log(2), 'tail_coef': c, 'mean_lambda': lam, 'twin_gaps': N,
           'H_data': H_data, 'H_law_round59': H_law, 'H_discrete_leading': H_disc, 'H_asym': H_asym,
           'data_minus_law': H_data - H_law, 'data_minus_asym': H_data - H_asym}
    print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in out.items()})
    json.dump(out, open(os.path.join(OUT, 'r60.json'), 'w'), indent=1)
