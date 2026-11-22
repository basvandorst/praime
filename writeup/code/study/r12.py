"""Round 12 (diagnostic of the open point of section 10): where in h does the variance excess at H = 1e4 come from?

Identity: for intervals I_j = [L + jH, L + (j+1)H), sum_j c_j^2 = sum_j c_j + 2 sum_{h=1}^{H-1} P_emp(h), with
P_emp(h) = number of prime pairs (p, p + h) inside one interval. Hardy-Littlewood (aligned grid, section 4) expects
P_HL(h) = sum_j c_H R(h) N_h / lam_j^2. The variance excess (measured - law) times the mean is then
  2/J sum_h (P_emp(h) - P_HL(h))  -  [(mean c)^2 - (H/lam)^2 terms],
and the cumulative sum over h shows which distances carry it. Window V = [1e9, 2e9), H = 1e4.
"""
import os, json, math
import numpy as np
import common as C

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
C2 = 0.66016181584686957393
H = 10_000
L0, W0 = 10 ** 9, 10 ** 9


def R_and_N():
    a = np.arange(H)
    u = (np.gcd(a, 10) == 1).astype(np.float64)
    Nh = np.rint(np.fft.irfft(np.abs(np.fft.rfft(u, 2 * H)) ** 2, 2 * H)[:H])     # Nh[h] = sum_a u(a) u(a+h)
    R = np.full(H, C2 / (0.6 / 0.64))
    for p in C.small_primes(H)[1:]:
        p = int(p)
        if p != 5:
            R[p::p] *= (p - 1) / (p - 2)
    return R, Nh


def main():
    R, Nh = R_and_N()
    c = 1 / (0.5 * 0.8) ** 2
    Pemp = np.zeros(H); Phl = np.zeros(H)
    sum_c = sum_c2 = sum_e = sum_e2 = 0.0; J = 0
    inv_lam2 = 0.0
    for s in range(100):
        lo = L0 + s * 10 ** 7
        M = C.primes_mask(lo, 10 ** 7).reshape(-1, H).astype(np.float64)       # 1000 intervals
        F = np.fft.rfft(M, 2 * H, axis=1)
        ac = np.fft.irfft(np.abs(F) ** 2, 2 * H, axis=1)[:, :H]                # ac[j, h] = pairs at distance h
        Pemp += np.rint(ac.sum(0))
        cnt = M.sum(1)
        mid = lo + np.arange(M.shape[0]) * H + H / 2
        lam = np.log(mid)
        Phl += c * R * Nh * np.sum(1 / lam ** 2)
        sum_c += cnt.sum(); sum_c2 += (cnt ** 2).sum(); J += len(cnt)
        e = H / lam; sum_e += e.sum(); sum_e2 += (e ** 2).sum()
    Pemp[0] = Phl[0] = 0
    mean_c = sum_c / J
    var_meas = sum_c2 / J - mean_c ** 2
    # the law: E c^2 = E c + 2 sum P_HL / J; its variance uses the law's mean H/lam
    var_law = sum_e / J + 2 * Phl.sum() / J - sum_e2 / J
    d = 2 * (Pemp - Phl) / J
    cum = np.cumsum(d)
    out = {'J': J, 'mean_c': mean_c, 'mean_e': sum_e / J, 'var_meas': var_meas, 'var_law': var_law,
           'ratio_meas': var_meas / (sum_e / J), 'ratio_law': var_law / (sum_e / J),
           'pair_part': float(d.sum()), 'mean_part': float((sum_c / J - sum_e / J) - (mean_c ** 2 - sum_e2 / J)),
           'cum_at': {str(h): float(cum[h]) for h in (2, 6, 30, 100, 300, 1000, 3000, 9999)},
           'rel_pair_excess_at': {str(h): float((Pemp[:h + 1].sum() - Phl[:h + 1].sum()) / Phl[:h + 1].sum())
                                  for h in (2, 6, 30, 100, 300, 1000, 3000, 9999)},
           'twins_emp': float(Pemp[2]), 'twins_hl': float(Phl[2])}
    # density correction: pair counts scale with the density squared, so with the mean count e + delta the law's
    # pairs become P_HL (1 + delta/e)^2. Residual = pair excess beyond that, cumulative in h, relative to sum P_HL.
    delta = mean_c - sum_e / J
    f = (1 + delta / (sum_e / J)) ** 2
    cumE = np.cumsum(Pemp); cumH = np.cumsum(Phl)
    out['delta'] = float(delta); out['rel_delta'] = float(delta / (sum_e / J))
    out['resid_rel_pair_excess_at'] = {str(h): float((cumE[h] - f * cumH[h]) / cumH[h]) for h in (30, 100, 300, 1000, 3000, 9999)}
    out['resid_pair_part'] = float(2 * (cumE[-1] - f * cumH[-1]) / J)
    out['excess_var_blocks_units'] = float(out['resid_pair_part'] + delta - delta ** 2)
    print(json.dumps(out, indent=1))
    json.dump(out, open(os.path.join(OUT, 'r12.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
