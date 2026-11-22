"""Round 12, with errors (review 2, point 2): the density-corrected pair residual per block of 1e8 (10 blocks), in the
window [L0, L0 + 1e9) for L0 = 1e9 and 3e9. Per block: delta = mean count - mean H/ln x, the law's pairs scaled by
(1 + delta/e)^2, residual = (pairs - scaled law pairs)/law pairs, for all h < 1e4 and for h <= 3000; mean and standard
error over the blocks. (The identity sum c^2 = sum c + 2 sum P_emp makes the variance excess equal to the total
residual by construction, so only the residual's own scatter is evidence.)"""
import os, json, math
import numpy as np
import common as C
from r12 import R_and_N, H

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
R, Nh = R_and_N()
c = 1 / (0.5 * 0.8) ** 2
res = {}
for L0 in (10 ** 9, 3 * 10 ** 9):
    tot, short = [], []
    for b in range(10):
        Pemp = np.zeros(H); Phl = np.zeros(H); sc = se = 0.0; J = 0
        for s in range(10):
            lo = L0 + b * 10 ** 8 + s * 10 ** 7
            M = C.primes_mask(lo, 10 ** 7).reshape(-1, H).astype(np.float64)
            F = np.fft.rfft(M, 2 * H, axis=1)
            Pemp += np.rint(np.fft.irfft(np.abs(F) ** 2, 2 * H, axis=1)[:, :H].sum(0))
            lam = np.log(lo + np.arange(M.shape[0]) * H + H / 2)
            Phl += c * R * Nh * np.sum(1 / lam ** 2)
            sc += M.sum(); se += np.sum(H / lam); J += M.shape[0]
        Pemp[0] = Phl[0] = 0
        f = (1 + (sc - se) / se) ** 2
        tot.append(float((Pemp.sum() - f * Phl.sum()) / Phl.sum()))
        short.append(float((Pemp[:3001].sum() - f * Phl[:3001].sum()) / Phl[:3001].sum()))
    rr = {}
    for name, v in (('all_h', tot), ('h_le_3000', short)):
        v = np.array(v)
        rr[name] = {'mean': float(v.mean()), 'se': float(v.std(ddof=1) / math.sqrt(10)), 't9': float(v.mean() / (v.std(ddof=1) / math.sqrt(10)))}
    res[str(L0)] = rr
    print(L0, rr, flush=True)
json.dump(res, open(os.path.join(OUT, 'r12b.json'), 'w'), indent=1)
