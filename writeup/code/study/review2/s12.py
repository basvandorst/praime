"""Round 12 check, per block: variance of counts in aligned intervals of H = 1e4 in [L0, L0 + 1e9), and the
decomposition sum c^2 = sum c + 2 sum_h Pemp(h), with block-level noise for the split at h = 3000."""
import math, sys, json, time
import numpy as np
from sieve import small_primes, mask

H = 10 ** 4
L0 = int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 9
PR = small_primes(10 ** 7)


def law_pairs():
    """Pair sum per interval times lambda^2: sum_{a<b} 6.25 R(b - a) u(a) u(b); returns per-h array."""
    a = np.arange(H)
    u = (np.gcd(a, 10) == 1).astype(float)
    Nh = np.rint(np.correlate(u, u, 'full')[H - 1:])          # Nh[h] = sum_a u(a) u(a + h)
    p = PR[(PR != 2) & (PR != 5)].astype(float)
    Cp = math.exp(np.sum(np.log1p(-1 / (p - 1) ** 2)))
    R = np.full(H, Cp)
    for q in PR[PR < H]:
        q = int(q)
        if q in (2, 5):
            continue
        R[q::q] *= (q - 1) / (q - 2)
    per_h = 6.25 * R * Nh
    per_h[0] = 0
    return per_h


def main():
    ph = law_pairs()
    S = ph.sum()
    blocks = []
    t0 = time.time()
    for b in range(10):
        Pemp = np.zeros(H); cs = []; mids = []; inv2 = 0.0
        for s in range(10):
            lo = L0 + b * 10 ** 8 + s * 10 ** 7
            M = mask(lo, 10 ** 7).reshape(-1, H)
            for r0 in range(0, M.shape[0], 250):
                X = M[r0:r0 + 250].astype(np.float32)
                F = np.fft.rfft(X, 2 * H, axis=1)
                ac = np.fft.irfft((F * F.conj()).real, 2 * H, axis=1)[:, :H]
                Pemp += np.rint(ac).sum(0)
                del F, ac
            cs.append(M.sum(1)); mids.append(lo + np.arange(M.shape[0]) * H + H / 2)
        c = np.concatenate(cs).astype(float); mid = np.concatenate(mids); lam = np.log(mid); e = H / lam
        Pemp[0] = 0
        Phl = ph * np.sum(1 / lam ** 2)
        J = len(c)
        # identity check: sum c^2 = sum c + 2 sum Pemp
        ident = float((c ** 2).sum() - c.sum() - 2 * Pemp.sum())
        meas = float(np.var(c - e, ddof=1) / e.mean())
        law = float(np.mean(1 + (2 * S - H * H) / (H * lam)))
        delta = c.mean() - e.mean(); fdc = (1 + delta / e.mean()) ** 2
        cE = np.cumsum(Pemp); cH = np.cumsum(Phl)
        res = {str(h): float((cE[h] - fdc * cH[h]) / cH[h]) for h in (100, 1000, 3000, 9999)}
        far = float(((cE[-1] - cE[3000]) - fdc * (cH[-1] - cH[3000])) / (cH[-1] - cH[3000]))
        blocks.append({'meas': meas, 'law': law, 'delta': float(delta), 'resid': res, 'far_3000_9999': far,
                       'identity_residual': ident, 'mean_e': float(e.mean())})
        print(b, round(time.time() - t0), round(meas, 4), round(law, 4), res, round(far, 7), ident, flush=True)
    out = {'blocks': blocks}
    def t(v):
        v = np.array(v); return float(v.mean()), float(v.std(ddof=1) / math.sqrt(len(v))), float(v.mean() / (v.std(ddof=1) / math.sqrt(len(v))))
    out['var_meas_minus_law'] = t([x['meas'] - x['law'] for x in blocks])
    out['mean_meas'] = float(np.mean([x['meas'] for x in blocks])); out['mean_law'] = float(np.mean([x['law'] for x in blocks]))
    for h in ('100', '1000', '3000', '9999'):
        out['resid_' + h] = t([x['resid'][h] for x in blocks])
    out['far'] = t([x['far_3000_9999'] for x in blocks])
    print(json.dumps({k: v for k, v in out.items() if k != 'blocks'}, indent=1))
    json.dump(out, open('/home/bas/scratch/study-review2/s12_%d.json' % L0, 'w'), indent=1)


if __name__ == '__main__':
    main()
