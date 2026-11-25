"""Measured spread of pair counts across the admissible classes mod Q, its linear-response part (regression on the
class prime errors of the two members, coefficient kappa of round 155) and the remainder, against the Hardy-Littlewood
predictions of each part.  usage: decomp.py win.npz [h ...]"""
import sys; sys.path.insert(0, '/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/twinclass')
import math, json
import numpy as np
from scipy.special import expi
from scipy.stats import chi2
import fourpoint as F

def li(x): return float(expi(math.log(x)))

def Lam(m):
    for p in F.PR[:30].tolist():
        if m % p == 0:
            while m % p == 0: m //= p
            return math.log(p) if m == 1 else 0.0
    return 0.0

def nsqrt(c, Q):
    out = 1
    for p in F.qfactors(Q):
        if p == 2: continue
        out *= 1 + pow(c % p, (p - 1) // 2, p) if c % p else 1
        if c % p and pow(c % p, (p - 1) // 2, p) == p - 1: out = 0
    return out

def fold(a, Q):
    return np.bincount(np.arange(a.shape[-1]) % Q, weights=a, minlength=Q)

def analyse(path, hs=(2,), Qs=(210, 2310, 30030), sub=None):
    D = np.load(path); L0 = int(D['L0']); W = int(D['W'])
    sl = slice(None) if sub is None else sub
    nsub = D['P'].shape[0]
    if sub is not None:                                # one sub-window
        L0 = L0 + (W * sub) // nsub; W = W // nsub
    P = D['P'][sl].reshape(-1, 30030).sum(0).astype(float); WP = D['WP'][sl].reshape(-1, 30030).sum(0)
    lam = math.log(L0 + W / 2)
    res = {'L0': L0, 'W': W}
    for Q in Qs:
        ph = F.phi(Q); units = [c for c in range(Q) if math.gcd(c, Q) == 1]
        PQ = fold(P, Q); WPQ = fold(WP, Q)
        Yu = PQ[units]
        vY_meas = float(np.var(Yu, ddof=1) / Yu.mean())
        vY_hl, mY = F.prime_class_v(Q, L0, W)
        # cross-class covariance (round 135 closed form) and the Chebyshev bias of the prime squares
        covd = np.array([vY_hl if d == 0 else -Lam(Q // math.gcd(d, Q)) / (F.phi(Q // math.gcd(d, Q)) * lam) for d in range(Q)])
        def acorr(S):                                   # cyclic autocorrelation of the indicator of S: pairs at difference d
            ind = np.zeros(Q); ind[S] = 1
            f = np.fft.rfft(ind); return np.rint(np.fft.irfft(f * np.conj(f), Q))
        cu = acorr(units); cu[0] = 0
        mean_cc = float(np.sum(cu * covd) / cu.sum())
        x1, x2 = L0, L0 + W
        bias = np.array([-(nsqrt(c, Q) - 1) * (li(math.sqrt(x2)) - li(math.sqrt(x1))) / (2 * ph) for c in units])
        vY_pred = vY_hl - mean_cc + float(np.var(bias, ddof=1)) / Yu.mean()
        res['primes_%d' % Q] = {'measured': vY_meas, 'HL_same_class': vY_hl, 'pred_with_cross_and_bias': vY_pred,
                                'bias_part': float(np.var(bias, ddof=1)) / Yu.mean()}
        for h in hs:
            key = 'X%d' % h
            X = fold(D[key][sl].reshape(-1, 30030).sum(0).astype(float), Q)
            adm = F.admissible(Q, h); K = len(adm)
            C = F.constants(Q, h)
            kap = C['S'] * ph / C['A']                   # kappa * ln p
            Xb = X[adm]; m = Xb.mean()
            Zb = kap * (WPQ[adm] + WPQ[[(b + h) % Q for b in adm]])
            Zb = Zb - Zb.mean() + m                      # only deviations across classes matter
            Rb = Xb - Zb
            u2X = float(np.var(Xb, ddof=1) / m); u2Z = float(np.var(Zb, ddof=1) / m); u2R = float(np.var(Rb, ddof=1) / m)
            cXZ = float(np.cov(Xb, Zb)[0, 1] / m)
            Dp = F.decompose(Q, h, L0, W)
            # across-class correction for the linear part: -mean_{b != b'} Cov(kappa Z_b, kappa Z_b')/m_X
            kb = Dp['kbar']
            ca = acorr(adm); ca[0] = 0
            tot = 2 * covd + np.roll(covd, -h) + np.roll(covd, h)     # d -> 2 cov(d) + cov(d + h) + cov(d - h)
            corr = -kb * float(np.sum(ca * tot) / ca.sum())
            zb_bias = kap / math.log(x2) * np.array([bias[units.index(b)] + bias[units.index((b + h) % Q)] for b in adm])
            bpart = float(np.var(zb_bias, ddof=1)) / m
            pX = Dp['vX'] + corr + bpart; pZ = Dp['vZ'] + corr + bpart; pR = Dp['vR']
            res['h%d_%d' % (h, Q)] = {
                'K': K, 'pairs': int(Xb.sum()), 'measured_u2X': u2X, 'pred_u2X': pX, 'ratio_X': u2X / pX,
                'p_lower_X': float(chi2.cdf((K - 1) * u2X / pX, K - 1)),
                'measured_u2Z': u2Z, 'pred_u2Z': pZ, 'measured_u2R': u2R, 'pred_u2R': pR, 'ratio_R': u2R / pR,
                'p_lower_R': float(chi2.cdf((K - 1) * u2R / pR, K - 1)),
                'measured_cov_XZ': cXZ, 'pred_cov_XZ': Dp['cXZ'] + corr + bpart,
                'beta_hat': cXZ / u2Z, 'beta_HL': Dp['beta'], 'se_beta_approx': math.sqrt(pR / ((K - 1) * pZ)),
                'HL': {k: Dp[k] for k in ('vX', 'vZ', 'cXZ', 'vR', 'diagX', 'off4', 'diagZ', 'offlin', 'off_same', 'off_cross', 'kbar')},
                'cross_corr': corr, 'bias_part': bpart, 'four_point_only': Dp['vX']}
    return res

if __name__ == '__main__':
    path = sys.argv[1]; hs = tuple(int(a) for a in sys.argv[2:]) or (2,)
    res = analyse(path, hs)
    for k, v in res.items():
        if isinstance(v, dict):
            print(k, json.dumps({a: (round(b, 4) if isinstance(b, float) else b) for a, b in v.items() if a != 'HL'}))
            if 'HL' in v: print('   HL', {a: round(b, 4) for a, b in v['HL'].items()})
    json.dump(res, open(path.replace('.npz', '_decomp.json'), 'w'), indent=1)
