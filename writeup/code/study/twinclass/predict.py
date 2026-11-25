"""Frozen-law predictions for one window and modulus (twins by default): the Hardy-Littlewood spread of the pair counts
across the admissible classes (pX), of the linear-response part kappa*Z (pZ), of the remainder R = X - kappa*Z (pR),
their covariance, the regression coefficient beta_HL, plus the two corrections of decomp.py (mean cross-class covariance
of the linear part from round 135's closed form; Chebyshev bias of the prime squares).  Only fourpoint.py is used."""
import math
import numpy as np
from scipy.special import expi
import fourpoint as F

def _li(x): return float(expi(math.log(x)))

def _lam(m):
    for p in F.PR[:30].tolist():
        if m % p == 0:
            while m % p == 0: m //= p
            return math.log(p) if m == 1 else 0.0
    return 0.0

def _nsqrt(c, Q):
    out = 1
    for p in F.qfactors(Q):
        if p == 2: continue
        if c % p == 0: return 0
        out *= 2 if pow(c % p, (p - 1) // 2, p) == 1 else 0
    return out

def _acorr(S, Q):
    ind = np.zeros(Q); ind[S] = 1
    f = np.fft.rfft(ind); return np.rint(np.fft.irfft(f * np.conj(f), Q))

def predict(L0, W, Q, h=2):
    L0 = int(L0); W = int(W)
    ph = F.phi(Q); units = [c for c in range(Q) if math.gcd(c, Q) == 1]
    lam = math.log(L0 + W / 2)
    vY, mY = F.prime_class_v(Q, L0, W)
    covd = np.array([vY if d == 0 else -_lam(Q // math.gcd(d, Q)) / (F.phi(Q // math.gcd(d, Q)) * lam) for d in range(Q)])
    bias = np.array([-(_nsqrt(c, Q) - 1) * (_li(math.sqrt(L0 + W)) - _li(math.sqrt(L0))) / (2 * ph) for c in units])
    adm = F.admissible(Q, h); K = len(adm)
    C = F.constants(Q, h); kap = C['S'] * ph / C['A']
    D = F.decompose(Q, h, L0, W)
    ca = _acorr(adm, Q); ca[0] = 0
    tot = 2 * covd + np.roll(covd, -h) + np.roll(covd, h)
    corr = -D['kbar'] * float(np.sum(ca * tot) / ca.sum())
    uidx = {c: i for i, c in enumerate(units)}
    zb = kap / math.log(L0 + W) * np.array([bias[uidx[b]] + bias[uidx[(b + h) % Q]] for b in adm])
    bpart = float(np.var(zb, ddof=1)) / D['mX']
    pX = D['vX'] + corr + bpart; pZ = D['vZ'] + corr + bpart; pR = D['vR']; pC = D['cXZ'] + corr + bpart
    return {'L0': L0, 'W': W, 'Q': Q, 'h': h, 'K': K, 'mX_per_class': D['mX'], 'pX': pX, 'pZ': pZ, 'pR': pR, 'pcov': pC,
            'beta_HL': pC / pZ, 'se_beta_nominal': math.sqrt(pR / ((K - 1) * pZ)), 'four_point_only': D['vX'],
            'vY_HL': vY, 'cross_corr': corr, 'bias_part': bpart, 'kbar': D['kbar']}


# closed form for windows whose exact mod-210 sum is too long to run here (reported readings only)
ASYM = {210: (0.013535, -0.1725, -0.1731), 2310: (0.011066, -0.1669, 0.0616), 30030: (0.009363, -0.1561, 0.1367)}

def predict_closed(L0, W, Q, h=2):
    """v_X from F(J) = -a ln^2 J + b ln J + c (a derived, b c from the exact HL partial sums, asym_fit.json) by summation
    by parts against the exact profile autocorrelation; P_lin = 2 kbar v_Y (round 135 law for one class); P_R = v_X - P_lin.
    Accuracy on the seen decades: v_X closed - exact = +0.005..+0.008 (210), P_lin +10..15 %."""
    assert h == 2
    L0 = int(L0); W = int(W)
    C = F.constants(Q, 2); a, b, c = ASYM[Q]
    N = W // Q; wf = F.profile_w(L0, Q, 2, C['c'], N)
    J = np.unique(np.concatenate([np.arange(1, 2000), np.logspace(np.log10(2000), np.log10(N - 1), 20000).astype(np.int64)]))
    FJ = -a * np.log(J) ** 2 + b * np.log(J) + c
    wJ = wf(J.astype(float)); wnext = np.append(wJ[1:], 0.0)
    off = 2 * float(np.sum(FJ * (wJ - wnext)))
    m = L0 + Q / 2 + Q * np.linspace(0, N - 1, 200001)
    r2 = C['c'] / (np.log(m) * np.log(m + 2)); r1 = (Q / F.phi(Q)) / np.log(m)
    mX = float(r2.mean() * N); vX = (float((r2 * (1 - r2)).mean() * N) + off) / mX
    kbar = float(np.sum(r2 * r2 / r1) / np.sum(r2))
    lam = math.log(L0 + W / 2); SQ = sum(math.log(p) / (p - 1) for p in F.qfactors(Q))
    vY = 1 - (math.log(W / Q) + SQ - (1 - 0.5772156649015329 - math.log(2 * math.pi))) / lam
    pZ = 2 * kbar * vY
    return {'L0': L0, 'W': W, 'Q': Q, 'h': 2, 'K': len(F.admissible(Q, 2)), 'mX_per_class': mX, 'pX': vX, 'pZ': pZ,
            'pR': vX - pZ, 'closed_form': True, 'kbar': kbar}
