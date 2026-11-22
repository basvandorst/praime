"""Round 99 (writeup/data/study/preregister/round99.txt): the tube of fact 2048 for an odd primitive Dirichlet character.
G_chi(x) = sum_{n<=x} Lambda(n) chi(n) (x/n - 1) + x L'/L(1,chi) = sum_rho x^rho/(rho(1-rho)) + L'/L(0,chi) - T(x),
T(x) = sum_{n odd} x^(-n)/(n(n+1)); under GRH |G - L'/L(0) + T| <= C sqrt(x), C = sum_rho 1/|rho|^2.
usage: chartube.py [X]"""
import os, sys, json, math
import numpy as np
import mpmath as mp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'study'))
import common as C

mp.mp.dps = 30
X = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10 ** 9


def Lfun(s, q, chi):
    return mp.power(q, -s) * mp.fsum(chi[a] * mp.zeta(s, mp.mpf(a) / q) for a in range(1, q))


def consts(q, chi):
    dl1 = mp.dirichlet(1, chi, 1) / mp.dirichlet(1, chi)
    dl0 = mp.dirichlet(0, chi, 1) / mp.dirichlet(0, chi)
    Ca = dl1 - dl0 + mp.log(2)
    Cb = 2 * mp.re(dl1) + mp.log(mp.mpf(q) / mp.pi) + mp.digamma(1)      # a = 1: psi((1+a)/2) = psi(1)
    return float(dl1), float(dl0), float(Ca), float(Cb)


def scan(q, chi, dl1, c0):
    """max |G(x) - c0|/sqrt(x) over [100, X] (T neglected there), G piecewise linear between prime powers."""
    A = 0.0; Bs = 0.0                      # A = sum a_n/n, B = sum a_n over n <= current breakpoint
    best = (0.0, None); seg = 10 ** 7
    pw = []                                # prime powers p^k, k >= 2
    for p in C.small_primes(int(math.isqrt(X)) + 2).tolist():
        v = p * p
        while v <= X:
            pw.append(v); v *= p
    pw = np.array(sorted(pw), dtype=np.int64)
    for s0 in range(0, X, seg):
        n = min(seg, X - s0)
        P = np.nonzero(C.primes_mask(max(s0, 2), n - (2 if s0 == 0 else 0)))[0].astype(np.int64) + max(s0, 2)
        extra = pw[(pw >= s0) & (pw < s0 + n)]
        nn = np.concatenate([P, extra]); nn.sort()
        # the prime p of each n and chi(n)
        base = nn.copy(); isp = np.isin(nn, P, assume_unique=False)
        for v in extra.tolist():
            pass
        chiv = np.zeros(len(nn)); logp = np.zeros(len(nn))
        pr = nn[isp]; chiv[isp] = np.array([chi[r] for r in (pr % q).tolist()]) if len(pr) < 0 else np.choose(pr % q, chi)
        logp[isp] = np.log(pr.astype(float))
        for idx in np.nonzero(~isp)[0].tolist():
            v = int(nn[idx]); p = None
            for k in range(39, 1, -1):
                r = round(v ** (1.0 / k))
                for cand in (r - 1, r, r + 1):
                    if cand > 1 and cand ** k == v:
                        p, kk = cand, k; break
                if p:
                    break
            chiv[idx] = chi[p % q] ** kk; logp[idx] = math.log(p)
        a = logp * chiv
        cumA = A + np.cumsum(a / nn); cumB = Bs + np.cumsum(a)
        x = nn.astype(float)
        # value at each breakpoint (G continuous): G(n) = n (A_n + dl1) - B_n, with A_n, B_n including n
        G = x * (cumA + dl1) - cumB
        r = np.abs(G - c0) / np.sqrt(x)
        m = x >= 100
        if m.any():
            i = int(np.argmax(np.where(m, r, 0)));
            if r[i] > best[0]:
                best = (float(r[i]), float(x[i]))
        # interior stationary points of (alpha x + beta)/sqrt(x) on [n_i, n_{i+1}): x* = beta/alpha... for f = (alpha x - beta')/sqrt x
        alpha = cumA + dl1; beta = -(cumB + c0)              # G - c0 = alpha x + beta
        xs = beta / np.where(alpha != 0, alpha, np.nan)
        nxt = np.concatenate([x[1:], [float(s0 + n)]])
        ok = (xs > x) & (xs < nxt) & (xs >= 100)
        if ok.any():
            rv = np.abs(alpha[ok] * xs[ok] + beta[ok]) / np.sqrt(xs[ok])
            j = int(np.argmax(rv))
            if rv[j] > best[0]:
                best = (float(rv[j]), float(xs[ok][j]))
        A = float(cumA[-1]) if len(cumA) else A; Bs = float(cumB[-1]) if len(cumB) else Bs
    return best


if __name__ == '__main__':
    out = {'X': X}
    for name, q, chi in (('chi4', 4, [0, 1, 0, -1]), ('chi3', 3, [0, 1, -1])):
        dl1, dl0, Ca, Cb = consts(q, chi)
        best = scan(q, chi, dl1, dl0)
        out[name] = {'Lp_over_L_1': dl1, 'Lp_over_L_0': dl0, 'C_from_x1': Ca, 'C_hadamard': Cb, 'C_diff': Ca - Cb,
                     'max_ratio': best[0], 'at_x': best[1], 'ratio_over_C': best[0] / Ca}
        print(name, {k: (round(v, 10) if isinstance(v, float) else v) for k, v in out[name].items()}, flush=True)
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'chartube.json'), 'w'), indent=1)
