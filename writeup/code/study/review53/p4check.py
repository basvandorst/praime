"""Independent check of the fourth-moment density P4 (frozen c4, Topacogullari (5508) as implemented) for real
primitive chi_D: CFKRS recipe for |L|^4 = L(a1)L(a2)L(b1)L(b2), six swap terms, Ramanujan identity with the local
factors at p | q removed, mpmath at 40 digits, shifts eps*(1, 1.37, 1.71, 2.13) -> 0 by polynomial extrapolation."""
import os, json, math, itertools, sys
import mpmath as mp
mp.mp.dps = 40
P = json.load(open(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/moment3/predictions_triple.json'))
SS = {r['key'][0]: r for r in P['singles']}


def primes_of(n):
    n = abs(n); f = []; p = 2
    while p * p <= n:
        if n % p == 0:
            f.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        f.append(n)
    return f


def Rfun(a, b, ps):
    v = mp.mpf(1)
    for ai in a:
        for bj in b:
            v *= mp.zeta(1 + ai + bj)
    v /= mp.zeta(2 + sum(a) + sum(b))
    for p in ps:
        loc = mp.mpf(1)
        for ai in a:
            for bj in b:
                loc *= (1 - mp.mpf(p) ** (-1 - ai - bj))
        loc /= (1 - mp.mpf(p) ** (-2 - sum(a) - sum(b)))
        v *= loc
    return v


def recipe4(D, t, e):
    q = abs(D); ps = primes_of(q); X = mp.mpf(q) * t / (2 * mp.pi)
    a0 = [e * 1, e * 1.37]; b0 = [e * 1.71, e * 2.13]
    tot = mp.mpf(0)
    for r in range(3):
        for S in itertools.combinations(range(2), r):
            for T in itertools.combinations(range(2), r):
                a = list(a0); b = list(b0)
                # swap: alpha_i (i in S) <-> -beta_j (j in T), matched in order
                for i, j in zip(S, T):
                    a[i], b[j] = -b0[j], -a0[i]
                fac = X ** (-(sum(a0[i] for i in S) + sum(b0[j] for j in T)))
                tot += fac * Rfun(a, b, ps)
    return tot


for D in [int(x) for x in sys.argv[1:]]:
    c4 = SS[D]['c4']
    for t in (100.0, 1000.0):
        epss = [mp.mpf('0.004'), mp.mpf('0.002'), mp.mpf('0.001'), mp.mpf('0.0005')]
        vals = [recipe4(D, mp.mpf(t), e) for e in epss]
        # Neville extrapolation to eps = 0
        n = len(epss); Pn = list(vals)
        for k in range(1, n):
            for i in range(n - k):
                Pn[i] = (epss[i + k] * Pn[i] - epss[i] * Pn[i + 1]) / (epss[i + k] - epss[i])
        ext = Pn[0]
        fro = sum(c * math.log(t) ** k for k, c in enumerate(c4))
        print('D %6d t %5.0f recipe (swap sum) %.8f  frozen c4 %.8f  rel %.2e' % (D, t, float(ext), fro, float(ext) / fro - 1), flush=True)
