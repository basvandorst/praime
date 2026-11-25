"""Independent check of two Rankin-Selberg cases: Dirichlet coefficients built in Python from ellap (own multiplicative expansion),
passed to PARI as a vector; lfuncheckfeq with eps = +1 and -1; theta check e^{2i theta}L/conj L = eps at t = 5.1, 17.3;
c by an own direct product (exact local factors p <= 5000, first order to 2e5, plain partial products and log-Cesaro)."""
import json, math, sys
import numpy as np, cypari2
sys.path.insert(0, '.')
from r1_offsets import theta
pari = cypari2.Pari(); pari.allocatemem(700 * 10 ** 6); pari.default('realprecision', 38)
C4 = '/home/bas/praime/writeup/code/study/lowzero/universal/blind4/constants4/'
NMAX = 200000
for key in ('RS_11x53_12', 'RS_14x37_19'):
    d = json.load(open(C4 + key + '.json')); s1, s2 = d['spec'][1], d['spec'][2]
    E1 = pari.ellinit(s1); E2 = pari.ellinit(s2); N1 = int(pari.ellglobalred(E1)[0]); N2 = int(pari.ellglobalred(E2)[0])
    a1 = pari.ellan(E1, NMAX); a2 = pari.ellan(E2, NMAX)
    pr = [int(p) for p in pari.primes([2, NMAX])]
    # local polynomial (arithmetic) of L(E1 x E2) at p
    def loc(p):
        a = int(a1[p - 1]); b = int(a2[p - 1])
        if N1 % p == 0: return [1, -a * b, a * a * p]
        if N2 % p == 0: return [1, -a * b, b * b * p]
        return [1, -a * b, p * a * a + p * b * b - 2 * p * p, -a * b * p * p, p ** 4]
    # Dirichlet coefficients up to 30000 (multiplicative), for the functional-equation check
    M = 30000; an = np.zeros(M + 1, dtype=object); an[1] = 1
    pp = {}
    for p in pr:
        if p > M: break
        P = loc(p); lam = [1]
        k = 1
        while p ** k <= M:
            lam.append(-sum(P[i] * lam[k - i] for i in range(1, min(k, len(P) - 1) + 1))); k += 1
        pp[p] = lam
    # sieve-style multiplicative fill
    for n in range(2, M + 1):
        m = n; val = 1
        for p in pr:
            if p * p > m: break
            if m % p == 0:
                e = 0
                while m % p == 0: m //= p; e += 1
                val *= pp[p][e]
        if m > 1: val *= pp[m][1]
        an[n] = val
    vec = pari('[' + ','.join(str(int(x)) for x in an[1:]) + ']')
    res = {}
    for eps in (1, -1):
        L = pari.lfuncreate([vec, 0, [-1, 0, 0, 1], 3, (N1 * N2) ** 2, eps])
        res[eps] = float(pari.lfuncheckfeq(L))
    L = pari.lfuncreate([vec, 0, [-1, 0, 0, 1], 3, (N1 * N2) ** 2, 1])
    fe = []
    for t in (5.1, 17.3):
        v = complex(pari.lfun(L, pari('3/2 + %s*I' % t)))
        fe.append(np.round(np.exp(2j * theta(t, (N1 * N2) ** 2, [0.0, 1.0, 1.0, 2.0])) * v / v.conjugate(), 8))
    # own c
    lg = 0.0; cum = []
    for p in pr:
        if p <= 5000:
            P = [c * p ** (-j) for j, c in enumerate(loc(p))]          # analytic: coefficient of X^j times p^{-j(k-1)/2}, k-1 = 2
            J = 60; lam = [1.0] + [0.0] * J
            for j in range(1, J + 1): lam[j] = -sum(P[i] * lam[j - i] for i in range(1, min(j, len(P) - 1) + 1))
            F = sum(lam[j] * P[j] * p ** (-j) for j in range(len(P)))
            lg += math.log(F / (1 - 1 / p))
        else:
            l = float(a1[p - 1]) * float(a2[p - 1]) / p
            lg += math.log((1 - l * l / p) / (1 - 1 / p))
        cum.append(lg)
    cum = np.array(cum); P = np.array(pr); sel = P >= NMAX // 10; w = 1.0 / P[sel]; corr = 3 / (np.sqrt(P[sel]) * np.log(P[sel]))
    ces = math.exp(np.sum((cum[sel] + corr) * w) / w.sum())
    print('%s: N1 %d N2 %d  lfuncheckfeq eps=+1 %.0f, eps=-1 %.0f ; theta check %s ; own c: partial 1e5 %.4f, 2e5 %.4f, log-Cesaro %.4f ; registered c %.4f (eps %+d)' % (
        key, N1, N2, res[1], res[-1], fe, math.exp(cum[P <= 100000][-1]), math.exp(cum[-1]), ces, d['c'], d['eps'][0]))
