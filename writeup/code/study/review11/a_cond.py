"""Review 11, script A: independent implementation of the conditioned expansion (rounds 47-48).
(1) exact identity check at small g: full HL inclusion-exclusion vs E_a[ ... Pi_a exact ];
(2) own o2 / o3 / overlap vs r28.json exact sums for g = 74, 90;
(3) shares for g = 150 at lam = ln 1e9, ln 1.5e9, ln 2e9 and the [1e9, 2e9) expected count."""
import math, json, sys, itertools
import numpy as np
import gmpy2

def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]

PR = primes_upto(10_000_000)
QMAX = 401
SMALLQ = [int(q) for q in PR if q < QMAX]
BIGQ = PR[PR >= QMAX].astype(float)
TAILK = {k: float(np.sum(np.log1p(-k / BIGQ) - k * np.log1p(-1 / BIGQ))) for k in range(1, 40)}

def lnS(pts, qmin):
    """ln of prod_{q > qmin} (1 - nu/q)/(1 - 1/q)^k (all q >= QMAX assumed all-distinct). -inf if inadmissible."""
    k = len(pts); tot = TAILK[k]
    for q in SMALLQ:
        if q <= qmin: continue
        nu = len({p % q for p in pts})
        if nu >= q: return -math.inf
        tot += math.log(1 - nu / q) - k * math.log(1 - 1 / q)
    return tot

def identity_check(g, y, lam):
    inter = list(range(2, g, 2))
    # full HL sum (all primes, qmin = 1)
    full = 0.0
    for r in range(len(inter) + 1):
        for B in itertools.combinations(inter, r):
            v = lnS([0, g, *B], 1)
            if v > -math.inf: full += (-1) ** r * math.exp(v) / lam ** (r + 1)
    sm = [q for q in SMALLQ if q <= y]; P = math.prod(sm)
    rho0 = 1 / (lam * math.prod(1 - 1 / q for q in sm))
    cls = [a for a in range(P) if math.gcd(a, P) == 1]
    S2 = math.exp(lnS([0, g], y))
    acc = 0.0
    cache = {}
    for a in cls:
        if math.gcd(a + g, P) != 1: continue
        Ra = [b for b in inter if math.gcd(a + b, P) == 1]
        Pi = 0.0
        for r in range(len(Ra) + 1):
            for B in itertools.combinations(Ra, r):
                key = B
                if key not in cache: cache[key] = lnS([0, g, *B], y)
                v = cache[key]
                if v > -math.inf: Pi += (-1) ** r * rho0 ** r * math.exp(v) / S2
        acc += rho0 * S2 * Pi
    cond = acc / len(cls)
    return full, cond

class Cond:
    def __init__(self, g, y=13):
        self.g = g
        sm = [q for q in SMALLQ if q <= y]; self.P = math.prod(sm)
        self.den = math.prod(1 - 1 / q for q in sm)
        self.pts = list(range(2, g, 2)); n = self.n = len(self.pts)
        l2 = lnS([0, g], y); self.S2 = math.exp(l2)
        self.s = np.array([math.exp(lnS([0, g, b], y) - l2) for b in self.pts])
        self.r = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                self.r[i, j] = self.r[j, i] = math.exp(lnS([0, g, self.pts[i], self.pts[j]], y) - l2) / (self.s[i] * self.s[j]) - 1
        self.K = np.zeros((n, n, n))
        for i, j, k in itertools.combinations(range(n), 3):
            t = math.exp(lnS([0, g, self.pts[i], self.pts[j], self.pts[k]], y) - l2) / (self.s[i] * self.s[j] * self.s[k]) - 1
            kap = t - self.r[i, j] - self.r[i, k] - self.r[j, k]
            for (a, b, c) in itertools.permutations((i, j, k)):
                self.K[a, b, c] = kap
        cop = [a for a in range(self.P) if math.gcd(a, self.P) == 1]
        self.ncls = len(cop)
        self.R = [np.array([i for i, b in enumerate(self.pts) if math.gcd(a + b, self.P) == 1]) for a in cop if math.gcd(a + g, self.P) == 1]

    def share(self, lam):
        rho0 = 1 / (lam * self.den); p = rho0 * self.s; w = p / (1 - p)
        z = np.outer(w, w) * self.r            # pair polymer activity (zero diagonal)
        out = {'o1': 0.0, 'o2': 0.0, 'o3': 0.0, 'o2ov': 0.0, 'o3ov': 0.0}
        for idx in self.R:
            L1 = np.sum(np.log1p(-p[idx]))
            Z = z[np.ix_(idx, idx)]
            pair = Z.sum() / 2
            T = self.K[np.ix_(idx, idx, idx)]; wi = w[idx]
            trip = -np.einsum('ijk,i,j,k->', T, wi, wi, wi) / 6
            # second Mayer term: -1/2 sum over ORDERED pairs of overlapping pair-polymers (incl. identical)
            rows = Z.sum(1)
            same = (Z ** 2).sum() / 2                      # identical polymer, each unordered pair once
            share1 = np.sum(rows ** 2 - (Z ** 2).sum(1))   # ordered distinct pairs sharing one point
            ov = -0.5 * (same + share1)
            for key, L in (('o1', L1), ('o2', L1 + pair), ('o3', L1 + pair + trip), ('o2ov', L1 + pair + ov), ('o3ov', L1 + pair + trip + ov)):
                out[key] += math.exp(L)
        return {k: rho0 * self.S2 * v / self.ncls for k, v in out.items()}

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'identity':
        for g, y, lam in ((30, 5, 17.8), (36, 7, 12.0), (34, 5, 20.0)):
            f, c = identity_check(g, y, lam)
            print('identity g=%d y=%d lam=%.1f full %.12e cond-exact %.12e rel %.2e' % (g, y, lam, f, c, c / f - 1), flush=True)
    elif mode == 'r28':
        ref = json.load(open('/home/bas/praime/writeup/data/study/r28.json'))
        for g in map(int, sys.argv[2].split(',')):
            C = Cond(g)
            for lam in (17.8, 23.0):
                ex = ref[str(g)]['lam%.1f' % lam]['exact']
                sh = C.share(lam)
                print(g, lam, 'exact %.6e' % ex, {k: '%+.5f%%' % (100 * (v / ex - 1)) for k, v in sh.items()}, flush=True)
    elif mode == 'g':
        g = int(sys.argv[2]); C = Cond(g)
        lams = [math.log(1e9), math.log(1.5e9), math.log(2e9)]
        res = {lam: C.share(lam) for lam in lams}
        for lam in lams: print(g, round(lam, 4), {k: '%.6e' % v for k, v in res[lam].items()}, flush=True)
        import mpmath
        u = np.linspace(lams[0], lams[2], 4001)
        ngaps = int(sys.argv[3]) if len(sys.argv) > 3 else 47374752
        norm = ngaps / float(mpmath.li(2e9) - mpmath.li(1e9))
        for key in ('o2', 'o3ov'):
            ls = np.interp(u, lams, [math.log(res[l][key]) for l in lams])
            print(g, key, 'expected count %.2f' % (float(np.trapezoid(np.exp(ls) * np.exp(u) / u, u)) * norm), 'norm %.7f' % norm, flush=True)
