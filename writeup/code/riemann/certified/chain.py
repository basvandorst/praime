"""Self-contained certified chain: every prime from 3 to X proved from zeta zeros, no prime as input.

Stages (end, a, beta, m, c, k): numbers below `end` are decided with kernel support a, order m,
Kaiser-Bessel beta, zeros up to T(x) = c*beta*x/a for the sample at x, k unknowns per window.
a = 1.9 with k = 1 is the local criterion (only n inside the support besides even numbers).
Nothing is ever guessed: a window whose feasible 0/1 patterns disagree on its first unknown stops the run.
Inputs: the zero list (Hypothesis Z(T, eps)), integer arithmetic. Lambda of even numbers (ln 2 at powers
of 2) and of perfect powers r^e (from the decision on r) need no primality decision.
"""
import math, sys, time, json
import numpy as np
from kbm import K, K0
from cert import sample

import os
S = os.environ.get('WORK', os.getcwd()) + '/'
ZFILE = os.environ.get('ZFILE', S + 'z2m.f64')
G = np.load(ZFILE) if ZFILE.endswith('.npy') else np.fromfile(ZFILE)
EPS = float(os.environ.get('EPS', '1e-5'))

def iroot(n, e):
    r = int(round(n ** (1.0 / e)))
    for q in (r - 1, r, r + 1):
        if q > 1 and q ** e == n:
            return q
    return 0

class Chain:
    def __init__(self, X):
        self.X = X
        self.lam = np.full(X + 400, np.nan)
        self.lam[0] = self.lam[1] = 0.0
        for n in range(2, X + 400, 2):
            self.lam[n] = math.log(2) if n & (n - 1) == 0 else 0.0
        self.powers = {}  # reduced root r -> list of r^e (e >= 2), odd
        for e in range(2, 40):
            r = 3
            while r ** e < X + 400:
                if not any(iroot(r, f) for f in range(2, 20)):
                    self.powers.setdefault(r, []).append(r ** e)
                r += 2
        self.is_pp = np.zeros(X + 400, dtype=bool)
        for r, lst in self.powers.items():
            for v in lst: self.is_pp[v] = True
        self.cache = {}
        self.maxT = 0.0
        self.log = []

    def commit(self, n, prime):
        self.lam[n] = math.log(n) if prime else 0.0
        for v in self.powers.get(n, []):
            self.lam[v] = math.log(n) if prime else 0.0

    def zs(self, x, a, beta, m, c):
        key = (x, a, beta, m, c)
        v = self.cache.get(key)
        if v is None:
            T = c * beta * x / a
            self.maxT = max(self.maxT, T)
            z, d, _ = sample(x, G, T, a, beta, m, EPS)
            v = self.cache[key] = (z, d)
        return v

    def unknowns(self, c, k):
        U = []
        n = c
        while len(U) < k + 1:
            if n % 2 == 1 and np.isnan(self.lam[n]) and not self.is_pp[n]:
                U.append(n)
            n += 2 if n % 2 else 1
        return U[:-1], U[-1]

    def step(self, c, a, beta, m, cc, k, h=1.0):
        U, w = self.unknowns(c, k)
        Ua = np.array(U, dtype=float)
        xs = []
        x = max(2.0, math.floor(c - a - 1))
        while x < w:
            if x * math.exp(-a / x) > 1.0 and abs(x * math.log(c / x)) < a and x * math.log(w / x) >= a:
                xs.append(float(x))
            x += h
        if not xs:
            return None, 'no samples'
        A = np.empty((len(xs), len(U))); y = np.empty(len(xs)); d = np.empty(len(xs))
        for i, x in enumerate(xs):
            z, dd = self.zs(x, a, beta, m, cc)
            lo, hi = math.floor(x * math.exp(-a / x)), math.ceil(x * math.exp(a / x))
            mm = np.arange(max(2, lo), hi + 1)
            v = x * np.log(mm / x)
            inside = np.abs(v) < a
            lv = self.lam[mm]
            unk = np.isnan(lv) & inside
            isU = np.isin(mm, Ua)
            if np.any(unk & ~isU):
                raise RuntimeError(f'unknown number inside the support of sample {x}')
            ok = inside & ~np.isnan(lv) & (lv > 0)
            known = float(np.sum(lv[ok] * np.sqrt(x / mm[ok]) * K(v[ok], a, beta, m)))
            A[i] = np.log(Ua) * np.sqrt(x / Ua) * K(x * np.log(Ua / x), a, beta, m)
            y[i] = z - known; d[i] = dd
        sols, trunc = self.enumerate(A, y, d)
        if trunc:
            return None, f'more than {len(sols)} feasible patterns'
        return (U, sols), None

    @staticmethod
    def enumerate(A, y, d, cap=200000):
        k = A.shape[1]
        lo_t, hi_t = y - d, y + d
        rest = np.hstack([np.cumsum(A[:, ::-1], axis=1)[:, ::-1], np.zeros((A.shape[0], 1))])
        sols = []; b = [0] * k; state = {'trunc': False}
        def dfs(j, acc):
            if state['trunc']: return
            if np.any(acc > hi_t) or np.any(acc + rest[:, j] < lo_t): return
            if j == k:
                sols.append(tuple(b))
                if len(sols) > cap: state['trunc'] = True
                return
            for v in (0, 1):
                b[j] = v; dfs(j + 1, acc + v * A[:, j])
            b[j] = 0
        dfs(0, np.zeros(A.shape[0]))
        return sols, state['trunc']

    def run(self, stages):
        c = 3; t0 = time.time(); maxsols = 0; decided = 0
        for (end, a, beta, m, cc, k) in stages:
            smax = 0
            while c < min(end, self.X):
                if not np.isnan(self.lam[c]) or self.is_pp[c] and not np.isnan(self.lam[c]):
                    c += 2; continue
                if self.is_pp[c]:
                    raise RuntimeError(f'perfect power {c} reached before its root was decided')
                res, err = self.step(c, a, beta, m, cc, k)
                if res is None:
                    return dict(stop=c, reason=err, decided=decided, secs=round(time.time() - t0))
                U, sols = res
                smax = max(smax, len(sols))
                j = 0
                while j < len(U) and len({s[j] for s in sols}) == 1:
                    self.commit(U[j], sols[0][j] == 1); decided += 1; j += 1
                if j == 0:
                    return dict(stop=c, reason=f'{len(sols)} feasible patterns disagree on {c}', decided=decided, secs=round(time.time() - t0))
                c = U[j - 1] + 2
            self.log.append(dict(end=end, a=a, beta=beta, m=m, c=cc, k=k, kappa=cc * beta / a, maxsols=smax))
        return dict(stop=None, decided=decided, maxT=self.maxT, secs=round(time.time() - t0))

if __name__ == '__main__':
    X = int(sys.argv[1]); stages = json.loads(sys.argv[2]); out = sys.argv[3] if len(sys.argv) > 3 else None
    C = Chain(X)
    r = C.run(stages)
    print(r); print(C.log)
    primes = [n for n in range(3, X) if n % 2 and C.lam[n] == math.log(n)]
    if out:
        with open(out, 'w') as f:
            f.write(json.dumps(dict(X=X, stages=stages, result=r, log=C.log, primes=primes)))
    print(f'{len(primes)} odd primes certified below {X}; zeros used up to height {C.maxT:.1f}')
