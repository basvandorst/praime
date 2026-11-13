"""Certified window decoding: proofs of primality from zeros below the local height.

Walk up over odd n. Everything below the current candidate c is already certified.
Unknowns U = the next k odd non-perfect-powers from c. Samples x (grid step h) are chosen so that
the kernel support (x e^{-a/x}, x e^{a/x}) contains c and no odd non-perfect-power beyond U.
Each sample gives  y(x) = sum_{j in U} b_j ln(u_j) sqrt(x/u_j) K(x ln(u_j/x)) + e(x),  |e(x)| <= delta(x).
All b in {0,1}^k with |A b - y| <= delta (componentwise) are enumerated (K >= 0, so interval pruning
is exact). c is certified if all of them agree on b_c; the leading run of agreeing unknowns is committed.
If they disagree, the walk stops: nothing is ever guessed.
"""
import math, sys, time
import numpy as np
from kbm import K, K0
from cert import sample
from sympy import isprime, perfect_power, integer_nthroot

import os
S = os.environ.get('WORK', os.getcwd()) + '/'
ZFILE = os.environ.get('ZFILE', S + 'z2m.f64')
G = np.load(ZFILE) if ZFILE.endswith('.npy') else np.fromfile(ZFILE)
EPS = float(os.environ.get('EPS', '1e-5'))

class Walk:
    def __init__(self, a, beta, m, c, k, h, start, truth_below=False):
        self.a, self.beta, self.m, self.c, self.k, self.h = a, beta, m, c, k, h
        self.k0 = K0(beta, m)
        self.lam = {}            # certified Lambda for odd numbers
        self.cache = {}
        self.start = start
        if truth_below:          # exploratory runs only: take the true values below start
            for n in range(3, start, 2):
                self.lam[n] = self.true_lam(n)
        self.nodes = 0

    @staticmethod
    def true_lam(n):
        pp = perfect_power(n)
        if pp:
            r, e = pp
            f = __import__('sympy').factorint(n)
            return math.log(next(iter(f))) if len(f) == 1 else 0.0
        return math.log(n) if isprime(n) else 0.0

    def known_lam(self, m):
        """Lambda(m) if it is determined without new primality decisions, else None."""
        if m % 2 == 0:
            return math.log(2) if m & (m - 1) == 0 else 0.0
        if m in self.lam:
            return self.lam[m]
        pp = perfect_power(m)
        if pp:
            r, e = pp
            # m = r^e with e maximal; Lambda(m) = ln q if r is a power of a prime q, else 0
            while True:
                pr = perfect_power(r)
                if not pr: break
                r = pr[0]
            if r in self.lam or r < 3:
                v = self.lam.get(r, 0.0)
                return math.log(r) if v == math.log(r) else 0.0
            return None
        return None

    def zs(self, x):
        v = self.cache.get(x)
        if v is None:
            T = self.c * self.beta * x / self.a
            z, d, _ = sample(x, G, T, self.a, self.beta, self.m, EPS)
            v = (z, d)
            self.cache[x] = v
        return v

    def unknown_list(self, c):
        U = []
        n = c
        while len(U) < self.k + 1:
            if n % 2 == 1 and self.known_lam(n) is None:
                U.append(n)
            n += 1
        return U[:-1], U[-1]

    def step(self, c):
        a = self.a
        U, w = self.unknown_list(c)
        # sample grid: c inside the support, w outside
        xs = []
        x = math.floor((c - a) / self.h) * self.h
        while x < w:
            if x > a + 2 and abs(x * math.log(c / x)) < a and x * math.log(w / x) >= a:
                xs.append(x)
            x += self.h
        if not xs:
            return None, 'no samples'
        rows, ys, ds = [], [], []
        for x in xs:
            z, d = self.zs(x)
            lo, hi = x * math.exp(-a / x), x * math.exp(a / x)
            known = 0.0
            for mm in range(max(2, math.floor(lo)), math.ceil(hi) + 1):
                if mm in U: continue
                lv = self.known_lam(mm)
                if lv is None:
                    if abs(x * math.log(mm / x)) < a:
                        raise RuntimeError(f'unknown {mm} inside support of {x}')
                    continue
                if lv:
                    known += lv * math.sqrt(x / mm) * float(K(x * math.log(mm / x), a, self.beta, self.m))
            row = [math.log(u) * math.sqrt(x / u) * float(K(x * math.log(u / x), a, self.beta, self.m)) for u in U]
            rows.append(row); ys.append(z - known); ds.append(d)
        A = np.array(rows); y = np.array(ys); d = np.array(ds)
        sols = self.enumerate(A, y, d)
        if len(sols) > 200000:
            return None, f'more than 200000 feasible patterns at {c}: undecided'
        return (U, sols, A, y, d, xs), None

    def enumerate(self, A, y, d, cap=200000):
        k = A.shape[1]
        lo_t, hi_t = y - d, y + d
        rest = np.cumsum(A[:, ::-1], axis=1)[:, ::-1]   # rest[:, j] = sum of columns j..k-1
        rest = np.hstack([rest, np.zeros((A.shape[0], 1))])
        sols = []
        b = [0] * k
        def dfs(j, acc):
            self.nodes += 1
            if len(sols) > cap: return
            if np.any(acc > hi_t) or np.any(acc + rest[:, j] < lo_t):
                return
            if j == k:
                sols.append(tuple(b)); return
            for v in (0, 1):
                b[j] = v
                dfs(j + 1, acc + v * A[:, j])
            b[j] = 0
        dfs(0, np.zeros(A.shape[0]))
        return sols

    def run(self, end, verbose=False):
        c = self.start | 1
        t0 = time.time(); decided = 0; maxsols = 0
        while c < end:
            if self.known_lam(c) is not None:
                c += 2; continue
            res, err = self.step(c)
            if res is None:
                return dict(stop=c, reason=err, decided=decided, maxsols=maxsols, secs=time.time() - t0)
            U, sols, A, y, dd, xs = res
            truth = tuple(1 if self.true_lam(u) > 0 else 0 for u in U)
            if truth not in sols:
                raise AssertionError(f'bound violated at {c}')
            maxsols = max(maxsols, len(sols))
            # leading run of unanimous unknowns
            j = 0
            while j < len(U) and len({s[j] for s in sols}) == 1:
                self.lam[U[j]] = math.log(U[j]) if sols[0][j] else 0.0
                decided += 1
                j += 1
            if j == 0:
                return dict(stop=c, reason=f'{len(sols)} patterns disagree on {c}', decided=decided, maxsols=maxsols, secs=time.time() - t0)
            c = U[j - 1] + 2 if j < len(U) else U[-1] + 2
        wrong = sum(1 for n, v in self.lam.items() if n >= self.start and (v > 0) != (self.true_lam(n) > 0))
        return dict(stop=None, decided=decided, wrong=wrong, maxsols=maxsols, nodes=self.nodes, secs=time.time() - t0)

if __name__ == '__main__':
    a, beta, m, c, k, h, start, end = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8])
    tb = len(sys.argv) > 9 and sys.argv[9] == 'truth'
    W = Walk(a, beta, m, c, k, h, start, truth_below=tb)
    r = W.run(end)
    kappa = c * beta / a
    print(f'a={a} beta={beta} m={m} c={c} k={k} h={h}: height T = {kappa:.3f} x = {kappa/math.pi:.3f} pi x; from {start} to {end}: {r}')
