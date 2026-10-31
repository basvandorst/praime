"""Expected number of runs of CONSECUTIVE primes that follow a fixed pattern.

A run at p means p + h is prime for every h in H, and no other number in (p, p + max H) is prime.
With --exact, the next prime after p + max H must also not be p + EXT (for chains like
Euler-Rodriguez, where the pattern would continue).

Model: an exact random sieve to Y (Y just above the span), with a coin rho_Y above it.
- p mod q is drawn uniformly from the classes that keep every p + h coprime to q, independently
  for every q <= Y (CRT).
- A number that survives the sieve to Y is prime with probability
  rho = 1 / (ln t * prod_{q<=Y} (1 - 1/q) * e^gamma * omega(ln t / ln Y)), where omega is Buchstab's function.
- The pattern elements get the Hardy-Littlewood factor for q > Y. P(nothing else prime) = E[(1 - rho)^m],
  where m is the number of survivors strictly inside the span.

Calibration below 1.4e11 (real counts from a C sieve):
gaps 2 and 6: 1.0003 and 1.0001; Euler-Rodriguez chains of exactly k = 3..8: 0.9996, 0.9971, 0.9966, 1.0103, 0.9766, 1.0137;
CPAP-4/5/6 with d = 30: 1.009, 0.979, 0.998.
Two traps: Y far below sqrt(t) adds a spurious Buchstab factor per element; averaging the HL factors over
p's classes instead of sampling them (Jensen) overcounts long spans (CPAP-6: +24 %).
Second order: for the sieving primes q > Y (Sigma = sum of 1/q^2 over Y < q <= sqrt t), an in-between survivor is less
likely prime next to the k pattern primes (+ m rho k Sigma/(1-rho)) and two in-between survivors are negatively
correlated (- C(m,2) rho^2 Sigma/(1-rho)^2). Without this term long spans were overcounted (CPAP-3, d = 60/90/120/150:
1.016/1.032/1.071/1.05). With it, and averaged over 5 seeds of 8000 samples (one run has +-3 % sampling noise at
span 240): CPAP-3 d = 60/90/120: 1.001 +- 0.003, 0.990 +- 0.009, 1.016 +- 0.011.

Examples:
  python3 facts_pattern_model.py --pattern 0,2,6,12,20,30,42,56,72,90,110,132,156,182,210 --to 1e17 1e18
  python3 facts_pattern_model.py --pattern 0,2,6,12,20 --exact 30 --to 1.4e11
  python3 facts_pattern_model.py --ap 7 --d 210,420,630 --to 1e20 1e21
"""
import argparse
import math

import numpy as np

GAMMA = 0.5772156649015329


def primes_upto(n):
    b = np.ones(n + 1, bool)
    b[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if b[i]:
            b[i * i::i] = False
    return np.nonzero(b)[0]


def omega_table(vmax=30.0, h=2e-4):
    n = int((vmax - 1) / h) + 2
    om = np.zeros(n)
    for i in range(n):
        u = 1 + i * h
        if u < 2:
            om[i] = 1 / u
        elif u <= 3:
            om[i] = (1 + math.log(u - 1)) / u
        else:
            om[i] = ((u - h) * om[i - 1] + h * om[int(round((u - 2) / h))]) / u
    return 1 + h * np.arange(n), om


class Model:
    def __init__(self, Y, tmax=1e30, samples=4000, seed=1):
        self.Y = Y
        self.qs = [int(q) for q in primes_upto(Y)]
        big = primes_upto(3_000_000).astype(float)
        self.big = big[big > Y]
        self.ns = samples
        self.rng = np.random.default_rng(seed)
        self.U = np.linspace(math.log(11), math.log(tmax), 1500)
        T = np.exp(self.U)
        lt = np.log(T)
        vg, om = omega_table()
        v = np.clip(lt / math.log(Y), 1.0, vg[-1])
        prodY = math.prod(1 - 1 / q for q in self.qs)
        self.T = T
        self.rho = np.minimum(1 / (lt * prodY * math.exp(GAMMA) * np.interp(v, vg, om)), 1.0)
        # sum of 1/q^2 over the sieving primes Y < q <= sqrt(t): drives the second-order (pair) corrections
        inv2 = np.concatenate([[0], np.cumsum(1 / self.big ** 2)])
        self.sig = inv2[np.searchsorted(self.big, np.sqrt(T), side="right")]

    def density(self, H, ext=None):
        """Expected runs per unit t at every grid point t."""
        H = np.array(sorted(set(H))) - min(H)
        k, D = len(H), int(H.max())
        inner = np.array([n for n in range(1, D) if n not in set(H.tolist())], int)
        tail = np.arange(D + 1, ext) if ext else np.array([], int)
        log_a, allowed = 0.0, {}
        for q in self.qs:
            bad = set(((-H) % q).tolist())
            al = np.array([a for a in range(q) if a not in bad])
            if len(al) == 0:
                return np.zeros_like(self.T)
            allowed[q], log_a = al, log_a + math.log(len(al) / q)
        log_big = float(np.sum(np.log1p(-k / self.big) - k * np.log1p(-1 / self.big)))
        groups = {}
        for _ in range(self.ns):
            alive = np.ones(D + 1 + (ext or 0), bool)
            for q, al in allowed.items():
                alive[(-int(al[self.rng.integers(len(al))])) % q::q] = False
            key = (int(alive[inner].sum()), bool(ext and alive[ext]), int(alive[tail].sum()) if ext else 0)
            groups[key] = groups.get(key, 0) + 1
        E = np.zeros_like(self.T)
        r, sg = self.rho, self.sig
        for (m, e, m2), c in groups.items():
            # second order in q > Y: each in-between survivor is a little less likely prime next to the k pattern
            # primes (+), and two in-between survivors are negatively correlated (-)
            corr = np.where(r < 0.99, m * r * k * sg / (1 - r + 1e-12) - m * (m - 1) / 2 * r ** 2 * sg / (1 - r + 1e-12) ** 2, 0.0)
            f = (1 - r) ** m * np.exp(corr)
            if e:
                f = f * (1 - self.rho * (1 - self.rho) ** m2)
            E += c / self.ns * f
        return math.exp(log_a + log_big) * self.rho ** k * E

    def cumulative(self, dens):
        f = dens * self.T
        return np.concatenate([[0], np.cumsum((f[1:] + f[:-1]) / 2 * np.diff(self.U))])

    def count(self, dens, X, lo=11.0):
        c = self.cumulative(dens)
        return float(np.interp(math.log(X), self.U, c) - np.interp(math.log(lo), self.U, c))

    def median_first(self, dens):
        c = self.cumulative(dens)
        i = np.searchsorted(c, math.log(2))
        return float(self.T[min(i, len(self.T) - 1)])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pattern", help="offsets, e.g. 0,2,6,12")
    ap.add_argument("--exact", type=int, help="offset whose primality would extend the run (exactly-k chains)")
    ap.add_argument("--ap", type=int, help="CPAP-k: k consecutive primes in arithmetic progression")
    ap.add_argument("--d", default=None, help="common differences for --ap, e.g. 210,420")
    ap.add_argument("--to", nargs="+", type=float, default=[1e12])
    ap.add_argument("--Y", type=int, default=None, help="exact sieve bound (default: just above the span, at least 200)")
    ap.add_argument("--samples", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    if a.ap:
        ds = [int(x) for x in a.d.split(",")]
        patterns = [[i * d for i in range(a.ap)] for d in ds]
        ext = None
    else:
        patterns = [[int(x) for x in a.pattern.split(",")]]
        ext = a.exact
    span = max(max(p) for p in patterns) + (ext or 0)
    Y = a.Y or max(200, span + 1)
    m = Model(Y, tmax=max(a.to) * 10, samples=a.samples, seed=a.seed)
    dens = sum(m.density(p, ext) for p in patterns)
    for X in a.to:
        print(f"expected runs below {X:.3g}: {m.count(dens, X):.4g}")
    print(f"median first occurrence: {m.median_first(dens):.3g}   (Y = {Y}, {a.samples} samples)")


if __name__ == "__main__":
    main()
