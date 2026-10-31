"""Maximal prime gaps predicted by the sieve + coin model, without fitting.

Model (FINDINGS: 'one explicit law for the next prime'): integers coprime to all q <= y = x^0.405 survive,
each survivor is prime with probability rho = 1/(ln x prod_{q<=y}(1-1/q) e^gamma omega(1/0.405)).
P(gap >= G after a prime n) = E[(1-rho)^S], S = survivors among n+1..n+G-1.
 - primes q < G: residues of n sampled with an exponentially tilted proposal (importance sampling);
 - primes G <= q <= 2e5: exact death process (each q removes at most ONE window position), backward DP;
 - primes above 2e5: independent survival (second-order error sum t^2/q^2 is negligible).
Usage: python3 facts_maxgap_is.py calib | python3 facts_maxgap_is.py <B> <G> ...  (counts primes in [B/100, B])
"""
import numpy as np, math, sys
GAMMA = 0.5772156649015329
def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]
PR = primes_upto(3 * 10**7)
LOGC = np.cumsum(np.log1p(-1.0 / PR))           # log prod_{q<=PR[i]} (1-1/q)
def logprod(lo, hi):                             # log prod_{lo<q<=hi}(1-1/q)
    a = np.searchsorted(PR, lo, 'right'); b = np.searchsorted(PR, hi, 'right')
    return LOGC[b-1] - (LOGC[a-1] if a > 0 else 0.0)
def omega(u):                                    # Buchstab, valid 2<=u<=3 here (enough for y=x^0.405)
    assert 2 <= u <= 3
    return (1 + math.log(u - 1)) / u
def T(G, x, ysh=0.405, nsamp=1500, seed=1, Q1=200000, exact=True):
    y = x ** ysh
    rho = 1.0 / (math.log(x) * math.exp(logprod(1, y)) * math.exp(GAMMA) * omega(1 / ysh))
    small = PR[PR < G]
    t = np.arange(G, dtype=float)
    if exact:
        hi = min(Q1, y)
        f = (1 - math.exp(logprod(hi, y)) * rho) ** t if hi < y else (1 - rho) ** t
        for q in PR[(PR >= G) & (PR <= hi)][::-1]:
            pk = np.minimum(t / (q - 1), 1.0)
            fn = (1 - pk) * f; fn[1:] += pk[1:] * f[:-1]; f = fn
        logf = np.log(f)
    else:
        logf = t * math.log1p(-math.exp(logprod(G, y)) * rho)
    beta = -(logf[-1] - logf[0]) / (G - 1)
    rng = np.random.default_rng(seed); pos0 = np.arange(1, G); vals = np.empty(nsamp)
    for s in range(nsamp):
        alive = pos0; logw = 0.0
        for q in small:
            q = int(q)
            k = np.bincount((-alive) % q, minlength=q)[1:]
            lw = beta * k; lw -= lw.max(); p = np.exp(lw); p /= p.sum()
            i = rng.choice(q - 1, p=p)
            logw += -math.log(q - 1) - math.log(p[i])
            alive = alive[(alive + i + 1) % q != 0]
        vals[s] = logw + logf[len(alive)]
    m = vals.max(); e = np.exp(vals - m)
    return e.mean() * math.exp(m), e.std() / math.sqrt(nsamp) * math.exp(m), rho
def count(G, A, B, npts=6, **kw):
    xs = np.geomspace(A, B, npts); lT = np.log([T(G, x, **kw)[0] for x in xs])
    grid = np.geomspace(A, B, 400)
    return np.trapezoid(np.exp(np.interp(np.log(grid), np.log(xs), lT)) / np.log(grid), grid)
def cramer(G, A, B):
    grid = np.geomspace(A, B, 2000); return np.trapezoid(np.exp(-G / np.log(grid)) / np.log(grid), grid)
if __name__ == "__main__":
    if sys.argv[1] == "calib":   # real counts of gaps >= G after primes in [1e8, 4e9] (facts round 14)
        real = {100: 915491, 200: 2825, 240: 260}
        for G, r in real.items():
            print(G, r, f"model/real {count(G, 1e8, 4e9) / r:.3f}", flush=True)
    else:
        B = float(sys.argv[1])
        for G in map(int, sys.argv[2:]):
            N = count(G, B / 100, B, npts=5, nsamp=800)
            print(f"below {B:.3g}: G={G} model N(>=G)={N:.3g} P(max<G)={math.exp(-N):.3f} Cramer N={cramer(G, B/100, B):.3g}", flush=True)
