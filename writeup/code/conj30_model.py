"""Maximal prime gaps predicted by the sieve + coin model, without fitting; extended past 3e7 with Mertens and run in log space, so it reaches x = 1e40.

Model (FINDINGS: 'one explicit law for the next prime'): integers coprime to all q <= y = x^0.405 survive,
each survivor is prime with probability rho = 1/(ln x prod_{q<=y}(1-1/q) e^gamma omega(1/0.405)).
P(gap >= G after a prime n) = E[(1-rho)^S], S = survivors among n+1..n+G-1.
 - primes q < G: residues of n sampled with an exponentially tilted proposal (importance sampling);
 - primes G <= q <= 2e5: exact death process (each q removes at most ONE window position), backward DP;
 - primes above 2e5: independent survival (second-order error sum t^2/q^2 is negligible).
Usage: python3 conj30_model.py trend | firoozbakht | csg | T <x> <G> [nsamp]
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
PMAX = float(PR[-1])
def logprod(lo, hi):
    if hi > PMAX:
        return logprod(lo, PMAX) + math.log(math.log(PMAX) / math.log(hi))
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
        base = math.log1p(-math.exp(logprod(hi, y)) * rho) if hi < y else math.log1p(-rho)
        logf = t * base
        for q in PR[(PR >= G) & (PR <= hi)][::-1]:
            pk = np.minimum(t / (q - 1), 1.0)
            with np.errstate(divide='ignore'):
                a = np.log1p(-pk) + logf
                b = np.log(pk) + np.concatenate(([-np.inf], logf[:-1]))
            logf = np.logaddexp(a, b)
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

def N_window(G, x, nsamp=200):
    return count(int(G), x / 1e4, x, npts=3, nsamp=nsamp)

def trend(xs=(1e16, 4e18, 1e20, 1e25, 1e30, 1e40), real={1e16: 1132, 4e18: 1476}, nsamp=200):
    for x in xs:
        L2 = math.log(x) ** 2
        lo, hi, cache = 0.65 * L2, 1.0 * L2, {}
        for _ in range(6):
            mid = int(0.5 * (lo + hi)); cache[mid] = N_window(mid, x, nsamp)
            if cache[mid] > math.log(2): lo = mid
            else: hi = mid
        Gs = np.array(sorted(cache)); a, b = np.polyfit(Gs, np.log([cache[g] for g in Gs]), 1)
        med = (math.log(math.log(2)) - b) / a; q10 = (math.log(-math.log(0.9)) - b) / a; q90 = (math.log(-math.log(0.1)) - b) / a
        extra = ""
        if x in real:
            extra = f"  real {real[x]} (ratio {real[x]/L2:.3f}) at model quantile {math.exp(-N_window(real[x], x, nsamp)):.2f}"
        print(f"x={x:.0e} ln^2x={L2:.0f}: median max gap {med:.0f} (ratio {med/L2:.3f}), 10-90% {q10/L2:.3f}-{q90/L2:.3f}, dlnN/dG {a:.4f} (Cramer {-1/math.log(x):.4f}){extra}", flush=True)

def density_integral(ts, Gfun, nsamp=200, label=""):
    ls, ld, cum = [], [], 0.0
    for t in ts:
        L = math.log(t); G = Gfun(L); v, e, rho = T(G, t, nsamp=nsamp)
        print(f"  t={t:.0e} G={G} (ratio {G/L/L:.4f}) P(gap>=G)={v:.2e} (rel {e/v:.2f}) dN/dln t={v*t/L:.2e}", flush=True)
        ls.append(L); ld.append(math.log(v * t / L))
        if len(ls) > 1:
            cum += 0.5 * (math.exp(ld[-2]) + math.exp(ld[-1])) * (ls[-1] - ls[-2])
            print(f"    cumulative {label} N(X=10^{L/math.log(10):.0f}) = {cum:.3e}", flush=True)
    return cum

def firoozbakht(nsamp=200):
    print("expected number of gaps g >= ln^2 p - ln p - 1 (Kourbatov's necessary condition for Firoozbakht) among primes above 4e18:")
    density_integral((4e18, 1e20, 1e22, 1e25, 1e30, 1e35, 1e40), lambda L: int(math.floor(L * L - L - 1)), nsamp, "violations")

def csg(nsamp=200):
    print("expected number of gaps with g/ln^2 p >= r among primes in [1e12, 4e18]:")
    for r in (0.80, 0.85, 0.90, 0.9206, 0.95, 1.00):
        n = density_integral((1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 4e18), lambda L, r=r: int(round(r * L * L)), nsamp, f"r={r}")
        print(f"r = {r}: N = {n:.3f}, P(none) = {math.exp(-n):.3f}", flush=True)

if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "trend": trend()
    elif cmd == "firoozbakht": firoozbakht()
    elif cmd == "csg": csg()
    elif cmd == "T":
        x = float(sys.argv[2]); G = int(sys.argv[3]); ns = int(sys.argv[4]) if len(sys.argv) > 4 else 300
        v, e, rho = T(G, x, nsamp=ns)
        print(f"x={x:.0e} G={G} ratio={G/math.log(x)**2:.3f}: P(gap >= G) = {v:.3e} +- {e:.1e}, rho = {rho:.4f}")
