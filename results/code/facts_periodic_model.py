"""Parameter-free prime model: every prime q <= sqrt x deletes one random residue class (periodic deletion),
and each survivor is prime with probability rho = 1/(ln x * prod_{q<=sqrt x}(1-1/q)) -> e^gamma/2 = 0.8905 (Mertens' paradox).
Usage: facts_periodic_model.py tail   (gaps >= G after primes in [1e8, 4e9], vs real)
       facts_periodic_model.py max B G...   (maximal gap below B)
       facts_periodic_model.py var    (Var/mean of window counts near 1e9, vs real)"""
import sys, math, numpy as np
import facts_maxgap_is as m
def T(G, x, nsamp=1500, seed=1):
    y = math.sqrt(x); rho = 1.0 / (math.log(x) * math.exp(m.logprod(1, y)))
    small = m.PR[m.PR < G]; t = np.arange(G, dtype=float); f = (1 - rho) ** t
    for q in m.PR[(m.PR >= G) & (m.PR <= y)][::-1]:          # q >= G: at most one window position per q
        pk = np.minimum(t / (q - 1), 1.0); fn = (1 - pk) * f; fn[1:] += pk[1:] * f[:-1]; f = fn
    logf = np.log(f); beta = -(logf[-1] - logf[0]) / (G - 1)
    rng = np.random.default_rng(seed); pos0 = np.arange(1, G); vals = np.empty(nsamp)
    for s in range(nsamp):                                   # q < G: tilted residue sampling
        alive = pos0; logw = 0.0
        for q in small:
            q = int(q); k = np.bincount((-alive) % q, minlength=q)[1:]
            lw = beta * k; lw -= lw.max(); p = np.exp(lw); p /= p.sum()
            i = rng.choice(q - 1, p=p); logw += -math.log(q - 1) - math.log(p[i])
            alive = alive[(alive + i + 1) % q != 0]
        vals[s] = logw + logf[len(alive)]
    mx = vals.max(); return np.exp(vals - mx).mean() * math.exp(mx), 0.0, rho
def count(G, A, B, npts=6, **kw):
    xs = np.geomspace(A, B, npts); lT = np.log([T(G, x, **kw)[0] for x in xs]); grid = np.geomspace(A, B, 400)
    return np.trapezoid(np.exp(np.interp(np.log(grid), np.log(xs), lT)) / np.log(grid), grid)
def var():
    X, L = 10**9, 10**8; r = math.isqrt(X + L); Hs = [100, 1000, 10**4, 10**5, 10**6]
    real = np.ones(L, bool); S = np.ones(L, bool); rng = np.random.default_rng(7)
    for q in m.PR[m.PR <= r]:
        q = int(q); real[(-X) % q::q] = False; S[rng.integers(q)::q] = False
    rho = 1.0 / (math.log(X) * math.exp(m.logprod(1, r)))
    def vm(a):
        s = a.astype(float); t = np.arange(len(s)); c = np.polyfit(t, s, 1); return np.var(s - np.polyval(c, t)), s.mean()
    for H in Hs:
        vr, mr = vm(real[:L//H*H].reshape(-1, H).sum(1)); vs, ms = vm(S[:L//H*H].reshape(-1, H).sum(1))
        print(f"H={H}: real {vr/mr:.3f}  model {(rho**2*vs + rho*(1-rho)*ms)/(rho*ms):.3f}", flush=True)
if __name__ == "__main__":
    if sys.argv[1] == "tail":
        for G, r in {100: 915491, 200: 2825, 240: 260}.items():
            print(G, r, f"model/real {count(G, 1e8, 4e9) / r:.3f}", flush=True)
    elif sys.argv[1] == "var": var()
    else:
        B = float(sys.argv[2])
        for G in map(int, sys.argv[3:]):
            N = count(G, B / 100, B, npts=5, nsamp=800); print(f"below {B:.3g}: G={G} N(>=G)={N:.3g} P(max<G)={math.exp(-N):.3f}", flush=True)
