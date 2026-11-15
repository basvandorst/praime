"""Sacks spiral: the streaks are the curves n = x^2 + d (they converge to the ray at angle 0, at
height pi*d) and n = x^2 + x + c (to the ray at angle pi, at height pi*(c - 1/4)).
Observed prime counts vs the Bateman-Horn prediction C(f) * sum_x 1/ln f(x)."""
import numpy as np, common as C, json
isp = C.sieve()
N = C.N
P = np.nonzero(C.sieve(10**6))[0][1:]   # odd primes up to 1e6 for the Euler product

def legendre_vec(D, p):
    """(D/p) for odd primes p (vector), D int."""
    a = np.mod(D, p).astype(np.int64); e = (p-1)//2; r = np.ones_like(p); b = a.copy()
    while e.any():
        r = np.where(e & 1, (r*b) % p, r); b = (b*b) % p; e >>= 1
    return np.where(a == 0, 0, np.where(r == 1, 1, -1))

def hl_const(disc, two_factor):
    leg = legendre_vec(disc, P)
    omega = np.where(leg == 0, 1, 1+leg)
    return two_factor*np.exp(np.sum(np.log1p(-omega/P) - np.log1p(-1/P)))

rows = []
X0 = 10
for side in ('right', 'left'):
    for c in range(-60, 61):
        if side == 'right':
            disc = -4*c; two = 1.0
            f = lambda x: x*x + c
            red = c <= 0 and int(round(np.sqrt(-c)))**2 == -c
        else:
            disc = 1-4*c; two = 0.0 if c % 2 == 0 else 2.0
            f = lambda x: x*x + x + c
            red = disc >= 0 and int(round(np.sqrt(disc)))**2 == disc
        x = np.arange(X0, 800, dtype=np.int64); v = f(x); v = v[(v >= 2) & (v <= N)]
        obs = int(isp[v].sum())
        Cf = 0.0 if red else hl_const(disc, two)
        exp_ = Cf*np.sum(1/np.log(v))
        rows.append(dict(side=side, c=c, n_values=len(v), obs=obs, C=float(Cf), exp=float(exp_), reducible=bool(red)))
obs = np.array([r['obs'] for r in rows]); ex = np.array([r['exp'] for r in rows])
ok = ex > 0
chi = ((obs[ok]-ex[ok])**2/ex[ok])
print(f"curves: {len(rows)}; reducible (C=0): {sum(r['reducible'] for r in rows)}; their primes: {obs[~ok].sum()}")
print(f"irreducible: sum obs {obs[ok].sum()}, sum BH {ex[ok].sum():.1f}; chi2/dof = {chi.sum():.1f}/{ok.sum()} ; corr = {np.corrcoef(obs[ok], ex[ok])[0,1]:.4f}")
print(f"slope obs~BH (through 0): {np.sum(obs[ok]*ex[ok])/np.sum(ex[ok]**2):.4f}")
z = (obs[ok]-ex[ok])/np.sqrt(ex[ok]); print(f"max |z| = {np.abs(z).max():.2f}")
for side in ('right', 'left'):
    rs = [r for r in rows if r['side'] == side]
    rs.sort(key=lambda r: -r['C'])
    print(side, 'richest:', [(r['c'], round(r['C'], 3), r['obs'], round(r['exp'], 1)) for r in rs[:4]])
    print(side, 'poorest irreducible:', [(r['c'], round(r['C'], 3), r['obs'], round(r['exp'], 1)) for r in rs if r['C'] > 0][-3:])
    print(side, 'C=0 curves:', [r['c'] for r in rs if r['C'] == 0][:12])
# left side parity: c even -> all values even -> no primes
L = [r for r in rows if r['side'] == 'left']
print('left c even: primes', sum(r['obs'] for r in L if r['c'] % 2 == 0), ' left c odd: primes', sum(r['obs'] for r in L if r['c'] % 2))
# brightness profile across the streaks (height above the ray) near both rays
n = np.arange(1, N+1); r = np.sqrt(n); th = 2*np.pi*r
for side, ph in (('right(0)', 0.0), ('left(pi)', np.pi)):
    dth = np.angle(np.exp(1j*(th-ph)))
    y = r*np.sin(dth)
    sel = (np.abs(dth) < 0.15) & (r > 300) & (np.abs(y) < 40)
    yb = y[sel]; pr = isp[n[sel]]
    # Fourier coefficient of prime fraction at period 2*pi in height (left: c parity) and pi
    for per in (np.pi, 2*np.pi):
        e = np.exp(2j*np.pi*yb/per)
        ap = np.abs(np.mean(e[pr])); ai = np.abs(np.mean(e))
        print(f"{side}: height period {per:.3f}px: |mean e| primes {ap:.3f}, all numbers {ai:.3f}")
json.dump(rows, open('out/sacks_curves.json', 'w'), indent=0)
