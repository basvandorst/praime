import sys, math, numpy as np
N = int(float(sys.argv[1]))
gam = 0.57721566490153286061
lim = 2 + gam - math.log(4 * math.pi)
s = np.ones(N + 1, dtype=bool); s[:2] = False
for p in range(2, int(N ** 0.5) + 1):
    if s[p]: s[p * p::p] = False
P = np.nonzero(s)[0]; del s
ns, ws = [P.astype(np.float64)], [np.log(P.astype(np.float64))]
for k in range(2, 64):
    q = P[P <= int(round(N ** (1 / k))) + 1].astype(np.int64) ** k
    pk = P[P <= int(round(N ** (1 / k))) + 1]
    m = q <= N
    if not m.any(): break
    ns.append(q[m].astype(np.float64)); ws.append(np.log(pk[m].astype(np.float64)))
n = np.concatenate(ns); w = np.concatenate(ws); o = np.argsort(n); n = n[o]; w = w[o]
A = np.cumsum((w / n).astype(np.longdouble)); psi = np.cumsum(w.astype(np.longdouble))
nl = n.astype(np.longdouble)
G = nl * A - psi - nl * np.log(nl) + (1 + gam) * nl
dev = (G - np.longdouble(math.log(2 * math.pi))) / np.sqrt(nl)
ystar = np.exp(A + gam)
nxt = np.append(nl[1:], np.longdouble(N))
inside = (ystar > nl) & (ystar < nxt)
Gs = ystar[inside] * A[inside] - psi[inside] - ystar[inside] * np.log(ystar[inside]) + (1 + gam) * ystar[inside]
devs = (Gs - np.longdouble(math.log(2 * math.pi))) / np.sqrt(ystar[inside])
big = nl >= 100
print(f"N={N:.0e} prime powers {len(n)}; bound (2+gamma-log 4pi) = {lim:.10f}")
print(f"y>=100: max (G-log2pi)/sqrt y at prime powers {float(dev[big].max()):+.6f} at y={float(nl[big][dev[big].argmax()]):.0f}; min {float(dev[big].min()):+.6f} at y={float(nl[big][dev[big].argmin()]):.0f}")
bs = ystar[inside] >= 100
print(f"interior maxima: max {float(devs[bs].max()):+.6f}; overall |dev| max {max(float(abs(dev[big]).max()), float(abs(devs[bs]).max())):.6f} (bound {lim:.6f})")
for y0 in [10, 100, 1000]:
    m = nl < y0
    print(f"  y<{y0}: max |dev| {float(abs(dev[m & (nl>=2)]).max()):.4f}")
