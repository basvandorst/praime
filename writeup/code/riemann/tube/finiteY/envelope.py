# Upper envelope of H(y)=G(y)-ln2pi+T(y) over blocks of u=ln y: c_j >= max H(y)/sqrt(y), y in block j.
# G is linear+(-y ln y) between prime powers; max of G on [n_i,n_{i+1}) at y*=exp(A+gamma) (concave).
# Rigorous-ish envelope: per interval, Gmax/sqrt(n_i) if Gmax-ln2pi>0 else /sqrt(n_{i+1}); T(y)<=1/y^2 added.
import sys, math, numpy as np
N = int(float(sys.argv[1])); W = 0.05
gam = 0.57721566490153286061; L2P = math.log(2*math.pi)
s = np.ones(N+1, dtype=bool); s[:2] = False
for p in range(2, int(N**0.5)+1):
    if s[p]: s[p*p::p] = False
P = np.nonzero(s)[0]; del s
ns, ws = [P.astype(np.float64)], [np.log(P.astype(np.float64))]
for k in range(2, 64):
    pk = P[P <= int(round(N**(1/k)))+1]; q = pk.astype(np.int64)**k; m = q <= N
    if not m.any(): break
    ns.append(q[m].astype(np.float64)); ws.append(np.log(pk[m].astype(np.float64)))
n = np.concatenate([[1.0]]+ns); w = np.concatenate([[0.0]]+ws); del ns, ws, P; o = np.argsort(n, kind='stable'); n=n[o]; w=w[o]
A = np.cumsum((w/n).astype(np.longdouble)); psi = np.cumsum(w.astype(np.longdouble)); nl = n.astype(np.longdouble)
nxt = np.append(nl[1:], np.longdouble(N))
del w, o
Gl = (nl*A - psi - nl*np.log(nl) + (1+gam)*nl - L2P).astype(np.float64)
Gr = (nxt*A - psi - nxt*np.log(nxt) + (1+gam)*nxt - L2P).astype(np.float64)
ys = np.clip(np.exp(A+gam), nl, nxt)
Gs = (ys*A - psi - ys*np.log(ys) + (1+gam)*ys - L2P).astype(np.float64)
del A, psi
ys = ys.astype(np.float64); nl = nl.astype(np.float64); nxt = nxt.astype(np.float64)
Gmax = np.maximum(np.maximum(Gl, Gr), Gs) + 1/nl**2
val = np.where(Gmax > 0, Gmax/np.sqrt(nl), Gmax/np.sqrt(nxt)); del Gmax
dev = np.maximum(Gl/np.sqrt(nl), Gs/np.sqrt(ys))
devmin = np.minimum(Gl/np.sqrt(nl), Gr/np.sqrt(nxt)); del Gl, Gr, Gs, ys
uL = np.log(nl.astype(np.float64)); j = np.floor(uL/W).astype(np.int64)
nb = int(math.log(N)/W)+1
env = np.full(nb, -np.inf); np.maximum.at(env, j, val.astype(np.float64))
# intervals spanning several blocks: propagate to blocks covered (rare, small y only)
for i in np.nonzero(np.floor(np.log(nxt.astype(np.float64))/W).astype(np.int64) > j)[0]:
    j2 = min(int(math.floor(math.log(float(nxt[i]))/W)), nb-1)
    env[j[i]:j2+1] = np.maximum(env[j[i]:j2+1], float(val[i]))
np.save(sys.argv[2], env)
for lo in [1, 100, 1e4, 1e6, 1e7]:
    m = nl >= lo
    print(f"y>={lo:.0e}: max H/sqrt y {float(dev[m].max()):+.6f} at {float(nl[m][dev[m].argmax()]):.0f}; min {float(devmin[m].min()):+.6f}; env max {env[int(math.log(lo)/W):].max():+.6f}")
