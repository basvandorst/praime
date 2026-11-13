"""Find zeros missed by the grid of lfam2.py for the characters whose argument-principle count disagrees.
A close pair of zeros inside one grid step gives no sign change but a local minimum of |Z| between neighbours of
the same sign (between two ordinary zeros |Z| has no interior minimum). For each such dip: evaluate Z at the vertex
of the parabola through the three grid values; if the sign flips there, locate both zeros by bisection.
Usage: lfam_fix.py zeros.npz out.npz"""
import sys, math, time
import numpy as np
from scipy.special import loggamma, bernoulli
zf, outf = sys.argv[1], sys.argv[2]
D = dict(np.load(zf)); q = int(D['q']); T = float(D['T']); h = float(D['h']); pw = D['pw']
mism = D['mism']; ks = np.array(sorted(set(mism.tolist())))
NT_, M_ = 30, 9
alpha = pw / q
LOG = np.log(np.arange(NT_)[None, :] + alpha[:, None]); Na = NT_ + alpha; LNa = np.log(Na)
B = bernoulli(2 * M_); coef = [float(B[2 * j]) / math.factorial(2 * j) for j in range(1, M_ + 1)]
def hurwitz_all(s):
    main = np.exp(-s * LOG).sum(axis=1); Nas = np.exp(-s * LNa)
    tail = Nas * Na / (s - 1) + 0.5 * Nas; rising = s; powr = Nas / Na
    for j in range(1, M_ + 1):
        tail = tail + coef[j - 1] * rising * powr; rising = rising * (s + 2 * j - 1) * (s + 2 * j); powr = powr / (Na * Na)
    return main + tail
jj = np.arange(q - 1)
CH = np.exp(2j * math.pi * np.outer(ks, jj) / (q - 1))       # chi_k(g^j) for the characters to repair
tau = CH @ np.exp(2j * math.pi * pw / q)
kap = ks % 2
eps = tau / ((1j) ** kap * math.sqrt(q)); rot = np.exp(-0.5j * np.angle(eps))
def Zk(t):
    s = 0.5 + 1j * t
    L = np.exp(-s * math.log(q)) * (CH @ hurwitz_all(s))
    th = 0.5 * t * math.log(q / math.pi) + np.imag(loggamma((0.5 + kap + 1j * t) / 2))
    return np.real(rot * np.exp(1j * th) * L)
t0 = time.time()
pad = 8
ts = (np.arange(int(math.ceil(T / h)) + 2 * pad + 1) - pad) * h
Zg = np.array([Zk(t) for t in ts])                            # (rows, len(ks))
print(f'grid for {len(ks)} characters: {time.time() - t0:.0f} s', flush=True)
new_k, new_g = [], []
for c_, k in enumerate(ks):
    z = Zg[:, c_]
    for i in range(1, len(ts) - 1):
        if z[i - 1] * z[i] > 0 and z[i] * z[i + 1] > 0 and abs(z[i]) < abs(z[i - 1]) and abs(z[i]) < abs(z[i + 1]) and 0 < ts[i] <= T + h:
            # parabola vertex
            den = z[i - 1] - 2 * z[i] + z[i + 1]
            tv = ts[i] + 0.5 * h * (z[i - 1] - z[i + 1]) / den
            zv = Zk(tv)[c_]
            if zv * z[i] < 0:
                found = []
                for (lo, hi) in ((ts[i - 1] if tv < ts[i] else ts[i], tv), (tv, ts[i + 1] if tv > ts[i] else ts[i])):
                    if lo > hi: lo, hi = hi, lo
                    zlo = Zk(lo)[c_]
                    if zlo * Zk(hi)[c_] > 0: continue
                    for _ in range(45):
                        mid = 0.5 * (lo + hi); zm = Zk(mid)[c_]
                        if zm * zlo > 0: lo, zlo = mid, zm
                        else: hi = mid
                    found.append(0.5 * (lo + hi))
                for gmm in found:
                    if 0 < gmm <= T:
                        new_k.append(int(k)); new_g.append(gmm)
                print(f'chi_{k}: dip at t = {ts[i]:.4f} (|Z| = {abs(z[i]):.2e}), pair found at {[round(v, 10) for v in found]}', flush=True)
D['k'] = np.concatenate([D['k'], np.array(new_k, dtype=D['k'].dtype)])
D['gamma'] = np.concatenate([D['gamma'], np.array(new_g)])
D['interr'] = np.concatenate([D['interr'], np.full(len(new_g), 1e-10, dtype=D['interr'].dtype)])
kbar = (-np.arange(q - 1)) % (q - 1)
found = np.bincount(D['k'], minlength=q - 1)
count = found + found[kbar]
D['count'] = count
still = [int(k) for k in range(1, q - 1) if int(round(D['NAP'][k])) != int(count[k])]
D['mism'] = np.array(still, dtype=int)
print(f'added {len(new_g)} zeros; characters whose count still differs: {still}; {time.time() - t0:.0f} s')
np.savez(outf, **D)
