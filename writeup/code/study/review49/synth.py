# Expected finite-band cross terms: var(dS_X) - (orthogonal-line sum on the same windows) for Sato-Tate random a_p
import sys, json, math, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import e2
X = int(sys.argv[1]); reps = int(sys.argv[2]); cid = sys.argv[3]
c = [x for x in e2.FS if x['id'] == cid][0]; N = c['N']
g = np.unique(np.sort(e2.Z[cid])); g = g[g > 1e-6]; b = g[(g >= 20) & (g <= 400)]; ub = e2.theta(N, b) / math.pi
u0, u1 = ub[0], ub[-1]; xg = u0 + 0.01 * np.arange(int((u1 - u0) / 0.01) + 2)
tt = np.linspace(19, 401, 400001); uu = e2.theta(N, tt) / math.pi; tg = np.interp(xg, uu, tt)
import cypari2
pari = cypari2.Pari()
pr = np.array([int(p) for p in pari.primes(pari.primepi(X))]); bad = np.array([N % p == 0 for p in pr])
rng = np.random.default_rng(5)
# Sato-Tate sampler by rejection
def st(n):
    out = []
    while len(out) < n:
        th = rng.uniform(0, math.pi, 4 * n); keep = rng.uniform(0, 1, 4 * n) < np.sin(th) ** 2; out.extend(th[keep])
    return np.array(out[:n])
LS = e2.LS; acc = []; accD = []; MW = None
for r in range(reps):
    th = st(len(pr)); terms = []
    for p, t_, bd in zip(pr, th, bad):
        k = 1
        while p ** k <= X:
            bk = ((-1) ** k) * p ** (-k / 2.0) if bd else 2 * math.cos(k * t_)
            terms.append((-bk / (math.pi * k * p ** (k / 2)), k * math.log(p))); k += 1
    T = np.array(terms); SX = np.zeros_like(xg)
    if MW is None:
        MW = []
        for l in LS:
            m = int(round(l / 0.01)); nw = len(np.arange(u0, u1 - l, 0.01)); w = tg[m:m + nw] - tg[:nw]
            MW.append(np.array([(2 * np.sin(w * f / 2) ** 2).mean() for f in T[:, 1]]))
    for a in range(0, len(T), 64):
        SX += (T[a:a + 64, 0][None, :] * np.sin(np.outer(tg, T[a:a + 64, 1]))).sum(1)
    v = []; d = []
    for l in LS:
        m = int(round(l / 0.01)); nw = len(np.arange(u0, u1 - l, 0.01))
        dd = SX[m:m + nw] - SX[:nw]; v.append(dd.var())
        d.append(float((T[:, 0] ** 2 * MW[len(d)]).sum()))
    acc.append(v); accD.append(d)
acc = np.array(acc); accD = np.array(accD); diff = acc - accD
pk = [0, 1, 3, 7, 11, 15, 19]
print('X', X, 'reps', reps, cid)
print('mean[var(dS_X) - orthogonal sum]', ' '.join('%+.5f' % x for x in diff.mean(0)[pk]))
print('se                              ', ' '.join('%.5f' % x for x in (diff.std(0, ddof=1) / math.sqrt(reps))[pk]))
print('sd per draw                     ', ' '.join('%.4f' % x for x in diff.std(0, ddof=1)[pk]))
