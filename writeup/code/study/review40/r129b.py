import numpy as np
exec(open('/home/bas/scratch/review-logic5/r129.py').read().split("# check integral")[0])
x = Nsm(T)
ls = [1, 2, 4.5, 6, 8]
mid = len(T)//2
for name, idx in [('lower', slice(0, mid)), ('upper', slice(mid, None))]:
    TT = T[idx]; xx = x[idx]
    Lm = np.mean(np.log(TT/(2*np.pi)))
    out = []
    for l in ls:
        d, _ = sigma2(xx, l, 200000)
        out.append(f'l={l}: d={d:.4f} B={berry(l,Lm):.4f} diff={d-berry(l,Lm):+.4f}')
    print(name, f'Lm={Lm:.3f}', ' | '.join(out))

# block error bars: blocks of ~ 1000 unfolded units
for l in ls:
    edges = np.arange(x[0], x[-1], 1340.6)
    vals = []
    for a, b in zip(edges[:-1], edges[1:]):
        xx = x[(x >= a) & (x < b + l + 1)]
        y = np.linspace(a, b, 20000, endpoint=False)
        n = np.searchsorted(xx, y + l, side='right') - np.searchsorted(xx, y, side='left')
        vals.append(np.mean((n - l)**2))
    vals = np.array(vals)
    print(f'l={l}: blocks={len(vals)} mean={vals.mean():.4f} se={vals.std(ddof=1)/np.sqrt(len(vals)):.4f}')

# cutoff sensitivity
Lm = np.mean(np.log(T/(2*np.pi)))
for l in ls:
    print(l, ' '.join(f'cut{c}:{berry(l,Lm,c):.4f}' for c in [0.8, 0.9, 1.0, 1.1, 1.25]),
          f' L=ln t: {berry(l, Lm+np.log(2*np.pi)):.4f}')
# diagonal only up to tau=1, no tail; and GUE-form-factor smooth version for comparison
from scipy.integrate import quad
for l in ls:
    tail = 2/np.pi**2*quad(lambda u: np.sin(np.pi*l*u)**2/u**2, 1, 200, limit=5000)[0]
    print(l, f'tail share {tail/berry(l,Lm):.3f}')
