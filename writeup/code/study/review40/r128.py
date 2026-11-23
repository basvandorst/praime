import numpy as np, glob
exec(open('/home/bas/scratch/review-logic5/r129.py').read().split("# check integral")[0])
xz = Nsm(T)
Lz = np.mean(np.log(T/(2*np.pi)))
sz, _ = sigma2(xz, 1.0)
print('zeta S2(1)=%.4f' % sz)

def s2_blocks(x, l, bs):
    # block estimate: split by zero count into blocks of bs zeros
    vals = []
    for i in range(0, len(x) - bs + 1, bs):
        xx = x[i:i+bs+5]
        a, b = x[i], x[min(i+bs, len(x)-1)]
        if b - a < bs*0.8: continue
        y = np.linspace(a, b - (l if i+bs >= len(x)-1 else 0), 20000, endpoint=False)
        n = np.searchsorted(xx, y+l, side='right') - np.searchsorted(xx, y, side='left')
        vals.append(np.mean((n-l)**2))
    return np.array(vals)

rows = []
for f in sorted(glob.glob('/home/bas/praime/writeup/data/riemann/lzeros/r128/D*.npy')):
    D = int(f.split('/D')[-1][:-4]); q = abs(D)
    g = np.sort(np.load(f))
    Lg = np.log(q*g/(2*np.pi))
    sgn = 1 if D < 0 else -1   # odd char: +1/8, even: -1/8
    x = g/(2*np.pi)*np.log(q*g/(2*np.pi*np.e)) + sgn/8.
    s, n = sigma2(x, 1.0, 200000)
    Lm = Lg.mean()
    pred = -2/np.pi**2*sum(np.sin(np.pi*k*np.log(q)/Lm)**2/(k**2*q**k) for k in range(1, 6))
    # berry-full prediction difference using same L
    bv = s2_blocks(x, 1.0, 200)
    se = bv.std(ddof=1)/np.sqrt(len(bv))
    # mean-count check
    rows.append((D, len(g), Lg.min(), Lg.max(), s, s - sz, pred, se, n.mean()-1, x[-1]-x[0]-len(g)+1))
    print('D=%4d n=%4d L=[%.3f,%.3f] S2=%.4f diff=%+.4f pred=%+.4f se_blk=%.4f meanN-1=%+.4f span-n=%+.2f' % rows[-1])
r = np.array(rows)
d, p, se = r[:, 5], r[:, 6], r[:, 7]
print('mean diff %.4f, mean pred %.4f' % (d.mean(), p.mean()))
print('se from scatter of diffs: %.4f; se of (diff-pred) scatter: %.4f; se from block se: %.4f' %
      (d.std(ddof=1)/np.sqrt(len(d)), (d-p).std(ddof=1)/np.sqrt(len(d)), np.sqrt(np.sum(se**2))/len(d)))
# expected sampling sd of S2(1) for n ~ 600 zeros using zeta blocks of the same size
for bs in [420, 600, 969]:
    v = []
    for i in range(0, len(xz)-bs, bs):
        xx = xz[i:i+bs]
        vv, _ = sigma2(xx, 1.0, 20000); v.append(vv)
    print('zeta chunks of %d zeros: sd of S2(1) = %.4f (n=%d)' % (bs, np.std(v, ddof=1), len(v)))
