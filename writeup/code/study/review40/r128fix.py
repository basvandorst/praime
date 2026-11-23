import numpy as np, glob
exec(open('/home/bas/scratch/review-logic5/r129.py').read().split("# check integral")[0])
sz,_ = sigma2(Nsm(T), 1.0)
add = {-47: [482.5284, 482.5419], -59: [344.3561, 344.3862]}
res = {}
for fix in [False, True]:
    d = []; p = []
    for f in sorted(glob.glob('/home/bas/praime/writeup/data/riemann/lzeros/r128/D*.npy')):
        D = int(f.split('/D')[-1][:-4]); q = abs(D)
        g = np.load(f)
        if fix and D in add: g = np.concatenate([g, add[D]])
        g = np.sort(g)
        x = g/(2*np.pi)*np.log(q*g/(2*np.pi*np.e)) + (1 if D<0 else -1)/8.
        s,_ = sigma2(x, 1.0, 400000)
        Lm = np.log(q*g/(2*np.pi)).mean()
        pr = -2/np.pi**2*sum(np.sin(np.pi*k*np.log(q)/Lm)**2/(k**2*q**k) for k in range(1,6))
        d.append(s-sz); p.append(pr)
        if D in add: print('fix' if fix else 'raw', D, 'S2=%.4f diff=%+.4f pred=%+.4f'%(s, s-sz, pr))
    d=np.array(d); p=np.array(p)
    print('fix' if fix else 'raw', 'mean diff %+.4f  se(scatter) %.4f  mean pred %+.4f  chi2(diff-pred, sd_i from scatter)' % (d.mean(), d.std(ddof=1)/np.sqrt(10), p.mean()))
