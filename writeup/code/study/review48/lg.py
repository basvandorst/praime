# Landau-Gonek on the stored zeros: mean over zeros of p^{-i gamma} vs -chi(p) ln p/(sqrt p L); bad primes: 0
import numpy as np, json, sympy
z=np.load('../all.npz'); d=[r for r in json.load(open('../perL.json')) if r['q']!=1153]
for p in [2,3,5,7]:
    acc={'bad':[], 'good':[]}
    for r in d:
        g=z['%s_%d_g'%(r['set'],r['D'])]; L=np.log(r['q']*g/(2*np.pi))
        c=sympy.jacobi_symbol(r['D']%p,p) if p>2 else None
        chi=int(sympy.ntheory.residue_ntheory.legendre_symbol(r['D']%p,p)) if (p>2 and r['q']%p) else (0 if r['q']%p==0 else (1 if r['D']%8 in (1,7) else -1))
        e=np.exp(-1j*g*np.log(p))
        if chi==0: acc['bad'].append(np.mean(e))
        else: acc['good'].append(np.mean(e*chi*L)/( -np.log(p)/np.sqrt(p)))   # should be ~1 (Landau-Gonek, first order)
    b=np.array(acc['bad']); gd=np.array(acc['good'])
    print('p=%d bad (n=%d): mean p^{-i gamma} = %.5f%+.5fi +- %.5f ; good (n=%d): chi*L*mean/(-ln p/sqrt p) = %.4f +- %.4f'%(
        p,len(b),b.mean().real,b.mean().imag,b.real.std()/np.sqrt(len(b)),len(gd),gd.mean().real,gd.real.std()/np.sqrt(len(gd))))
