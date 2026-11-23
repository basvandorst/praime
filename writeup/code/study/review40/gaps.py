import numpy as np, glob
for f in sorted(glob.glob('/home/bas/praime/writeup/data/riemann/lzeros/r128/D*.npy')):
    D = int(f.split('/D')[-1][:-4]); q = abs(D)
    g = np.sort(np.load(f))
    sgn = 1 if D < 0 else -1
    x = g/(2*np.pi)*np.log(q*g/(2*np.pi*np.e)) + sgn/8.
    s = np.diff(x)
    big = np.where(s > 2.4)[0]; small = np.where(s < 0.02)[0]
    # compare count with N(T) at ends: x should be ~ index + offset
    res = x - np.arange(len(x)); 
    print(D, 'max sp %.3f min sp %.4f' % (s.max(), s.min()), 'big', [(round(g[i],2), round(s[i],2)) for i in big],
          'small', [(round(g[i],3), round(s[i],4)) for i in small], 'res start %.2f end %.2f' % (res[:20].mean(), res[-20:].mean()),
          'dup', np.sum(np.diff(g) < 1e-6))
