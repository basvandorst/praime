import numpy as np, sympy, time
X0, BLK, NB, H = 10**9, 10**8, 10, 10**6
small = np.array(list(sympy.primerange(3, int((X0+NB*BLK)**0.5)+2)), dtype=np.int64)

def count_odd_primes(a, b):
    # number of primes in [a, b), a odd-aligned handled generically
    n = b - a
    s = np.ones(n, dtype=bool)
    s[(a % 2)::2] = False  # even numbers off (a even -> index0 even)
    for p in small:
        st = (-a) % p
        if a + st == p: st += p
        s[st::p] = False
    return s

t0 = time.time()
counts = []
for j in range(NB):
    a = X0 + j*BLK
    c = []
    for k in range(0, BLK, 10**7):
        s = count_odd_primes(a+k, a+k+10**7)
        c.extend(s.reshape(-1, H).sum(1) if H <= 10**7 else [])
    counts.append(np.array(c))
print('sieve time', time.time()-t0, flush=True)
counts = np.array(counts)  # NB x (BLK/H)
meas = []
for j in range(NB):
    meas.append(counts[j].var(ddof=0)/counts[j].mean())
print('measured Var/mean (ddof0) %.4f' % np.mean(meas), ' ddof1 %.4f' % np.mean([counts[j].var(ddof=1)/counts[j].mean() for j in range(NB)]))
np.save('/home/bas/scratch/review-logic5/counts127.npy', counts)

g = np.loadtxt('/home/bas/praime/writeup/data/riemann/zeros_40000_flint.txt')
rho = 0.5 + 1j*g
def lisum(x):
    z = rho*np.log(x)
    ser = 1 + 1/z + 2/z**2 + 6/z**3 + 24/z**4
    return np.sum(np.exp(z)/z*ser)
zv = []
tails = []
for j in range(NB):
    a = X0 + j*BLK
    xs = a + H*np.arange(BLK//H + 1)
    S = np.array([lisum(float(x)) for x in xs])
    dev = -2*np.real(S[1:] - S[:-1])
    mid = xs[:-1] + H/2
    main = H/np.log(mid)
    pred = main + dev
    # diagonal random-phase tail above last zero for this block (x at block middle)
    xm = a + BLK/2; lx = np.log(xm); h = np.log1p(H/xm)
    gg = np.linspace(g[-1], 2e7, 4000001)
    f = 8*xm*np.sin(gg*h/2)**2/(gg**2*lx**2)*np.log(gg/(2*np.pi))/(2*np.pi)
    tail = np.trapezoid(f, gg) + 8*xm/(2*lx**2)*np.log(2e7/(2*np.pi))/(2*np.pi)/2e7  # sin^2 avg 1/2 beyond
    zv.append((pred.var() + tail)/pred.mean()); tails.append(tail/pred.mean())
    if j == 0:
        print('block0 corr(pred,counts)=%.3f' % np.corrcoef(pred, counts[0])[0, 1])
print('zeros Var/mean %.4f  tail share %.3f' % (np.mean(zv), np.mean(tails)/np.mean(zv)))
