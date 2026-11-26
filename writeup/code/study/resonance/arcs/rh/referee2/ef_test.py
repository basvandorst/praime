# Own test of Lemma EF for L(s, lambda^m), K = Q(i), q = 1, k = 4m (zeros from PARI; GRH used only numerically).
# Psi = sum_a Lambda(a) xi(a) V(Na) should equal -sum_rho Vt(rho) - L'/L(0) + (1/2pi) int Vt(-1/4+it)(-L'/L)(-1/4+it) dt,
# with -L'/L(-1/4+it) = 2 log A + psi(-1/4+it+|k|/2) + psi(5/4-it+|k|/2) + L'/L(5/4-it)  (FE; L real-coefficient, A = 1/pi).
import cypari2, math, cmath, sys, time
import mpmath as mp
pari = cypari2.Pari(); pari.allocatemem(800*10**6); pari.set_real_precision(25)
m = int(sys.argv[1]); x = float(sys.argv[2]); k = 4*m
H = math.log(x)**2; Tz = 5.0*H
t0 = time.time()
pari('bnf=bnfinit(y^2+1,1); gc=gcharinit(bnf,1); L=lfuncreate([gc,[%d,0]~])' % m)
a5 = complex(pari('lfunan(L,5)[5]'))
print(f'k={k} x={x} H={H:.1f}; a_5 PARI {a5.real:.6f} expected {2*(((2+1j)/math.sqrt(5))**k).real:.6f}')
def V(u): return 0.5*math.erfc(H*math.log(u/x)/math.sqrt(2))
umax = x*math.exp(9/H); R = math.isqrt(int(umax))+2
def fac(n):
    f = {}; d = 2
    while d*d <= n:
        while n % d == 0: f[d] = f.get(d,0)+1; n//=d
        d += 1
    if n > 1: f[n] = f.get(n,0)+1
    return f
psi = 0j
for a in range(1, R):
    for b in range(0, R):
        n = a*a+b*b
        if n > umax: break
        f = fac(n)
        if len(f) != 1: continue
        p, e = next(iter(f.items()))
        if p == 2: lam = math.log(2)
        elif p % 4 == 3:
            if e % 2: continue
            lam = 2*math.log(p)
        else:
            if math.gcd(a,b) != 1: continue
            lam = math.log(p)
        z = complex(a,b)/math.sqrt(n)
        psi += lam*(z**k)*V(n)
print(f'LHS done {time.time()-t0:.1f}s')
pari('Ld=lfuninit(L,[1/2,3/4,%d])' % int(Tz+5))
print(f'lfuninit done {time.time()-t0:.1f}s')
zs = [float(t) for t in pari('lfunzeros(Ld,[%f,%f])' % (-Tz, Tz))]
print(f'zeros done {time.time()-t0:.1f}s')
def Vt(s): return cmath.exp(s*math.log(x) + s*s/(2*H*H))/s
zsum = sum(Vt(0.5+1j*g) for g in zs)
lg = lambda t: float(mp.im(mp.loggamma(0.5+abs(k)/2+1j*t)))
th = lambda t: t*math.log(1/math.pi) + lg(t)
print(f'zeros found {len(zs)} in [-{Tz:.0f},{Tz:.0f}], smooth count {(th(Tz)-th(-Tz))/math.pi:.1f}')
res0 = -complex(pari('lfun(L,0,1)/lfun(L,0)'))
h = 0.1; N = int(2*Tz/h); line = 0j; logA = math.log(1/math.pi)
for j in range(N+1):
    t = -Tz + j*h
    s = complex(-0.25, t)
    ld = complex(pari('my(s=5/4-I*(%r)); lfun(Ld,s,1)/lfun(Ld,s)' % t))
    mlog = 2*logA + complex(mp.digamma(s+abs(k)/2)) + complex(mp.digamma(1-s+abs(k)/2)) + ld
    w = h*(0.5 if j in (0,N) else 1.0)
    line += w*Vt(s)*mlog
line /= 2*math.pi
rhs = -zsum + res0 + line
print(f'Psi = {psi:.6f}')
print(f'-sum Vt(rho) = {-zsum:.6f}   (sqrt x = {math.sqrt(x):.2f})')
print(f'res0 = {res0:.6f}  line = {line:.6f}')
print(f'RHS = {rhs:.6f}   |LHS-RHS| = {abs(psi-rhs):.2e}')
ld1 = abs(complex(pari("lfun(Ld,1,1)/lfun(Ld,1)")))
print(f'|R| = |Psi + sum Vt| = {abs(psi+zsum):.4f}; EF bound without the C x^(-1/4) term: '
      f'{math.log(x)+2*math.log(math.pi)+2*math.log(abs(k)+4)+20+ld1:.2f}')
print(f'time {time.time()-t0:.1f}s')
