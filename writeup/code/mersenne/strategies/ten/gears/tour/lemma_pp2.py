"""Prime squares on Pascal's wheel, by roots of unity.
p h(x) = (x+1)^p - x^p - 1.  Proved in the text: p^2 comes home on the ring of r cells  <=>  h = 0 mod (p, x^r - 1).
Here: (1) h(1) = (2^p - 2)/p, so ring 1 needs a Wieferich prime; (2) h is divisible by x(x+1)(x^2+x+1) as a polynomial with
integer coefficients, so the rings of 2 and 3 cells add nothing; (3) the ring of 4 cells adds nothing either for the two known
Wieferich primes; (4) which rings of up to 200 cells the squares of 1093 and 3511 fool."""
from sympy import symbols, Poly, primerange, GF
x=symbols('x')
ok=True
for p in primerange(5,60):
    h=Poly(((x+1)**p-x**p-1)/p,x,domain='ZZ')
    q,rem=divmod(h,Poly(x*(x+1)*(x**2+x+1),x,domain='ZZ'))
    ok&=rem.is_zero and h.eval(1)==(2**p-2)//p
print("(1),(2) for the primes 5..59: h(1) = (2^p - 2)/p and x(x+1)(x^2+x+1) divides h:",ok)
def h_mod(p,r):
    M=p*p
    def mul(a,b):
        t=[0]*r
        for i,u in enumerate(a):
            if u:
                for j,v in enumerate(b):
                    if v: t[(i+j)%r]=(t[(i+j)%r]+u*v)%M
        return t
    base=[0]*r; base[0]=1; base[1%r]=(base[1%r]+1)%M
    res=[0]*r; res[0]=1; e=p; b=base
    while e:
        if e&1: res=mul(res,b)
        b=mul(b,b); e>>=1
    res[p%r]=(res[p%r]-1)%M; res[0]=(res[0]-1)%M
    return [v//p for v in res] if all(v%p==0 for v in res) else None
for p in (1093,3511):
    z=[r for r in range(1,201) if p%r and (lambda v:v is not None and all(c%p==0 for c in v))(h_mod(p,r))]
    print("(3),(4) p = %d (2^(p-1) = 1 mod p^2: %s): h vanishes mod (p, x^r - 1) for r = %s among r <= 200"%(p,pow(2,p-1,p*p)==1,z))
w=[p for p in primerange(3,10**6) if pow(2,p-1,p*p)==1]
print("Wieferich primes below 1e6:",w)
