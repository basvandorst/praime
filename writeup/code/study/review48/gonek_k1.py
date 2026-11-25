# Own main term of S(T) = sum_{0<gamma<=T} L'(rho,chi)L'(1-rho,chi), chi real primitive mod q,
# by Gonek's method (residue theorem + functional equation + Gonek's lemma):
#   S(T) = -(1/phi(q)) Res_{s=1} [ (A'' - 3 B' + 2 C)(s) Y^s / s ],  Y = qT/2pi,
#   A = zq*zq', B = zq'^2, C = zq'^3/zq, zq = zeta(s) prod_{p|q}(1-p^-s).
# Truncated Laurent series in u = s-1 with mpmath coefficients.
import mpmath as mp, sympy
mp.mp.dps=30
K=12
class LS:
    def __init__(s,v,a): s.v=v; s.a=[mp.mpf(x) for x in a][:K]+[mp.mpf(0)]*max(0,K-len(a))
    def __mul__(s,o):
        c=[mp.mpf(0)]*K
        for i in range(K):
            if s.a[i]==0: continue
            for j in range(K-i): c[i+j]+=s.a[i]*o.a[j]
        return LS(s.v+o.v,c)
    def __add__(s,o):
        v=min(s.v,o.v); c=[mp.mpf(0)]*K
        for i in range(K):
            if i+s.v-v<K: c[i+s.v-v]+=s.a[i]
            if i+o.v-v<K: c[i+o.v-v]+=o.a[i]
        return LS(v,c)
    def scale(s,x): return LS(s.v,[x*t for t in s.a])
    def deriv(s):
        # d/du of sum a_k u^(v+k)
        c=[(s.v+k)*s.a[k] for k in range(K)]
        # new valuation v-1 (if v=0 the k=0 term vanishes: fine, stays as 0 coefficient)
        return LS(s.v-1,c)
    def inv(s):
        a0=s.a[0]; b=[mp.mpf(0)]*K; b[0]=1/a0
        for n in range(1,K):
            b[n]=-sum(s.a[k]*b[n-k] for k in range(1,n+1))/a0
        return LS(-s.v,b)
    def coef(s,n):  # coefficient of u^n
        k=n-s.v
        return s.a[k] if 0<=k<K else mp.mpf(0)
def zeta_ls():
    a=[mp.mpf(1)]+[(-1)**n*mp.stieltjes(n)/mp.factorial(n) for n in range(K-1)]
    return LS(-1,a)
def E_ls(q):
    e=LS(0,[1])
    for p in sympy.factorint(q):
        lp=mp.log(p); t=[ (1 if k==0 else 0) - mp.mpf(p)**-1*(-lp)**k/mp.factorial(k) for k in range(K)]
        e=e*LS(0,t)
    return e
Z1=zeta_ls()
def poly_coeffs(q):
    """returns coefficients P[m], m=0..4, with S(T) = (T/2pi) sum_m P[m] x^m, x = ln(qT/2pi)"""
    z=Z1*E_ls(q); zp=z.deriv()
    A=z*zp; B=zp*zp; C=zp*zp*zp*z.inv()
    G=A.deriv().deriv()+B.deriv().scale(-3)+C.scale(2)
    # Res_u Y e^{ux}/(1+u) G(1+u):  e^{ux}/(1+u) = sum_j u^j sum_{i<=j} x^i/i! (-1)^{j-i}
    P=[mp.mpf(0)]*5
    for j in range(0,5):
        g=G.coef(-1-j)
        for i in range(j+1):
            P[i]+=g*(-1)**(j-i)/mp.factorial(i)
    cq=mp.mpf(1)
    for p in sympy.factorint(q): cq*=1-mp.mpf(1)/p
    return [-x/cq for x in P]     # -(1/phi(q)) * Y = -(T/2pi)/c_q
def S(q,T):
    P=poly_coeffs(q); x=mp.log(q*T/(2*mp.pi))
    return T/(2*mp.pi)*sum(P[m]*x**m for m in range(5))
if __name__=='__main__':
    P=poly_coeffs(1)
    g=[mp.stieltjes(n) for n in range(4)]; g0,g1,g2,g3=g
    mil=[ (6+6*g0*(5*g1+4*g2-2)+6*g0**2*(g0+g0**2+6*g1+1)-12*g1+42*g1**2+3*g2+10*g3)/3,
          -2+4*g0-2*g0**2-2*g0**3-10*g0*g1+4*g1-g2, 1-2*g0+g0**2-2*g1, (2*g0-1)/3, mp.mpf(1)/12]
    print('q=1 own   :',[mp.nstr(x,12) for x in P])
    print('Milinovich:',[mp.nstr(x,12) for x in mil])
