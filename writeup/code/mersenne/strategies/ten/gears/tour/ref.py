"""Reference implementation in plain Python (no Montgomery, list arithmetic) of the wheel tests, to check tour.c."""
import sys, struct
def ring_test(n,r,a=1):
    def mul(p,q):
        t=[0]*r
        for i,x in enumerate(p):
            if x:
                for j,y in enumerate(q):
                    if y: t[(i+j)%r]=(t[(i+j)%r]+x*y)%n
        return t
    base=[0]*r; base[0]=a%n; base[1%r]=(base[1%r]+1)%n
    res=[0]*r; res[0]=1
    e=n; b=base
    while e:
        if e&1: res=mul(res,b)
        b=mul(b,b); e>>=1
    want=[0]*r; want[0]=a%n; want[n%r]=(want[n%r]+1)%n
    return res==want
def poly_test(n,c,a=1):
    d=len(c)
    def red(t):
        for k in range(len(t)-1,d-1,-1):
            v=t[k]
            if v:
                for j in range(d): t[k-d+j]=(t[k-d+j]+v*c[j])%n
        return t[:d]
    def mul(p,q):
        t=[0]*(2*d-1)
        for i,x in enumerate(p):
            if x:
                for j,y in enumerate(q):
                    if y: t[i+j]=(t[i+j]+x*y)%n
        return red(t)
    def pw(b):
        res=[0]*d; res[0]=1; e=n
        while e:
            if e&1: res=mul(res,b)
            b=mul(b,b); e>>=1
        return res
    if d==1: bx=[c[0]%n]; b1=[(c[0]+a)%n]
    else:
        bx=[0]*d; bx[1]=1; b1=[0]*d; b1[1]=1; b1[0]=a%n
    p1=pw(b1); p2=pw(bx); p2[0]=(p2[0]+a)%n
    return p1==p2
if __name__=='__main__':
    d=open('psp2_1e8.bin','rb').read(); L=struct.unpack('<%dQ'%(len(d)//8),d)
    for r in (5,7,12):
        mine=[n for n in L if ring_test(n,r)]
        print('ring',r,'liars on the list to 1e8:',len(mine),mine[:3])
    print('poly x^2=x+1 liars below 20000:',[n for n in range(9,20000,2) if poly_test(n,[1,1]) and any(n%p==0 for p in range(3,int(n**.5)+1,2))])
    print('poly x^3=x+1 on the list to 1e8:',[n for n in L if poly_test(n,[1,1,0])][:5])
    print('ring 5 direct below 300000 (no prefilter):',[n for n in range(9,300000,2) if n%5 and pow(2,n,n)==2 and ring_test(n,5) and any(n%p==0 for p in range(3,int(n**.5)+1,2))])
