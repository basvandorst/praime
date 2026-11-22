import numpy as np, math
def base_primes(n):
    s=np.ones(n+1,bool); s[:2]=False
    for i in range(2,int(n**0.5)+1):
        if s[i]: s[i*i::i]=False
    return np.flatnonzero(s)
_BP={}
def primes_in(lo, hi):
    """primes in [lo, hi), lo odd-aligned handled; hi-lo <= 1e8"""
    r=int(math.isqrt(hi))+1
    if r not in _BP: _BP.clear(); _BP[r]=base_primes(r)[1:]  # odd primes
    bp=_BP[r]
    if lo%2==0: lo0=lo+1
    else: lo0=lo
    n=(hi-lo0+1)//2  # odd numbers lo0, lo0+2, ...
    m=np.ones(n,bool)
    for p in bp.tolist():
        if p*p>=hi: break
        st=max(p*p, ((lo0+p-1)//p)*p)
        if st%2==0: st+=p
        m[(st-lo0)//2::p]=False
    P=lo0+2*np.flatnonzero(m).astype(np.int64)
    if lo<=2<hi: P=np.concatenate([[2],P])
    return P[P<hi]
if __name__=='__main__':
    import time, gmpy2
    t=time.time(); P=primes_in(10**12,10**12+10**8); print(len(P), time.time()-t)
    # check against gmpy2 on a piece
    x=10**12; Q=[]
    while True:
        x=int(gmpy2.next_prime(x))
        if x>=10**12+10**6: break
        Q.append(x)
    print(np.array_equal(P[P<10**12+10**6],np.array(Q)))
    t=time.time(); P=primes_in(10**13,10**13+10**8); print(len(P), time.time()-t)
    x=10**13+10**8-10**6; Q=[]
    while True:
        x=int(gmpy2.next_prime(x))
        if x>=10**13+10**8: break
        Q.append(x)
    print(np.array_equal(P[P>10**13+10**8-10**6],np.array(Q)))
