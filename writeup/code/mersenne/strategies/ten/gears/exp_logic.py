import numpy as np, math, sys
exec(open('machines.py').read().split("RULES=[")[0])
RULES=[('div','Wheels'),('square','Square jumps'),('double','Doubling'),('fib','Fibonacci'),('accel','Accelerating'),('flip','Back and forth'),('pclock','Prime clock'),('slip','Slipping')]
for rule,name in RULES:
    H=np.zeros((NMAX+1,T+1),dtype=bool)
    for k in range(1,NMAX+1): H[k]=home(rule,k)
    H[:,0]=False
    if rule=='pclock': H[:,PI==0]=False
    line=[]
    for N in NS:
        ks=np.arange(2,N+1); pk=[k for k in ks if isp[k]]; ck=[k for k in ks if not isp[k]]; ek=[k for k in EXPS if k<=N]; nek=[k for k in ks if k not in ek]
        chh=H[ck].sum(0); eh=H[ek].sum(0); neh=H[nek].sum(0)
        new=np.where(chh==0,eh,0); old=np.where(neh==0,eh,0); anyc=eh
        tn=int(np.argmax(new)); to=int(np.argmax(old)); ta=int(np.argmax(anyc))
        line.append("N=%3d: only exponents %d (tick %d) | other primes allowed %d %s (tick %d) | anything allowed %d (tick %d) of %d"%(N,old.max(),to,new.max(),[int(k) for k in pk if H[k][tn]],tn,anyc.max(),ta,len(ek)))
    print(name); [print("   "+l) for l in line]
