"""c_E = B/L(1,Sym^2 E) for small curves given by a-invariants; |K| = (2 sqrt2/3) c N^(-1/4); root numbers."""
import cypari2, math, json
P = cypari2.Pari(); P.allocatemem(300*10**6); P.set_real_precision(30)
cur = {'11a1':[0,-1,1,-10,-20],'14a1':[1,0,1,4,-6],'15a1':[1,1,1,-10,-10],'17a1':[1,-1,1,-1,-14],'19a1':[0,1,1,-9,-15],
       '21a1':[1,0,0,-4,-1],'26a1':[1,0,1,-5,-8],'26b1':[1,-1,1,-3,3],'30a1':[1,0,1,1,2],'33a1':[1,1,0,-11,0],'35a1':[0,1,1,9,1],
       '37a1':[0,0,1,-1,0],'37b1':[0,1,1,-23,-50],'43a1':[0,1,1,0,0],'53a1':[1,-1,1,0,0],'57a1':[0,-1,1,-2,2],'58a1':[1,-1,0,-1,1],
       '61a1':[1,0,0,-2,1],'65a1':[1,0,0,-1,0],'67a1':[0,1,1,-12,-21],'79a1':[1,1,1,-2,0],'83a1':[1,1,1,1,0],'89a1':[1,1,1,-1,0]}
out = {}
for nm, ai in cur.items():
    E = P.ellinit(ai); N = int(P.ellglobalred(E)[0]); eps = int(P.ellrootno(E))
    Ls = P.lfunsympow(E, 2); L1 = float(P.real(P.lfun(Ls, 2)))
    bad = [int(x) for x in P.factor(N)[0]]
    B = 1.0
    for p in P.primes(30000):
        p = int(p); X = 1.0/p; ap = int(P.ellap(E, p))
        if p in bad:
            cp = p/(p-1.0) if ap == 0 else 1 + 1.0/p
            loc = float(P('my(E=ellinit(%s),f=lfuneuler(lfunsympow(E,2),%d)); subst(f,variable(f),1/%d^2)' % (ai, p, p)))
            B *= cp * loc
        else:
            lam2 = ap*ap / p
            B *= (1 - (lam2 - 1)*X) / ((1 - X)*(1 - (lam2 - 2)*X + X*X))
    c = B / L1
    K = (2*math.sqrt(2)/3) * c * N**-0.25
    aps = [int(P.ellap(E, p)) for p in (2,3,5,7)]
    out[nm] = dict(ai=ai, N=N, eps=eps, c=c, absK=K, bad=bad, ap2357=aps, L1sym2=L1)
    print('%-6s N=%3d eps=%+d c=%.5f |K|=%.4f bad=%s a2357=%s' % (nm, N, eps, c, K, bad, aps), flush=True)
json.dump(out, open('curves_c.json','w'), indent=1)
