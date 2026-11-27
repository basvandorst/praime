# Independent exact check of a certificate tree written by bb3.py.
# Plain Python integers only: no solver, no numpy, no floating point.  Shares no code with bb3.py.
#
# CLAIM (level n, L = 2^n*45, B = all divisors m >= MMIN of L, default MMIN = 6):
#   let r assign to every m in B a residue class r_m mod m, and let y >= 0 assign a real number to every
#   residue class modulo 2^n*d, d | 45.  If every b mod L lies in a class r_m or has
#   sum_d y(b mod 2^n d) >= 1, then  sum_r y(r mod 2^n d) >= 1  for at least one d.
#
# TREE.  A node carries restrictions  r_m = rho_m (mod g_m), g_m | m  (root: g_m = 1 for all m).
#   internal {"m":m,"q":q,"ch":[c_0..c_{q-1}]}: g_m*q | m; child j has g_m -> g_m*q, rho_m -> rho_m + g_m*j.
#       The q children exhaust the node.
#   child {"sym":[j2,perm]}: perm is a permutation of Z/q^e (q^e = exact power of q in L).  It defines
#       phi on Z/L (act on b mod q^e, keep b mod L/q^e).  Checked here: phi is a bijection of Z/L that maps every
#       residue class modulo every divisor of L onto a residue class of the same modulus; phi maps every
#       restriction set {b = rho_m mod g_m} onto itself and the set of child j onto the set of child j2;
#       child j2 is a real subtree.  Then phi carries every (r,y) of child j to one of child j2 with the same sums.
#   leaf {"w":{b:W_b}}: integers W_b >= 0.  With U = residues outside the fixed classes (g_m = m):
#       N = W(U) - sum_{m not fixed} max_{r = rho_m mod g_m} W(U & r mod m),
#       D = sum_{d|45} max_r W(U & r mod 2^n d);   checked: D > 0 and N >= D.
#       (Weak duality: W(U) <= sum_{m not fixed} W(U & r_m) + sum_d sum_r y_{d,r} W(U & class) <= (W(U)-N) + t*D.)
import sys, json
cert=json.load(open(sys.argv[1])); n=cert['n']; L=2**n*45
MMIN=int(sys.argv[2]) if len(sys.argv)>2 else 6
PERTURB=len(sys.argv)>3          # mutation test: lower one weight in every leaf; the check must then fail somewhere
DIVS=[m for m in range(1,L+1) if L%m==0]
BASE=[m for m in DIVS if m>=MMIN]
TOP=[2**n*d for d in (1,3,5,9,15,45)]
QE={2:2**n,3:9,5:5}
cnt={'leaves':0,'internal':0,'sym':0,'maxW':0,'tight':0,'maxdepth':0}
phi_ok={}
def make_phi(q,perm):
    key=(q,tuple(perm))
    if key in phi_ok: return phi_ok[key]
    N=QE[q]; rest=L//N
    assert sorted(perm)==list(range(N))
    phi=[None]*L
    for b in range(L):
        target=perm[b%N]
        # unique b2 with b2 = target (mod N) and b2 = b (mod rest)
        b2=[c for c in range(b%rest,L,rest) if c%N==target]
        assert len(b2)==1
        phi[b]=b2[0]
    assert sorted(phi)==list(range(L))
    for m in DIVS:                       # classes mod m go onto classes mod m
        for r in range(m):
            img={phi[b]%m for b in range(r,L,m)}
            assert len(img)==1, ('not class preserving',q,perm,m,r)
    phi_ok[key]=phi
    return phi
def check(node,st,depth):
    cnt['maxdepth']=max(cnt['maxdepth'],depth)
    if 'w' in node:
        U=[b for b in range(L) if not any(g==m and b%m==rho%m for m,(g,rho) in st.items())]
        W={int(b):v for b,v in node['w'].items()}
        assert all(type(v) is int and v>=0 for v in W.values())
        if PERTURB:
            k=max(W,key=lambda b:(W[b],b)); W[k]-=1
        def cls(m,r): return sum(W.get(b,0) for b in U if b%m==r)
        Nn=sum(W.get(b,0) for b in U)
        for m,(g,rho) in st.items():
            if g==m: continue
            Nn-=max(cls(m,r) for r in range(m) if r%g==rho%g)
        Dd=sum(max(cls(m,r) for r in range(m)) for m in TOP)
        if not (Dd>0 and Nn>=Dd):
            raise AssertionError(('leaf fails',st,Nn,Dd))
        cnt['leaves']+=1; cnt['maxW']=max(cnt['maxW'],max(W.values())); cnt['tight']+=(Nn==Dd)
        return
    m=node['m']; q=node['q']; ch=node['ch']
    assert m in st and q in (2,3,5)
    g,rho=st[m]
    assert m%(g*q)==0 and len(ch)==q
    cnt['internal']+=1
    for j,c in enumerate(ch):
        if 'sym' in c:
            j2,perm=c['sym']
            assert 0<=j2<q and j2!=j and 'sym' not in ch[j2]
            phi=make_phi(q,perm)
            for m2,(g2,rho2) in st.items():
                S={b for b in range(L) if b%g2==rho2%g2}
                assert {phi[b] for b in S}==S, ('restriction not preserved',st,m2)
            A={b for b in range(L) if b%(g*q)==(rho+g*j)%(g*q)}
            B={b for b in range(L) if b%(g*q)==(rho+g*j2)%(g*q)}
            assert {phi[b] for b in A}==B, ('child not mapped',st,m,q,j,j2)
            cnt['sym']+=1
        else:
            st2=dict(st); st2[m]=(g*q,rho+g*j); check(c,st2,depth+1)
sys.setrecursionlimit(10000)
try:
    check(cert['tree'],{m:(1,0) for m in BASE},0)
except AssertionError as e:
    print('CERTIFICATE REJECTED:',str(e)[:300]); sys.exit(1)
print('level',n,'L',L,'base moduli',BASE,'top moduli',TOP)
print('CERTIFICATE VERIFIED:',cnt,'distinct symmetries used',len(phi_ok))
