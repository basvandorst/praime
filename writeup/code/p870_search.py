import itertools, sys
from sympy import isprime
S=[[0,1,2,3],[2,3,0,1],[3,2,1,0],[1,0,3,2]]   # orthogonal diagonal Latin squares of order 4
T=[[0,1,2,3],[3,2,1,0],[1,0,3,2],[2,3,0,1]]
assert len({(S[i][j],T[i][j]) for i in range(4) for j in range(4)})==16
for M in (S,T):
    assert all(sorted(r)==[0,1,2,3] for r in M) and all(sorted(M[i][j] for i in range(4))==[0,1,2,3] for j in range(4))
    assert sorted(M[i][i] for i in range(4))==[0,1,2,3] and sorted(M[i][3-i] for i in range(4))==[0,1,2,3]
def search(rule):
    # inner entries 3-digit primes a+b, a even in [0,500), b odd in [101,500)
    good=lambda s: 100<=s<1000 and isprime(s) and isprime(rule(s))
    As=range(0,500,2); Bs=range(101,500,2)
    nb={b:{a for a in As if good(a+b)} for b in Bs}
    best=None
    for b4 in itertools.combinations(Bs,4):
        common=nb[b4[0]]&nb[b4[1]]
        if len(common)<4: continue
        common&=nb[b4[2]]
        if len(common)<4: continue
        common&=nb[b4[3]]
        for a4 in itertools.combinations(sorted(common),4):
            sums=[a+b for a in a4 for b in b4]
            if len(set(sums))<16: continue
            sq=[[a4[S[i][j]]+b4[T[i][j]] for j in range(4)] for i in range(4)]
            return sq
rules={'R  (n&3)':lambda s:10*s+3,'R  (n&1)':lambda s:10*s+1,'L  (1&n)':lambda s:1000+s,'L  (3&n)':lambda s:3000+s,
       'LR (1&n&3)':lambda s:(1000+s)*10+3,'LR (2&n&1)':lambda s:(2000+s)*10+1}
for name,r in rules.items():
    sq=search(r); print(name, sq, [[r(x) for x in row] for row in sq] if sq else None, flush=True)
