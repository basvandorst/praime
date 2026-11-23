import cypari2, json, cmath
pari = cypari2.Pari()
z=json.load(open('zeros131.json'))
print(z['5']['chars']['2'][:3])
q=7; G=pari.znstar(q,1); g=3
for k in [1,2]:
    L=pari.lfuncreate([G,[k]])
    v=complex(pari.lfun(L,2))
    # my chi: chi(g^j)=exp(2 pi i k j/6)
    dl={pow(g,j,q):j for j in range(6)}
    s=sum(cmath.exp(2j*cmath.pi*k*dl[m%q]/6)/m**2 for m in range(1,200000) if m%q)
    print(k,v,s)
