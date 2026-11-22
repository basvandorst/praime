import mpmath as mp
mp.mp.dps = 40
def L(s, chi, d=0): return mp.dirichlet(s, chi, d)
for name, q, chi in (('chi4',4,[0,1,0,-1]),('chi3',3,[0,1,-1]),('chi5c',5,[0,1,1j,-1j,-1])):
    d1 = L(1,chi,1)/L(1,chi); d0 = L(0,chi,1)/L(0,chi)
    Ca = d1 - d0 + mp.log(2)
    Cb = 2*mp.re(d1) + mp.log(mp.mpf(q)/mp.pi) - mp.euler
    print(name, 'L0', mp.nstr(L(0,chi),15), 'd1', mp.nstr(d1,15), 'd0', mp.nstr(d0,15), 'Ca', mp.nstr(Ca,15), 'Cb', mp.nstr(Cb,15))
# independent closed form for chi4: L'(1)/L(1) = gamma + 2 ln 2 + 3 ln pi - 4 ln Gamma(1/4)
print('chi4 closed d1', mp.nstr(mp.euler + 2*mp.log(2) + 3*mp.log(mp.pi) - 4*mp.log(mp.gamma(0.25)),15))
# d0 independent via derivative of Hurwitz: L(s)=4^-s (zeta(s,1/4)-zeta(s,3/4))
f = lambda s: mp.power(4,-s)*(mp.zeta(s,0.25)-mp.zeta(s,0.75))
print('chi4 d0 hurwitz', mp.nstr(mp.diff(f,0)/f(0),15), 'd1 hurwitz', mp.nstr(mp.diff(f,1)/f(1),15))
g = lambda s: mp.power(3,-s)*(mp.zeta(s,mp.mpf(1)/3)-mp.zeta(s,mp.mpf(2)/3))
print('chi3 d0 hurwitz', mp.nstr(mp.diff(g,0)/g(0),15), 'd1', mp.nstr(mp.diff(g,1)/g(1),15))
