# Own test of Lemma ZC: #{rho: |gamma-T|<=1} <= 1.25 log Nf + 2.5 log(|k|+|T|+2) + 3.35
# on Groessencharacters of Q(i) with large infinity type k, via PARI gchar + lfunzeros.
import cypari2, math, mpmath
pari = cypari2.Pari()
pari.allocatemem(600*10**6)
pari.set_real_precision(38)
pari('bnf=bnfinit(y^2+1,1); gc=gcharinit(bnf,1); gc2=gcharinit(bnf,idealhnf(bnf,2+y))')
cases = [('gc',100,1,400), ('gc',250,1,1000), ('gc2',401,5,401), ('gc2',-999,5,-999)]
T0 = 25.0
for g, e, Nf, k in cases:
    pari('L=lfuncreate([%s,[%d,0]~])' % (g, e))
    cond = pari('L[5]')   # conductor of the L-data
    gam = pari('L[3]')
    z = pari('lfunzeros(L,[%f,%f])' % (-T0-1.5, T0+1.5))
    zs = sorted(float(t) for t in z)
    # smooth main term: theta(t) = t log A + Im logGamma(1/2+|k|/2+it), A = sqrt(Nf)/pi
    A = math.sqrt(Nf)/math.pi
    def th(t):
        return t*math.log(A) + float(mpmath.im(mpmath.loggamma(0.5+abs(k)/2+1j*t)))
    expected = (th(T0+1.5)-th(-T0-1.5))/math.pi
    worst = 0.0; worstT=None
    T = -T0
    while T <= T0+1e-9:
        cnt = sum(1 for t in zs if abs(t-T) <= 1)
        bound = 1.25*math.log(Nf) + 2.5*math.log(abs(k)+abs(T)+2) + 3.35
        r = cnt/bound
        if r > worst: worst, worstT = r, (T, cnt, bound)
        assert cnt <= bound
        T += 0.25
    # also the root number |W| = 1 check by PARI
    chk = float(pari('lfuncheckfeq(L)'))
    print(f"{g} chi=[{e}]: k={k} Nf={Nf} PARI cond={cond} gammaV={gam} lfuncheckfeq={chk:.0f}")
    print(f"   zeros on line in [-{T0+1.5},{T0+1.5}]: {len(zs)}, smooth main term {expected:.2f}")
    print(f"   max window count/bound {worst:.3f} at T={worstT[0]} (count {worstT[1]}, bound {worstT[2]:.2f})")
