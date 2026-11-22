import sympy as sp
d = sp.symbols('delta', positive=True)
D = -sp.EulerGamma - sp.log(2*sp.pi) + 2*sp.log(1+d) + ((2+d)/d)*sp.log(1+d/2)
HL = 1 - sp.EulerGamma - sp.log(2*sp.pi) + (1+d)*sp.log(1+d)/d - 1
print('D-HL series:', sp.series(D-HL, d, 0, 6))
# pieces: pole term (2/d) h(i/2) /... : Weil: D = (1/d)[2 h(i/2) - g0 ln pi + dig] ; 2h(i/2)/d = 2 ln(1+d)
print('2 ln(1+d) series:', sp.series(2*sp.log(1+d), d, 0, 4))
# "+1 in N" by parts: (2/d) h(0) with h(0) = (sqrt(1+d)-1)^2*4
plus1 = (2/d)*4*(sp.sqrt(1+d)-1)**2
print('+1 term (2/d)h(0):', sp.series(plus1, d, 0, 4))
print('-1/8 term:', sp.series(-plus1/8, d, 0, 4))
