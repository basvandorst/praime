import mpmath as mp
mp.mp.dps=30
# M(s)=int_0^inf t^(s-2)(1-cos t)dt = -Gamma(s-1)cos(pi(s-1)/2), -1<Re s<1; M(0)=pi/2, M'(0)=int ln t (1-cos t)/t^2
M=lambda s: -mp.gamma(s-1)*mp.cos(mp.pi*(s-1)/2)
e=mp.mpf("1e-12"); print(M(e), mp.pi/2); print(mp.diff(M,mp.mpf("1e-8")), mp.pi/2*(1-mp.euler))
# independent numeric check of M(0.3) by direct quadrature
s=mp.mpf('0.3'); f=lambda t: t**(s-2)*(1-mp.cos(t))
print(mp.quad(f,[0,1])+mp.quadosc(f,[1,mp.inf],period=2*mp.pi), M(s))
