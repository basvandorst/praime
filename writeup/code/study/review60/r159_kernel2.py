import mpmath as mp
mp.mp.dps = 20
exec(open(__file__.replace("r159_kernel2.py", "r159_kernel.py")).read().split("mp.mp.dps = 20")[0])
for w, k in ((mp.mpf(-0.6), 2), (mp.mpf(-1.3), 4)):
    X = 400 * mp.pi; f = lambda x: x ** (w - mp.mpf(1) / 2) * mp.besselj(k - 1, x) * mp.cos(x)
    body = mp.quad(f, mp.linspace(0, X, 801))
    phi = (k - 1) * mp.pi / 2 + mp.pi / 4
    tail = mp.cos(phi) / mp.sqrt(2 * mp.pi) * X ** w / (-w)
    a = I_closed(w, k); b = body + tail
    print("I(%s) k=%d closed %s  quad+tail %s  rel %.1e" % (mp.nstr(w, 3), k, mp.nstr(a, 10), mp.nstr(b, 10), abs(a / b - 1)))
