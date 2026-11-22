# own explicit-formula prediction for blocks of 1e8 in [1e9,2e9)
import mpmath as mp, json, math, time
mp.mp.dps = 25
Z = [mp.mpf(l.strip()) for l in open('/home/bas/praime/writeup/data/riemann/zeros_hp2419.txt') if l.strip()]
T = float(Z[-1])
def mob(n):
    r, m, p = 1, n, 2
    while p*p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0: return 0
            r = -r
        p += 1
    return -r if m > 1 else r
def R(y, N=200):
    ly = mp.log(y)
    return sum(mob(n)*mp.ei(ly/n)/n for n in range(1, N) if mob(n))   # li(y^(1/n)) = Ei(ln y / n)
def zterm(y, frac=1):
    # sum over zeros (both conjugates) of li(y^(rho*frac)) = 2 Re Ei(rho*frac*ln y)
    ly = mp.log(y)
    return sum(2*mp.re(mp.ei(mp.mpc(0.5, g)*frac*ly)) for g in Z)
# branch check: Ei(rho ln y) vs asymptotic y^rho/(rho ln y) * sum k!/(rho ln y)^k for one zero
y = mp.mpf(10)**9; rho = mp.mpc(0.5, Z[0]); w = rho*mp.log(y)
asym = mp.exp(w)/w*sum(mp.factorial(k)/w**k for k in range(12))
print('branch check zero 1:', mp.ei(w), asym)
edges = [10**9 + k*10**8 for k in range(11)]
t0 = time.time()
out = {}
for e in edges:
    out[e] = dict(li=float(mp.li(e)), R=float(R(e)), Z1=float(zterm(e)), Z2=float(zterm(e, mp.mpf(1)/2)), Z3=float(zterm(e, mp.mpf(1)/3)))
    print(e, out[e], round(time.time()-t0), flush=True)
json.dump({str(k): v for k, v in out.items()}, open('explicit.json', 'w'), indent=1)
