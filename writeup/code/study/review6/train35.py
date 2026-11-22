import mpmath as mp, math, json, numpy as np
exec(open('explicit.py').read().split("# branch check")[0])   # reuse my own R, zterm, mob (this directory)
from sieve_blocks_lib import count
edges = [10**7 + k*9*10**6 for k in range(11)]
E = {e: (float(mp.li(e)), float(R(e)), float(zterm(e)), float(zterm(e, mp.mpf(1)/2))) for e in edges}
def sig_rp(x): return (2*math.sqrt(x)/math.log(x))*math.sqrt(0.5*(math.log(T/(2*math.pi))+1)/(2*math.pi*T))
chi = 0
for k in range(10):
    a, b = edges[k], edges[k+1]; c = count(a, b)
    pred = (E[b][1]-E[a][1]) - (E[b][2]-E[a][2]); dli = E[b][0]-E[a][0]
    s = math.hypot(sig_rp(a), sig_rp(b)); z = (c-pred)/s; chi += z*z
    print(a, c, round(c-dli,1), round(pred-dli,1), round(c-pred,1), 'sigma_rp %.1f z %.2f' % (s, z), 'n=2 term %.1f' % (0.5*(E[b][3]-E[a][3])))
print('chi2 (rand-phase) %.1f' % chi)
