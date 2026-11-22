import mpmath as mp, math, json
mp.mp.dps = 20
q = 4
def L(s): return mp.power(4,-s)*(mp.zeta(s,0.25)-mp.zeta(s,0.75))
def Z(t):
    s = mp.mpf(0.5)+1j*t
    th = mp.im(mp.loggamma((s+1)/2)) + t/2*mp.log(q/mp.pi)
    v = mp.exp(1j*th)*L(s)
    return v
T = 200.0; h = 0.04
ts = [i*h for i in range(1, int(T/h)+1)]
vals = []
maxim = 0
for t in ts:
    v = Z(t); maxim = max(maxim, abs(mp.im(v))/(abs(v)+1e-30)); vals.append(float(mp.re(v)))
print('max rel imag part of Z (should be ~0):', maxim)
zs = []
for i in range(len(ts)-1):
    if vals[i]*vals[i+1] < 0:
        zs.append(float(mp.findroot(lambda t: mp.re(Z(t)), (ts[i], ts[i+1]), solver='anderson')))
N = len(zs)
Nexp = T/(2*math.pi)*math.log(q*T/(2*math.pi*math.e))
print('zeros found', N, 'first', zs[:4], 'smooth N(T)', Nexp)
S = 2*sum(1/(0.25+g*g) for g in zs)
tail = 2*(1/(2*math.pi*T))*(math.log(q*T/(2*math.pi))+1)
print('sum to T', S, 'tail est', tail, 'total', S+tail, 'claimed 0.155568')
json.dump(zs, open('/home/bas/scratch/study-review33/zeros4.json','w'))
