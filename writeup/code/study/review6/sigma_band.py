# empirical check of the truncation-noise model: zeros in a band [T1,T2] play the role of the omitted zeros
import numpy as np, math
g = np.array([float(l) for l in open('/home/bas/praime/writeup/data/riemann/zeros_hp2419.txt') if l.strip()])
rng = np.random.default_rng(1)
x = rng.uniform(1e9, 1.9e9, 4000)
def S(x, gs):
    lx = np.log(x)[:, None]; rho = 0.5 + 1j*gs[None, :]; w = rho*lx
    v = np.exp(w)/w*(1 + 1/w + 2/w**2)
    return 2*v.real.sum(1)
for T1 in (800.0, 1500.0, 2200.0):
    gs = g[g > T1]; T2 = gs[-1]
    s1 = np.concatenate([S(x[i:i+500], gs) for i in range(0, len(x), 500)])
    s2 = np.concatenate([S(x[i:i+500]+1e8, gs) for i in range(0, len(x), 500)])
    amp2 = (2*np.sqrt(x)/np.log(x))**2
    emp1 = np.mean(s1**2/amp2); empd = np.mean((s2-s1)**2/amp2)
    pred_sum = np.sum(0.5/(0.25+gs**2))                      # random phases, per-zero
    f = lambda T: (math.log(T/(2*math.pi))+1)/(2*math.pi*T)  # int_T^inf ln(t/2pi)/(2pi t^2)
    pred_int = 0.5*(f(T1)-f(T2))
    rep = lambda T: math.log(T/(2*math.pi))/(2*math.pi*T)    # repo's per-endpoint variance factor (no 1/2, no +1)
    print('band %.0f-%.0f (%d zeros): emp var/amp2 one end %.3e, increment %.3e; rand-phase per-zero %.3e, integral %.3e; repo-style %.3e'
          % (T1, T2, len(gs), emp1, empd, pred_sum, pred_int, rep(T1)-rep(T2)))
T=g[-1]; print('tail beyond T: rand-phase %.3e, repo %.3e, ratio of sigmas %.3f' % (0.5*(math.log(T/2/math.pi)+1)/(2*math.pi*T), math.log(T/2/math.pi)/(2*math.pi*T), math.sqrt(2*math.log(T/2/math.pi)/(math.log(T/2/math.pi)+1))))
