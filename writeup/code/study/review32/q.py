"""Independent Q: positive ordinates only, explicit sign bookkeeping."""
import numpy as np, math, sys
Z = np.loadtxt('/home/bas/praime/writeup/data/riemann/zeros_10000_flint.txt')
DS = (0.1, 0.03, 0.01)
def w(g, d):
    rho = 0.5 + 1j*g
    return np.abs(((1+d)**rho - 1)/rho)**2
def Qparts(g, ds=DS, ch=1000):
    W = np.stack([w(g, d) for d in ds], 1)       # n x 3
    same = np.zeros(len(ds)); opp = np.zeros(len(ds))
    for i0 in range(0, len(g), ch):
        a = g[i0:i0+ch, None]
        Km = 1/(4 + (a - g[None, :])**2)
        idx = np.arange(i0, min(i0+ch, len(g))); Km[idx-i0, idx] = 0
        same += np.einsum('ij,jd,id->d', Km, W, W[i0:i0+ch])
        Kp = 1/(4 + (a + g[None, :])**2)
        opp += np.einsum('ij,jd,id->d', Kp, W, W[i0:i0+ch])
    selfc = np.sum(W**2/(4+4*g[:, None]**2), 0)   # (+a,-a) pairs, per sign
    return W, same, opp, selfc
if __name__ == '__main__':
    W, same, opp, selfc = Qparts(Z)
    Q = 2*same + 2*opp      # both signs, j != k
    for i, d in enumerate(DS):
        ms = 4/d**2*(2*Q[i] - 2*selfc[i])
        print(d, 'Q', Q[i], 'same', 2*same[i], 'opp', 2*opp[i], 'rms(2/d)sqrt(2Q)', 2/d*math.sqrt(2*Q[i]), 'exact-LI rms', math.sqrt(ms), 'rms(2/d)sqrtQ', 2/d*math.sqrt(Q[i]))
