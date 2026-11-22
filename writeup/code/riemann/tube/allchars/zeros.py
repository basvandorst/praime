"""Zeros of L(s,chi) on the critical line with |t| <= TMAX (sign changes of the rotated completed L-function)."""
import json, sys, mpmath as mp
from chars import chars
mp.mp.dps = 15
TMAX = float(sys.argv[1]) if len(sys.argv) > 1 else 20.0
STEP = 0.1
out = {}
only = sys.argv[2:] if len(sys.argv) > 2 else None
for name, q, v, a, real in chars():
    if only and name not in only:
        continue
    chi = [mp.mpc(z) for z in v]
    def Lam(t):
        s = mp.mpc(0.5, t)
        return mp.power(mp.mpf(q) / mp.pi, s / 2) * mp.gamma((s + a) / 2) * mp.dirichlet(s, chi)
    h0 = Lam(0.37)
    ph = mp.exp(-1j * mp.arg(h0))
    Z = lambda t: mp.re(Lam(t) * ph)
    ts = [-TMAX + STEP * i for i in range(int(2 * TMAX / STEP) + 1)]
    if real:
        ts = [t for t in ts if t >= -1e-9]
    vals = [Z(t) for t in ts]
    # sanity: imaginary part after rotation must be ~0
    imax = max(abs(mp.im(Lam(t) * ph)) / (abs(Lam(t)) + 1e-30) for t in ts[::25])
    zs = []
    for i in range(len(ts) - 1):
        if vals[i] == 0 or vals[i] * vals[i + 1] < 0:
            zs.append(float(mp.findroot(Z, (ts[i], ts[i + 1]), solver='anderson')))
    if real:
        zs = sorted(set([z for z in zs] + [-z for z in zs]))
    out[name] = dict(q=q, a=a, real=real, zeros=zs, rot_imag_rel=float(imax))
    print(name, 'rot check %.1e' % float(imax), 'lowest', sorted(zs, key=abs)[:4], 'n=%d' % len(zs), flush=True)
json.dump(out, open('zeros_%g.json' % TMAX, 'w'), indent=1)
