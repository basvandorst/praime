import numpy as np, sys, time
from flow import load, integrate

M = int(sys.argv[1]); c = float(sys.argv[2]); dtmax = float(sys.argv[3]); out = sys.argv[4]
tmax = float(sys.argv[5]) if len(sys.argv) > 5 else 0.6
x0, Y0 = load(M)
tsnap = np.round(np.arange(0, tmax + 1e-9, 0.05), 10)
t0 = time.time()
snaps, nsteps = integrate(x0, Y0, tsnap, c, dtmax, log=lambda s: print(s, '%.0fs' % (time.time() - t0), flush=True))
np.savez(out, t=tsnap, x=snaps, Y0=Y0, M=M, c=c, dtmax=dtmax, nsteps=nsteps)
print('done', nsteps, 'steps', '%.1f s' % (time.time() - t0))
