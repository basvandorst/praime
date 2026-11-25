"""Second route (task 3), light version run with unbuffered output: (1) CUE(N) <r~> exactly by Janossy densities;
(2) the wave closure by brute force on Haar CUE(300) spectra (paired, with and without p's wave; random phase per
matrix), against the deterministic closure (percell.json); the dilation comparator on the same spectra."""
import sys, os, math, json, time
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H + '/lib')
import numpy as np, fred, closure
t0 = time.time()
r_inf, _ = fred.ratio_mean(fred.sinc_kernel)
print('(1) CUE(N) <r~> exact; sine kernel %.6f' % r_inf, flush=True)
cue = {}
for N in (3, 4, 6, 8, 12):
    K = lambda x, N=N: np.where(np.abs(x) < 1e-12, 1.0, np.sin(np.pi * x) / (N * np.sin(np.pi * x / N) + (np.abs(x) < 1e-12)))
    r, nrm = fred.ratio_mean(K, nT=32, nf=16, Tmax=min(5.0, N - 1e-9), m=24); cue[N] = (r, nrm)
    print('   N=%2d  <r~> %.6f  norm %.5f   N^2 (r - r_inf) %+.5f   N^4 (r - r_inf) %+.4f' % (N, r, nrm, N * N * (r - r_inf), N ** 4 * (r - r_inf)), flush=True)
rng = np.random.default_rng(20261125)
def cue_angles(N):
    Z = (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))) / math.sqrt(2)
    Q, R = np.linalg.qr(Z); d = np.diagonal(R); Q = Q * (d / np.abs(d))
    return np.sort(np.angle(np.linalg.eigvals(Q)))
rr = lambda s: np.minimum(s[1:], s[:-1]) / np.maximum(s[1:], s[:-1])
N = 300; M = 120; L = 9.4; res = {}
print('(2) brute force: %d Haar CUE(%d), interior 240 points, L = %.1f (missing minus present)' % (M, N, L), flush=True)
specs = [cue_angles(N) * N / (2 * math.pi) for _ in range(M)]
for p in (2, 3, 5, 13):
    w = closure.Wave(p, L); acc = []
    for y in specs:
        y = y[30:270]; phi = rng.uniform(0, 2 * math.pi); x = y.copy()
        for _ in range(40): x = y - w.dN(x, phi)
        s0 = np.diff(y); s1 = np.diff(x); sl = s0 / w.rho(x[:-1], phi)
        acc.append([np.mean((s1 - s1.mean()) ** 2) - np.var(s0), rr(s1).mean() - rr(s0).mean(), np.var(sl) - np.var(s0), rr(sl).mean() - rr(s0).mean()])
    a = np.array(acc); m = -a.mean(0); se = a.std(0) / math.sqrt(M)
    res[p] = dict(dvar=m[0], dvar_se=se[0], dr=m[1], dr_se=se[1], dvar_slow=m[2], dr_slow=m[3])
    print('   p=%2d  dVar %+.5f +- %.5f   d<r~> %+.5f +- %.5f  | dilation: dVar %+.5f d<r~> %+.5f' % (p, m[0], se[0], m[1], se[1], m[2], m[3]), flush=True)
pc = {c['p']: c for c in json.load(open(H + '/percell.json'))['cells'] if c['L'] == 9.4}
print('   deterministic closure: ' + '; '.join('p=%d dVar %+.5f d<r~> %+.5f (slow %+.5f)' % (p, pc[p]['dvar_closure'], pc[p]['dr_closure'], pc[p]['dvar_slow']) for p in (2, 3, 5, 13)))
json.dump({'cue_r': {str(k): v for k, v in cue.items()}, 'r_inf': r_inf, 'sim': {str(k): v for k, v in res.items()}}, open(H + '/sim_check.json', 'w'), indent=1)
print('time %.0f s' % (time.time() - t0))
