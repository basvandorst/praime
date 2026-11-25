"""g2_check.py: can the degree-4 line of a genus-2 curve be seen on primes p <= X (and so check PARI's conductor)?
For y^2 + Q(x) y = P(x): lambda(p) = a_p / sqrt p, a_p = p + 1 - #C(F_p) (PARI hyperellcharpoly), bad p skipped.
Statistic (degree d = 4): S(Q') = sum_p lambda(p) ln p p^{-3/8} W(p^{1/4}) e(-4 (p/Q')^{1/4}), W = Hann in v = p^{1/4} over
[v0, v1]; unit main term (c = 1, eps = 1): (5/2) K1 int W dv, K1 = e^{-i 5pi/4} 0.8 Q^{-1/8}.  Reported: z = |S|/sigma at
the PARI conductor, the reading S/(unit main term), and the largest z over Q' in [Q/3, 3Q] (grid in Q' fine enough:
the line width in Q' is ~ 4 Q / (alpha (v1 - v0)) ), plus the predicted S/sigma for c = 1.
usage: g2_check.py X"""
import sys, os, json, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import spectro_lib as S

X = int(float(sys.argv[1])); P = S.pari()
CURVES = {'249a': ('x^2+x', 'x^3+1'), '277a': ('-x^2-x', 'x^3+x^2+x+1'), '295a': ('-x^4-x^3', 'x^3+1'), '349a': ('-x^2', 'x^3+x^2+x+1'), '353a': ('x^2', 'x^3+x+1')}
out = {}
for name, (Pp, Qp) in CURVES.items():
    t0 = time.time()
    red = P('genus2red([%s, %s])' % (Pp, Qp)); N = int(red[0])
    v = P('my(N=%d); [Vecsmall(primes([3,%d])), Vecsmall(apply(p->if(N%%p, my(f=hyperellcharpoly(Mod(1,p)*[%s, %s])); -polcoef(f,3), 0), primes([3,%d])))]' % (N, X, Pp, Qp, X))
    p = np.array(list(v[0]), np.int64); a = np.array(list(v[1]), np.float64); pf = p.astype(np.float64); lam = a / np.sqrt(pf)
    vv = pf ** 0.25; v0, v1 = 1e3 ** 0.25, X ** 0.25
    w = lam * np.log(pf) * pf ** -0.375 * S.hann(vv, v0, v1); sig = math.sqrt(float(np.dot(w, w)) / 2)
    uu = np.linspace(v0, v1, 20001); IW = np.trapezoid(S.hann(uu, v0, v1), uu)
    def unit(Qc): return 2.5 * complex(math.cos(-5 * math.pi / 4), math.sin(-5 * math.pi / 4)) * 0.8 * Qc ** -0.125 * IW
    def Ssum(Qc): return complex(np.dot(w, np.exp(-1j * S.TWO_PI * 4 * (pf / Qc) ** 0.25)))
    St = Ssum(N); zt = abs(St) / sig; rd = St / unit(N)
    # scan Q' (log grid; line width in log Q' ~ 4 / (alpha (v1 - v0)), alpha = 4 Q^{-1/4})
    al = 4 * N ** -0.25; dlog = 4.0 / (al * (v1 - v0)) / 4
    grid = np.exp(np.arange(math.log(N / 3), math.log(3 * N), dlog))
    zs = np.array([abs(Ssum(q)) / sig for q in grid])
    out[name] = dict(N=N, z_at_N=zt, reading_c1=[rd.real, rd.imag], pred_z_c1=abs(unit(N)) / sig, max_z_scan=float(zs.max()),
                     at=float(grid[int(np.argmax(zs))]), n_grid=len(grid), sigma=sig)
    print('%s N=%d: z at N = %.2f (predicted for c = 1: %.2f), reading S/unit = %+.2f %+.2fi; max z over %d values of Q in [N/3, 3N] = %.2f at %.0f  %.0fs'
          % (name, N, zt, abs(unit(N)) / sig, rd.real, rd.imag, len(grid), zs.max(), grid[int(np.argmax(zs))], time.time() - t0), flush=True)
json.dump(out, open(os.path.join(S.HERE, 'out_g2_%d.json' % X), 'w'), indent=1)
