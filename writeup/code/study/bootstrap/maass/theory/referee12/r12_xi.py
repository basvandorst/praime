# Referee code, stage 1 (own code, no repo imports).
# Xi_Delta(t) = Lambda(1/2+it, Delta) from the q-expansion Mellin integral, split at y = 1:
#   Xi(t) = 4 sum_n tau(n) Re[(2 pi n)^{-(6+it)} Gamma(6+it, 2 pi n)]
# tau(n) from my own product q prod (1-q^n)^24.  Values at Gauss-Legendre nodes on [0,T]
# (unit panels, NG nodes each) are stored for stage 2.  Also: normalisation checks.
import mpmath as mp, json, sys, time
mp.mp.dps = 60
T = 56
NTERMS = 26

def tau_list(N):
    poly = [0]*N; poly[0] = 1
    for n in range(1, N):
        for _ in range(24):
            for k in range(N-1, n-1, -1):
                poly[k] -= poly[k-n]
    return [0] + poly  # tau[k] = coeff of q^k

TAU = tau_list(402)

def xi(t):
    w = mp.mpc(6, t); s = mp.mpf(0)
    for n in range(1, NTERMS+1):
        a = 2*mp.pi*n
        s += TAU[n]*mp.re(a**(-w)*mp.gammainc(w, a))
    return 4*s

def lam_sym(s):
    # Lambda(s) via the split Mellin integral (uses the functional equation)
    s = mp.mpc(s); tot = 0
    for n in range(1, NTERMS+1):
        a = 2*mp.pi*n
        tot += TAU[n]*(a**(-(s+mp.mpf(11)/2))*mp.gammainc(s+mp.mpf(11)/2, a)
                       + a**(-(mp.mpf(13)/2-s))*mp.gammainc(mp.mpf(13)/2-s, a))
    return 2*tot

def gl_nodes(n):
    xs, ws = [], []
    for i in range(1, n+1):
        x = mp.cos(mp.pi*(i-mp.mpf(1)/4)/(n+mp.mpf(1)/2))
        for _ in range(100):
            p0, p1 = mp.mpf(1), x
            for k in range(2, n+1):
                p0, p1 = p1, ((2*k-1)*x*p1-(k-1)*p0)/k
            dp = n*(x*p1-p0)/(x*x-1)
            dx = p1/dp; x -= dx
            if abs(dx) < mp.mpf(10)**(-mp.mp.dps+5): break
        p0, p1 = mp.mpf(1), x
        for k in range(2, n+1):
            p0, p1 = p1, ((2*k-1)*x*p1-(k-1)*p0)/k
        dp = n*(x*p1-p0)/(x*x-1)
        xs.append(x); ws.append(2/((1-x*x)*dp*dp))
    return xs, ws

if __name__ == '__main__':
    import cypari2
    pari = cypari2.Pari(); pari.set_real_precision(60)
    out = {}
    # checks of tau and normalisation
    ptau = [int(pari('ramanujantau(%d)' % n)) for n in range(1, 221)]
    out["tau_ok"] = (ptau == TAU[1:221])
    # Lambda(10): Dirichlet series vs split Mellin
    s = mp.mpf(10)
    ds = 2*(2*mp.pi)**(-(s+5.5))*mp.gamma(s+5.5)*mp.fsum(TAU[n]*mp.mpf(n)**(-(s+5.5)) for n in range(1, 221))
    out['lam10_dirichlet'] = mp.nstr(ds, 25); out['lam10_mellin'] = mp.nstr(lam_sym(s), 25)
    # Lambda(s) = Lambda(1-s) is built into lam_sym; check the FE through the Dirichlet series at s=10
    # against lam_sym(1-10) (same formula, but the Dirichlet side is independent)
    out['lam_m9_mellin'] = mp.nstr(lam_sym(1-s), 25)
    for t in [0, 5.5, 20, 40]:
        mine = xi(mp.mpf(t))
        pv = pari('lfunlambda(lfunetaquo(Mat([1,24])), 6+%s*I)' % t)
        pv = mp.mpf(str(pari.real(pv)).replace(' E', 'E'))
        out['xi_%s' % t] = [mp.nstr(mine, 30), mp.nstr(pv, 30), mp.nstr(abs(mine-pv), 3)]
    print(json.dumps(out, indent=1)); sys.stdout.flush()
    for NG in [20, 24]:
        t0 = time.time()
        xs, ws = gl_nodes(NG)
        nodes, wts, vals = [], [], []
        for k in range(T):
            for x, w in zip(xs, ws):
                t = k + (x+1)/2
                nodes.append(t); wts.append(w/2); vals.append(xi(t))
        json.dump({'NG': NG, 'T': T, 'dps': mp.mp.dps,
                   't': [mp.nstr(v, 70) for v in nodes], 'w': [mp.nstr(v, 70) for v in wts],
                   'xi': [mp.nstr(v, 70) for v in vals]},
                  open('xi_nodes_%d.json' % NG, 'w'))
        print('NG', NG, 'done', time.time()-t0, 'xi(T)~', mp.nstr(vals[-1], 5)); sys.stdout.flush()
