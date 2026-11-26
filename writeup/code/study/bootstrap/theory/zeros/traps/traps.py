"""D1: trap data for every stored degree and node (numerics, arb ~ 100 digits; NOT a proof).

For the optimum gt = g_d (hp_d<d>.json, re-polished by polish_hp.py; gt(0) = 1, kt = W(gt) = lambda) and the
interpolation element h_k (h(u_j) = h'(u_j) = 0 for j != k, h(u_k) = 1, h'(u_k) = 0, top coefficient 0):
  eta_0 = lambda / W(h_k)                      (threshold of (ii) when eta h(0) <= 1)
  eta_1 = min_{t real, B_h(t) > 0} B_gt(t)/B_h(t)   (sup of the eta for which (i) holds; inf if B_h <= 0)
  (iii) holds for every eta >= 0 (b_d(h_k) = 0 < b_d(gt)).
For eta = 2 eta_0 and 10 eta_0 (when < eta_1): S = {u > 0 : w = eta h - gt > 0}.
  main component (a, b) around u_k: safeguarded Newton from u_k +- sqrt(2 eta / gt''(u_k));
  sup_S w = w(u_k) if w'' < 0 at 9 equally spaced points of [a, b] (u_k is the critical point), else the max of w
  at 41 equally spaced points and u_k;
  other components: sign of w on the grid u = i/1000 (i <= 4000), 4 + i/100 (i <= 400) (grid.pkl), the limits
  rho_j = eta h''(u_j)/gt''(u_j) at the other nodes (w ~ (eta h'' - gt'')(u-u_j)^2/2 there), rho(0) = eta h(0),
  and the top-coefficient sign at infinity (w -> -gt);
  margin = max of rho = eta h/gt over grid points with u < a - max(b-a, 0.03) or u > b + max(b-a, 0.03);
  other_logs_in_S = integers m != n_k with log m in [a, b] (the trap does not isolate n_k if non-empty).
  weight bound = (eta W(h) - lambda (1 + (eta h(0) - 1)^+)) / sup_S w.
Usage: python traps.py d [d ...]   -> traps_d<d>.json and lines in traps_d<d>.out"""
import sys, json, pickle, math, time
from fam import *
import mpmath as mp

mp.mp.dps = 60
PP = {}
for q in range(2, 200):
    f = [p for p in range(2, q + 1) if q % p == 0 and all(p % r for r in range(2, p))]
    if len(f) == 1:
        PP[q] = f[0]
GRID = pickle.load(open("grid.pkl", "rb"))
ctx.prec = 450
GU = [arb(fmpq(*map(int, s.split("/")))) for s in GRID["u"]]
GUf = [float(Fraction(s)) for s in GRID["u"]]
GROWS = [[arb(x) for x in row] for row in GRID["rows"]]
ctx.prec = 500


def fl(x):
    return float(x.mid().str(20, radius=False)) if isinstance(x, arb) else float(x)


def eta1_of(bg, bh):
    """min over s = t^2 >= 0 with q_h(s) > 0 of q_g(s)/q_h(s)."""
    qg = [mp.mpf(x.mid().str(130, radius=False)) for x in bg][::-1]
    qh = [mp.mpf(x.mid().str(130, radius=False)) for x in bh][::-1]
    dqg = mp.polyval  # noqa
    ts = [mp.mpf(0)] + [mp.mpf(10) ** (mp.mpf(i) / 400) for i in range(-800, 1601)]
    vals = []
    for t in ts:
        s = t * t
        h = mp.polyval(qh, s)
        vals.append(mp.polyval(qg, s) / h if h > 0 else mp.inf)
    best = min(vals)
    if best == mp.inf:
        return mp.inf, None
    # refine local minima (golden section in t)
    cands = []
    for i in range(len(ts)):
        if vals[i] == mp.inf:
            continue
        l = vals[i - 1] if i > 0 else mp.inf
        r = vals[i + 1] if i + 1 < len(ts) else mp.inf
        if vals[i] <= l and vals[i] <= r:
            cands.append(i)
    res = []
    for i in cands:
        lo = ts[max(i - 1, 0)]
        hi = ts[min(i + 1, len(ts) - 1)]
        f = lambda t: (mp.polyval(qg, t * t) / mp.polyval(qh, t * t)) if mp.polyval(qh, t * t) > 0 else mp.inf
        gr = (mp.sqrt(5) - 1) / 2
        a, b = lo, hi
        c1 = b - gr * (b - a); c2 = a + gr * (b - a)
        f1, f2 = f(c1), f(c2)
        for _ in range(200):
            if f1 < f2:
                b, c2, f2 = c2, c1, f1
                c1 = b - gr * (b - a); f1 = f(c1)
            else:
                a, c1, f1 = c1, c2, f2
                c2 = a + gr * (b - a); f2 = f(c2)
            if b - a < mp.mpf(10) ** -40 * (1 + b):
                break
        res.append((min(f1, f2, vals[i]), (a + b) / 2))
    v, t = min(res, key=lambda x: x[0])
    return v, t


def run(d):
    t0 = time.time()
    R = json.load(open("hp_d%d.json" % d))
    K = R["K"]
    n = d + 1
    F = Fam(d, 2)
    bg = [arb(s) for s in R["b"]]
    us = [arb(s) for s in R["u"]]
    lam = arb(R["lam"])
    mu = [arb(s) for s in R["mu"]]
    gpp = [arb(s) for s in R["gpp"]]
    Pg = [F.comb(bg, r) for r in range(3)]
    g_grid = [dot(bg, row[:n]) for row in GROWS]
    neg_g = [i for i, v in enumerate(g_grid) if v < 0]
    out = dict(d=d, deg=2 * d, K=K, gt_negative_on_grid=len(neg_g), nodes=[])
    lines = []
    for k in range(K):
        H = R["h"][k]
        bh = [arb(s) for s in H["b"]]
        W = arb(H["W"]); h0 = arb(H["h0"])
        uk = us[k]
        e_u = math.exp(fl(uk))
        nk = min(PP, key=lambda q: abs(math.log(q) - fl(uk)))
        logn = arb(nk).log()
        lamn = 2 * arb(PP[nk]).log() / arb(nk).sqrt()
        eta0 = lam / W
        e1, t1 = eta1_of(bg, bh)
        eta1 = arb(mp.nstr(e1, 40)) if e1 != mp.inf else None
        hpp = [arb(s) for s in H["hpp"]]
        rho_nodes = [fl(x) for x in [None] * 0]
        Ph = [F.comb(bh, r) for r in range(3)]
        h_grid = [dot(bh, row[:n]) for row in GROWS]
        rec = dict(k=k + 1, n=nk, e_u=e_u, u_minus_logn=astr(uk - logn, 20), mu=astr(mu[k], 30),
                   cost_weight=astr(lamn, 30), mu_minus_costweight=astr(mu[k] - lamn, 10), gpp=astr(gpp[k], 30),
                   W=astr(W, 30), kkt_resid=astr(W - lam * h0 - mu[k], 5), h0=astr(h0, 20), eta0=astr(eta0, 20),
                   eta1=(astr(eta1, 20) if eta1 is not None else "inf"),
                   eta1_over_eta0=(fl(eta1 / eta0) if eta1 is not None else "inf"),
                   eta1_at_t=(mp.nstr(t1, 8) if t1 is not None else None), etas=[])
        for fac in (2, 10):
            eta = eta0 * fac
            E = dict(factor=fac, eta=astr(eta, 20))
            if eta1 is not None and not (eta < eta1):
                E["status"] = "(i) fails: eta >= eta_1"
                rec["etas"].append(E)
                continue

            def wv(u, r=0):
                v = F.eval_polys([Pg[r], Ph[r]], u.exp())
                return eta * v[1] - v[0]

            ends = []
            r0 = (2 * eta / gpp[k]).sqrt()
            okb = True
            for s in (-1, 1):
                L = r0 * arb("1.5")
                while wv(uk + s * L) > 0 and fl(L) < 0.5:
                    L = L * 2
                if wv(uk + s * L) > 0:
                    okb = False
                    ends.append(None)
                    continue
                lo, hi = arb(0), L            # w(uk + s lo) > 0, w(uk + s hi) < 0
                x = r0 if fl(r0) < fl(L) else L / 2
                for it in range(80):
                    ux = mid(uk + s * x)
                    w0 = wv(ux)
                    w1 = wv(ux, 1) * s        # derivative in x
                    if w0 > 0:
                        lo = x
                    else:
                        hi = x
                    step = -w0 / w1 if not w1.contains(0) else None
                    xn = mid(x + step) if step is not None else None
                    if xn is None or not (xn > lo and xn < hi):
                        xn = mid((lo + hi) / 2)
                    if abs(fl(xn - x)) < 1e-40 * fl(r0):
                        x = xn
                        break
                    x = xn
                ends.append(mid(uk + s * x))
            if not okb:
                E["status"] = "main component wider than 0.5 (not bracketed)"
                rec["etas"].append(E)
                continue
            a, b = ends
            conc = [fl(wv(mid(a + (b - a) * arb(i) / 8), 2)) for i in range(9)]
            if all(x < 0 for x in conc):          # concave at the samples: max at the critical point u_k (w'(u_k) = 0)
                supw = wv(uk)
            else:
                samp = [mid(a + (b - a) * arb(i) / 40) for i in range(41)] + [uk]
                supw = max([wv(x) for x in samp], key=lambda x: fl(x))
            # other components on the grid
            extra = []
            run_ = None
            rmax, rmax_u = -1e300, None
            for i, (gv, hv) in enumerate(zip(g_grid, h_grid)):
                u = GUf[i]
                w = eta * hv - gv
                inside = (GU[i] > a) and (GU[i] < b)
                if w > 0 and not inside:
                    if run_ is None:
                        run_ = [u, u]
                    else:
                        run_[1] = u
                else:
                    if run_ is not None:
                        extra.append(run_); run_ = None
                if (u < fl(a) - max(fl(b - a), 0.03) or u > fl(b) + max(fl(b - a), 0.03)) and gv > 0:
                    rho = fl(eta * hv / gv)
                    if rho > rmax:
                        rmax, rmax_u = rho, u
            if run_ is not None:
                extra.append(run_)
            rho_j = [(j + 1, fl(eta * hpp[j] / gpp[j])) for j in range(K) if j != k]
            rho_j_max = max([x[1] for x in rho_j], default=-1e300)
            rho0 = fl(eta * h0)
            # tail: ratio on the last grid points (u in (6, 8])
            tail = [fl(eta * h_grid[i] / g_grid[i]) for i in range(len(GU)) if GUf[i] > 6]
            num = eta * W - lam * (1 + max(fl(eta * h0 - 1), 0))
            wlb = num / supw
            R1 = max(abs(fl(a - logn)), abs(fl(b - logn)))
            useful = (not extra) and rho_j_max < 1 and rho0 < 1 and rmax < 1
            others = [mm for mm in range(2, 200) if mm != nk and fl(a) <= math.log(mm) <= fl(b)]
            E["other_logs_in_S"] = others
            E["contains_logn"] = bool(fl(a) <= math.log(nk) <= fl(b))
            E.update(a=astr(a, 40), b=astr(b, 40), a_minus_logn=astr(a - logn, 8), b_minus_logn=astr(b - logn, 8),
                     a_minus_uk=astr(a - uk, 8), b_minus_uk=astr(b - uk, 8), radius_logn=R1,
                     quad_radius=fl(r0), sup_w=astr(supw, 10), sup_w_over_eta=fl(supw / eta), w2_on_ab=conc,
                     extra_components_grid=extra, rho_other_nodes_max=rho_j_max, rho_at_0=rho0,
                     rho_max_away=rmax, rho_max_away_at=rmax_u, rho_tail_max=max(tail), rho_at_8=tail[-1],
                     weight_lb=astr(wlb, 12), weight_lb_over_mu=fl(wlb / mu[k]),
                     status=("useful: S = one interval (a, b)" if useful else "NOT single interval / check"))
            rec["etas"].append(E)
        out["nodes"].append(rec)
        e2 = rec["etas"][0]
        e10 = rec["etas"][1]

        def part(E):
            if "a" not in E:
                return "%-46s" % E["status"]
            return "R=%.3e lb/mu=%.3f extra=%d rhomax=%.1e rho_j=%.1e%s" % (
                E["radius_logn"], E["weight_lb_over_mu"], len(E["extra_components_grid"]), E["rho_max_away"],
                E["rho_other_nodes_max"], (" logs%s" % E["other_logs_in_S"]) if E["other_logs_in_S"] else "")
        def verdict(E):
            if "a" not in E:
                return "n/a"
            if not E["status"].startswith("useful"):
                return "NOT-single"
            return "isolates" if not E["other_logs_in_S"] else "single-but-wide"

        line = ("deg %2d k %2d n %3d u-logn %10.3e mu %.10f 2L/sqn %.10f g'' %.4e W %.10f kkt %.1e h0 %.3e "
                "eta0 %.4e eta1/eta0 %s | x2: %s | x10: %s | %s") % (
            2 * d, k + 1, nk, fl(uk - logn), fl(mu[k]), fl(lamn), fl(gpp[k]), fl(W), fl(W - lam * h0 - mu[k]), fl(h0),
            fl(eta0), ("%.3e" % rec["eta1_over_eta0"]) if eta1 is not None else "inf", part(e2), part(e10),
            verdict(e2) + "/" + verdict(e10))
        lines.append(line)
        print(line, flush=True)
    out["time_s"] = time.time() - t0
    json.dump(out, open("traps_d%d.json" % d, "w"), indent=1)
    open("traps_d%d.out" % d, "w").write("\n".join(lines) + "\n")


for d in map(int, sys.argv[1:]):
    run(d)
