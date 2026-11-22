"""Round 15 (preregister/round15.txt): the joint law of three consecutive gaps (four-point correlations)."""
import os, json, math, time
import numpy as np
import theory2 as T2
import theory_r3 as R3
import common as C

OUT = R3.OUT
HS = list(range(2, 19, 2))
rng = np.random.default_rng(15)


def main():
    t0 = time.time()
    pats = {}
    for a in HS:
        for b in HS:
            for c in HS:
                B = (0, a, a + b, a + b + c)
                pats[(a, b, c)] = T2.Pattern(list(B), 4)
    print('patterns built', round(time.time() - t0), flush=True)
    out = {}
    for win, (L, W) in R3.WIN.items():
        u = np.linspace(math.log(L), math.log(L + W), 9)
        wts = np.exp(u) / u
        P = np.zeros((9, 9, 9))
        for (a, b, c), pt in pats.items():
            if pt.SB == 0:
                continue
            vals = np.array([pt.SB / lam ** 3 * math.exp(pt.log_avoid(lam, 4)) for lam in u])
            P[a // 2 - 1, b // 2 - 1, c // 2 - 1] = np.trapezoid(wts * vals, u) / np.trapezoid(wts, u)
        M = C.primes_mask(L, W); pr = np.nonzero(M)[0].astype(np.int64) + L
        g = np.diff(pr)
        n3 = len(g) - 2
        sel = (g[:-2] <= 18) & (g[1:-1] <= 18) & (g[2:] <= 18)
        idx = ((g[:-2][sel] // 2 - 1) * 9 + (g[1:-1][sel] // 2 - 1)) * 9 + (g[2:][sel] // 2 - 1)
        obs = np.bincount(idx, minlength=729).astype(float)
        exp = P.ravel() * n3

        def chi(o, e):
            big = e >= 5
            o2 = np.append(o[big], o[~big].sum()); e2 = np.append(e[big], e[~big].sum())
            if e2[-1] < 5:
                o2, e2 = o2[:-1], e2[:-1]
            c2 = float(np.sum((o2 - e2) ** 2 / e2)); d = len(o2) - 1
            return c2, d
        c2, d = chi(obs, exp)
        ctrl = []
        for _ in range(10):
            s = rng.multinomial(int(obs.sum()), exp / exp.sum()).astype(float)
            cc, dd = chi(s, exp * obs.sum() / exp.sum()); ctrl.append(cc / dd)
        out[win] = {'triples_all_le18_data': float(sel.mean()), 'law': float(P.sum()), 'chi2': c2, 'dof': d,
                    'chi2_per_dof': c2 / d, 'z': (c2 / d - 1) / math.sqrt(2 / d),
                    'control_mean': float(np.mean(ctrl)), 'control_sd': float(np.std(ctrl, ddof=1))}
        print(win, {k: round(v, 4) for k, v in out[win].items()}, flush=True)
    json.dump(out, open(os.path.join(OUT, 'r15.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
