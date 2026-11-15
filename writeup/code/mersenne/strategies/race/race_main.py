"""Layout 2 (chebrace): the squares-vs-non-squares prime race mod P, P = 3..40.

R = squares among the units mod P, Nn = the other units.
h(x) = |R| * #{p <= x, p mod P in Nn} - |Nn| * #{p <= x, p mod P in R}   (primes p | P ignored).
Prediction: E[h(x)] = (1/2)|Nn| Li(sqrt x), Li(y) = li(y) - li(2).
For X = 10^7 and X = 5*10^5: average of h/((1/2)|Nn|Li(sqrt x)) over ln x in [10^3, X]
(geometric grid, 10^5 points), share of primes p <= X after which h < 0, first such p,
logarithmic density of {h < 0} on the grid, number of sign changes.
Output: race_main.json, race_main.txt.
"""
import json
import math
import numpy as np
from scipy.special import expi

XMAX = 10 ** 7
GRID = 100_000


def sieve(n):
    s = np.ones(n + 1, bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return s


primes = np.nonzero(sieve(XMAX))[0].astype(np.int64)
LI2 = expi(math.log(2.0))


def Li(y):
    return expi(np.log(y)) - LI2


out = []


def say(s=""):
    print(s)
    out.append(s)


say(f"pi(10^7) = {len(primes)}, pi(5*10^5) = {int((primes <= 500000).sum())}")
res = {}
for P in range(3, 41):
    units = [a for a in range(P) if math.gcd(a, P) == 1]
    phi = len(units)
    R = sorted({b * b % P for b in units})
    Nn = [a for a in units if a not in R]
    roots = {a: sum(1 for b in units if b * b % P == a) for a in units}
    r = roots[1]
    assert all(roots[a] == r for a in R) and all(roots[a] == 0 for a in Nn)
    assert sum(roots[a] for a in R) == phi and len(R) * r == phi
    step_tab = np.zeros(P, np.int64)
    step_tab[R] = -len(Nn)
    step_tab[Nn] = len(R)
    h = np.cumsum(step_tab[primes % P])      # h after each prime
    entry = {"phi": phi, "R": len(R), "N": len(Nn), "r": r, "Rset": R}
    # Riemann prime-power count Pi*(y) = sum_{p^k <= y, p !| P} 1/k, for the refined prediction (1/2)|Nn| Pi*(sqrt x)
    pp = sorted((int(p) ** k, 1.0 / k) for p in primes[primes <= 4000] if P % int(p)
                for k in range(1, 13) if int(p) ** k <= 4000)
    ppv = np.array([0] + [a for a, _ in pp], np.int64)
    ppw = np.cumsum([0.0] + [w for _, w in pp])
    for X in (XMAX, 500_000):
        sel = primes <= X
        hp = h[sel]
        neg = hp < 0
        first = int(primes[sel][np.argmax(neg)]) if neg.any() else None
        neg_big = neg & (primes[sel] > 1000)
        first_big = int(primes[sel][np.argmax(neg_big)]) if neg_big.any() else None
        last = int(primes[sel][np.nonzero(neg)[0][-1]]) if neg.any() else None
        # sign changes of h (ignoring zeros): count switches between h>0 and h<0 regimes
        nz = hp[hp != 0]
        changes = int((np.sign(nz[1:]) != np.sign(nz[:-1])).sum())
        xs = np.geomspace(1e3, X, GRID)
        idx = np.searchsorted(primes, xs, side="right") - 1
        hx = h[idx]
        ratio = hx / (0.5 * len(Nn) * Li(np.sqrt(xs)))
        PiS = ppw[np.searchsorted(ppv, np.sqrt(xs), side="right") - 1]
        ratio_pi = hx / (0.5 * len(Nn) * PiS)
        entry[str(X)] = {
            "avg_ratio": float(ratio.mean()),
            "avg_ratio_Pi": float(ratio_pi.mean()),
            "avg_PiS_over_Li": float((PiS / Li(np.sqrt(xs))).mean()),
            "share_primes_neg": float(neg.mean()),
            "n_primes_neg": int(neg.sum()),
            "first_neg": first,
            "first_neg_after_1000": first_big,
            "last_neg": last,
            "logdens_neg": float((hx < 0).mean()),
            "min_h": int(hp.min()),
            "h_end": int(hp[-1]),
            "pred_end": float(0.5 * len(Nn) * Li(math.sqrt(X))),
            "sign_changes": changes,
        }
    res[P] = entry

say()
say("X = 10^7 (ratio averaged over ln x in [10^3, 10^7])")
say(" P phi |R| |N|  r  avg h/pred  h(X)    pred(X)  share<0    #p<0   logdens<0  first h<0  first>1000   last h<0  changes")
for P, e in res.items():
    d = e[str(XMAX)]
    say(f"{P:2d} {e['phi']:3d} {e['R']:3d} {e['N']:3d} {e['r']:2d}  {d['avg_ratio']:8.4f}  {d['h_end']:6d} {d['pred_end']:9.1f}  "
        f"{d['share_primes_neg']:.5f} {d['n_primes_neg']:7d}  {d['logdens_neg']:.5f}  {str(d['first_neg']):>9s}  "
        f"{str(d['first_neg_after_1000']):>9s}  {str(d['last_neg']):>9s}  {d['sign_changes']:5d}")
say()
say("X = 5*10^5 (ratio averaged over ln x in [10^3, 5*10^5])")
say(" P  avg h/pred  h(X)   pred(X)  share<0  #p<0  logdens<0  first h<0  first>1000  last h<0")
for P, e in res.items():
    d = e["500000"]
    say(f"{P:2d}  {d['avg_ratio']:8.4f}  {d['h_end']:5d} {d['pred_end']:8.1f}  {d['share_primes_neg']:.5f} {d['n_primes_neg']:5d}  "
        f"{d['logdens_neg']:.5f}  {str(d['first_neg']):>9s}  {str(d['first_neg_after_1000']):>9s}  {str(d['last_neg']):>9s}")

say()
div24 = [P for P in res if 24 % P == 0]
say(f"P | 24 in range: {div24}; R = {{1}} for exactly these P: "
    f"{sorted(P for P, e in res.items() if e['R'] == 1) == div24}")
nev = [P for P, e in res.items() if e[str(XMAX)]['first_neg'] is None]
say(f"P with h >= 0 at every prime <= 10^7: {nev}")
nev5 = [P for P, e in res.items() if e['500000']['first_neg'] is None]
say(f"P with h >= 0 at every prime <= 5*10^5: {nev5} ({len(nev5)} moduli)")
rank = sorted(res, key=lambda P: res[P][str(XMAX)]['share_primes_neg'])
say("P ranked by share of primes with h < 0 (10^7): " +
    ", ".join(f"{P}:{res[P][str(XMAX)]['share_primes_neg']:.4f}" for P in rank))
for X in (XMAX, 500_000):
    allr = np.array([res[P][str(X)]['avg_ratio'] for P in res])
    allp = np.array([res[P][str(X)]['avg_ratio_Pi'] for P in res])
    alls = np.array([res[P][str(X)]['avg_PiS_over_Li'] for P in res])
    say(f"X = {X}: avg of h/((1/2)|N|Li(sqrt x)) over all P: mean {allr.mean():.4f}, median {np.median(allr):.4f}, "
        f"min {allr.min():.4f}, max {allr.max():.4f}")
    say(f"          avg of h/((1/2)|N|Pi*(sqrt x)):           mean {allp.mean():.4f}, median {np.median(allp):.4f}, "
        f"min {allp.min():.4f}, max {allp.max():.4f}")
    say(f"          avg of Pi*(sqrt x)/Li(sqrt x):              mean {alls.mean():.4f}, min {alls.min():.4f}, max {alls.max():.4f}")
say()
say("Refined ratio h/((1/2)|N|Pi*(sqrt x)) per P (10^7 | 5*10^5):")
say("  " + "; ".join(f"{P}: {res[P][str(XMAX)]['avg_ratio_Pi']:.3f} | {res[P]['500000']['avg_ratio_Pi']:.3f}" for P in res))

json.dump({str(k): v for k, v in res.items()}, open("race_main.json", "w"), indent=1)
open("race_main.txt", "w").write("\n".join(out) + "\n")
