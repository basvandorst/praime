#!/usr/bin/env python3
"""Tables from results.json: calibration of the nulls, the six frozen hypotheses, the
positive control, the broad feature family (discovery -> holdout), and power."""
import json, math, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "results.json")))
T = R["tests"]
fam = R["families"]
EG = math.exp(0.5772156649015329)

print("== data ==")
print(f"families {len(fam)} (discovery {len(R['disc'])}, holdout {len(R['hold'])})")
print(f"universe n (log2 R >= 64): {sum(f['n'] for f in fam)}, PRP hits {sum(f['O'] for f in fam)}, "
      f"sieve survivors {sum(f['surv'] for f in fam)}")
print(f"null S calibration: observed {sum(f['O'] for f in fam)} vs heuristic e^gamma ln Q/ln R "
      f"{sum(f['thS'] for f in fam):.1f}  (ratio {sum(f['O'] for f in fam) / sum(f['thS'] for f in fam):.3f})")
print("discovery:", " ".join(R["disc"]))
print("holdout:  ", " ".join(R["hold"]))


def row(name, nulls=("G", "C", "S")):
    out = []
    for N in nulls:
        oD, eD, puD, plD = T[name][N + "D"]; oH, eH, puH, plH = T[name][N + "H"]
        out.append(f"{N}: D {oD:.0f}/{eD:.1f}={oD / max(eD, 1e-9):.2f} p+{puD:.3f} | "
                   f"H {oH:.0f}/{eH:.1f}={oH / max(eH, 1e-9):.2f} p+{puH:.3f}")
    return out


print("\n== calibration: residue classes under each null (pooled all families, O/E) ==")
for m in (3, 4, 8):
    for r in range(1, m):
        if math.gcd(r, m) != 1:
            continue
        nm = f"n={r} mod {m}"
        s = []
        for N in ("G", "C", "S"):
            o = T[nm][N + "D"][0] + T[nm][N + "H"][0]; e = T[nm][N + "D"][1] + T[nm][N + "H"][1]
            s.append(f"{N} {o:.0f}/{e:.1f}={o / e:.2f}")
        print(f"  {nm:12s} " + "  ".join(s))

print("\n== positive control (theorem: q = 2n+1 prime, q not dividing b-/+1, (+-b|q) = 1  =>  q | R) ==")
nm = "P0 control: 2n+1 prime, (+-b|2n+1)=1"
for N in ("G", "C", "S"):
    oD, eD, puD, plD = T[nm][N + "D"]; oH, eH, puH, plH = T[nm][N + "H"]
    print(f"  {N}: D {oD:.0f}/{eD:.2f} p- {plD:.2e} | H {oH:.0f}/{eH:.2f} p- {plH:.2e}")

print("\n== six frozen hypotheses (one-sided, excess) ==")
for nm in [k for k in T if k.startswith("H")]:
    print(nm)
    for line in row(nm):
        print("    " + line)

print("\n== broad family: residues n = r mod m (3 <= m <= 60) and digit statistics ==")
broad = [k for k in T if k.startswith("n=") or k.startswith("dig ")]
for N in ("G", "C", "S"):
    pD = np.array([2 * min(T[k][N + "D"][2], T[k][N + "D"][3]) for k in broad]).clip(max=1)
    pH = np.array([2 * min(T[k][N + "H"][2], T[k][N + "H"][3]) for k in broad]).clip(max=1)
    sel = [k for k, p in zip(broad, pD) if p < 0.01]
    rep = []
    for k in sel:
        oD, eD, puD, plD = T[k][N + "D"]; oH, eH, puH, plH = T[k][N + "H"]
        up = oD > eD
        ph = puH if up else plH
        rep.append((k, oD, eD, oH, eH, ph))
    nrep = sum(1 for r in rep if r[5] < 0.05 / max(len(sel), 1))
    print(f"null {N}: {len(broad)} features; p_D<0.05: {np.mean(pD < 0.05):.3f} of them, p_H<0.05: "
          f"{np.mean(pH < 0.05):.3f}; selected on D (p<0.01): {len(sel)}; replicate in H "
          f"(same direction, p < 0.05/{max(len(sel), 1)}): {nrep}")
    for k, oD, eD, oH, eH, ph in sorted(rep, key=lambda r: r[5])[:8]:
        print(f"    {k:28s} D {oD:.0f}/{eD:.1f}={oD / eD:.2f}  H {oH:.0f}/{eH:.1f}={oH / max(eH, 1e-9):.2f}  p_H {ph:.3g}")
    # best holdout replications regardless of selection
print("\n== all 79 families pooled (one-sided excess p) ==")
import numpy as np
for nm in [k for k in T if k.startswith("H")]:
    s = []
    for N in ("G", "C", "S"):
        s.append(f"{N} {T[nm][N + 'A'][0]:.0f}/{T[nm][N + 'A'][1]:.1f}={T[nm][N + 'A'][0] / T[nm][N + 'A'][1]:.2f} p+{T[nm][N + 'A'][2]:.3f}")
    print(f"  {nm:42s} " + "  ".join(s))
print("\n== power (holdout, null C): smallest O/E found at one-sided 5% with 80% power ==")
for nm in [k for k in T if k.startswith("H")]:
    eH = T[nm]["CH"][1]
    print(f"  {nm:42s} E_H {eH:6.1f}  detectable O/E ~ {1 + 2.49 / math.sqrt(eH):.2f}")
