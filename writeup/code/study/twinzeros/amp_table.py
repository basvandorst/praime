"""amp_table.py -- markdown tables of the predicted coefficients (all primitive characters of conductor <= 30)."""
import numpy as np, theory
from collections import defaultdict

chars = theory.primitive_characters(30)
out = []
for name in ("h2", "h4", "h6", "h10", "h30", "SG", "T026", "T046"):
    S, rows = theory.prediction_table(name, chars)
    forms = theory.PATTERNS[name]
    sc = theory.scales(name)
    zeta = {a: S / a * sum(1 for a2, b in forms if a2 == a) for a in sc}
    out.append("\n**%s** forms %s, S = %.6f; zeta: coefficient on psi(a x) - a x: %s\n" % (
        name, forms, S, ", ".join("a=%d: %.4f" % (a, v) for a, v in zeta.items())))
    out.append("| conductor | characters present (label: parity, order, coefficient beta_a on psi(a x, chi)) | loudness |beta|/zeta |")
    out.append("|---|---|---|")
    byq = defaultdict(list)
    absent = defaultdict(list)
    for r in rows:
        (byq if r["present"] else absent)[r["q"]].append(r)
    for q in sorted(byq):
        cells = []
        seen = set()
        for r in byq[q]:
            if r["label"] in seen:
                continue
            seen.add(r["label"]); seen.add(r.get("conj"))
            bs = "; ".join("a=%d: %s" % (a, ("%.4f" % v.real) if abs(v.imag) < 1e-9 else ("%.4f%+.4fi" % (v.real, v.imag))) for a, v in r["beta"].items() if abs(v) > 1e-12)
            pair = "" if r["real"] else " (and conjugate %s)" % r["conj"]
            cells.append("%s%s: %s, order %d, %s" % (r["label"], pair, "even" if r["parity"] == 1 else "odd", r["order"], bs))
        loud = max(abs(v) / zeta[a] for r in byq[q] for a, v in r["beta"].items() if abs(v) > 1e-12)
        out.append("| %d | %s | %.4f |" % (q, "<br>".join(cells), loud))
    ab = ", ".join("%d (%d)" % (q, len(v)) for q, v in sorted(absent.items()))
    out.append("\nAbsent (conductor (number of characters)): %s\n" % ab)
open("amp_table.md", "w").write("\n".join(out))
print("\n".join(out[:40]))
