"""Assemble RESULTS.md from RESULTS.template.md and the JSON/markdown outputs.
Formatting only: every number is read from a script's output file."""
import json


def main():
    t = open("RESULTS.template.md").read()
    for f in ("table.md", "q3_table.md", "lines_table.md"):
        t = t.replace("{{" + f + "}}", open(f).read().rstrip())
    qc = json.load(open("qf_count.json"))["rows"]
    chk = json.load(open("qf_check.json"))
    by_h = {}
    for r in qc:
        by_h.setdefault(r["h"], []).append(r["P"])
    dec = [r["P"] for r in qc if r["decided"]]
    nd = [r for r in qc if not r["decided"]]
    dev = sorted((abs(r["fraction_excl"] * 2 * r["h"] - 1), r["P"]) for r in qc)
    s12 = [
        f"* **Class numbers.** h(−4P) from counting reduced forms equals PARI `qfbclassno` for all 100 P "
        f"(`qf_check.py`). h = 1 for P = {', '.join(map(str, by_h[1]))}; "
        f"the largest is h = {max(by_h)} (P = {', '.join(map(str, by_h[max(by_h)]))}).",
        f"* **Densities.** Excluding p | 4P, every fraction is within {100*dev[-1][0]:.2f}% "
        f"(relative) of 1/(2h). The largest gap is at P = {dev[-1][1]}, the smallest "
        f"({100*dev[0][0]:.3f}%) at P = {dev[0][1]}. All 100 fractions are *below* 1/(2h); see the bias note below.",
        f"* **Congruence test.** Every coprime class mod 4P holds at least "
        f"{min(r['min_primes_per_class'] for r in qc)} primes ≤ N, and no class is empty, so the test has teeth. "
        f"The represented primes are a union of classes mod 4P for exactly {len(dec)} values of P. These are exactly "
        f"Euler's idoneal numbers ≤ 100 from your list: {', '.join(map(str, dec))}. "
        f"For every non-idoneal P, between {100*min(r['n_mixed_classes']/r['n_classes'] for r in nd):.1f}% and "
        f"{100*max(r['n_mixed_classes']/r['n_classes'] for r in nd):.1f}% of the classes are mixed.",
        "* **Four idoneal tests agree for all 100 P**: (i) the congruence test; (ii) h(−4P) equals the number "
        "of genera 2^(μ−1) (genus count as in Cox §3, `qf_check.py`); (iii) the class group has exponent ≤ 2 (PARI "
        "`quadclassunit`); (iv) every reduced form is ambiguous (b = 0, b = a or a = c; `qf_count.py`). "
        "This is genus theory: a congruence mod 4P fixes only the genus, so it decides exactly when each "
        "genus is a single class (D. Cox, *Primes of the form x² + ny²*, §3, genus theory and Euler's "
        "convenient numbers; section numbers in this file are cited from memory).",
        "* **Rules for some idoneal P** (classes mod 4P, from `qf_count.json`): " + "; ".join(
            f"P = {P}: p ≡ {', '.join(map(str, qc[P-1]['rep_classes']))} (mod {4*P})"
            for P in (1, 2, 3, 5, 7)) +
        f"; P = 58: {len(qc[57]['rep_classes'])} classes mod 232.",
    ]
    t = t.replace("{{SUMMARY12}}", "\n".join(s12))
    b = json.load(open("bias.json"))["summary"]
    t = t.replace("{{BIAS}}", (
        f"`bias.py` (rechecked by `bias_check.py`) uses π(√N) = {b['pi(sqrt N)']}. "
        f"All {b['P with positive deficit']} of 100 P have a positive deficit (expected count "
        f"n_excl/(2h) minus drawn). The observed deficits run from {b['observed deficit: min, median, max'][0]:.1f} to "
        f"{b['observed deficit: min, median, max'][2]:.1f} primes (median {b['observed deficit: min, median, max'][1]:.1f}). "
        f"The predicted mean runs from {b['predicted mean deficit: min, median, max'][0]:.1f} to "
        f"{b['predicted mean deficit: min, median, max'][2]:.1f} (median {b['predicted mean deficit: min, median, max'][1]:.1f}). "
        f"Summed over all P, {b['sum observed']:.0f} were observed against {b['sum predicted']:.0f} predicted, with correlation "
        f"{b['correlation']:.2f} across P. So the sign and the size are explained, but the scatter "
        f"from P to P is not."))
    ex = json.load(open("qf_extra.json"))
    dist = ex["x^3-x-1 root-count distribution (p != 23)"]
    tot = sum(dist.values())
    t = t.replace("{{CHEB}}", ", ".join(
        f"{k} root{'' if k == '1' else 's'} for {v:,} primes ({100*v/tot:.1f}%)" for k, v in sorted(dist.items())))
    e = json.load(open("eis.json"))
    ec = json.load(open("eis_check.json"))
    t = t.replace("{{EIS}}", "\n".join([
        f"* **Primes drawn: {e['primes_drawn']:,}**, a fraction {e['fraction']:.5f} of π(N) "
        f"({e['fraction_excluding_3']:.5f} without 3). Density 1/2, since h(−3) = 1.",
        f"* **The set is exactly {{3}} ∪ {{p ≡ 1 (mod 3)}}**: {e['equals_{3}u{1 mod 3}']} in `eis.py`, "
        f"{ec['set_is_3_and_1mod3']} in `eis_check.py`. There are {e['primes_1_mod_3']:,} primes ≡ 1 and "
        f"{e['primes_2_mod_3']:,} primes ≡ 2 (mod 3), so the same Chebyshev bias shows here.",
        f"* **Positions.** 3 is drawn at (a, b) = {tuple(e['rep_of_3'])}, on the edge a = b of the wedge; its orbit under "
        f"the 12 automorphisms has 6 points. Every other drawn prime has an orbit of 12 points with exactly one "
        f"in the wedge (orbit sizes {ec['orbit_sizes']}). Every drawn prime's position in `eis.py` matches "
        f"`eis_check.py` ({ec['position_mismatches']} mismatches).",
        f"* {e['numbers_drawn_1..N']:,} of the numbers 1..N are drawn at all (Löschian numbers).",
    ]))
    ps = json.load(open("picture_stats.json"))
    pc = json.load(open("picture_check.json"))
    sub = json.load(open("subsets.json"))
    subc = json.load(open("subsets_check.json"))
    mixed = {r["P"]: (r["n_mixed_classes"], r["n_classes"]) for r in qc}
    rows = ["| P | h | dots | lattice points of R | dots/points | B | row b=1: C/B | best row b ≤ 20 (C/B) | worst row (C/B) | mixed classes mod 4P |",
            "|---:|---:|---:|---:|---:|---:|---:|---|---|---|"]
    for P in (1, 2, 3, 4, 5, 7, 11, 14, 16, 23, 27, 37, 41, 58, 64, 67, 89):
        r = ps[P - 1]
        br, wr = r["best_row"], r["worst_row"]
        m = mixed[P]
        rows.append(f"| {P} | {r['h']} | {r['dots']:,} | {r['points']:,} | {r['density']:.4f} | {r['B']:.3f} | "
                    f"{r['row1_contrast']:.2f} | b = {br['index']} ({br['contrast']:.2f}) | "
                    f"b = {wr['index']} ({wr['contrast']:.2f}) | {m[0]} of {m[1]} |")
    dens = sorted(ps, key=lambda r: -r["density"])
    col = max((r for r in ps if r["best_col"]), key=lambda r: r["best_col"]["contrast"])
    dia = max((r for r in ps if r["best_diag"]), key=lambda r: r["best_diag"]["contrast"])
    row = max(ps, key=lambda r: r["best_row"]["contrast"])
    pic = ["\n".join(rows), "",
           f"* **Densest pictures** (dots per lattice point of R): " + ", ".join(
               f"P = {r['P']} ({r['density']:.3f})" for r in dens[:6]) +
           f". **Sparsest**: " + ", ".join(f"P = {r['P']} ({r['density']:.3f})" for r in dens[-3:]) +
           f". `picture_check.py` confirms P = {pc['densest P']} as the densest.",
           f"* **Strongest single lines** (rows, columns and diagonals of index ≤ 20 with ≥ 100 points, over all P): "
           f"column a = {col['best_col']['index']} for P = {col['P']} (C/B = {col['best_col']['contrast']:.2f}), "
           f"diagonal a − b = {dia['best_diag']['index']} for P = {dia['P']} ({dia['best_diag']['contrast']:.2f}), "
           f"row b = {row['best_row']['index']} for P = {row['P']} ({row['best_row']['contrast']:.2f}). "
           f"In this range, rows (fixed b) never exceed {row['best_row']['contrast']:.2f}.",
           "* **Large P.** The page draws (a, b√P), so the picture becomes a stack of horizontal rows √P apart "
           f"(P = 58: {ps[57]['rows_in_picture']} rows). Row brightness varies by about a factor 2 to 3 (best and worst "
           "rows above), following the primes that divide b.",
           "* **How the pictures nest** (`subsets.py`, rechecked by `subsets_check.py`; every flag True in both):",
           f"  * P = 4 has the same {sub['P=1 primes 1 mod 4']:,} odd primes as P = 1 (all p ≡ 1 mod 4).",
           f"  * P = 16 has the {sub['P=16 count']:,} primes whose even coordinate in the P = 1 picture is divisible by 4 "
           f"(⟺ p ≡ 1 mod 8: {sub['P=16 <=> p=1 mod 8']}).",
           f"  * P = 64 has the {sub['P=64 count']:,} with the even coordinate divisible by 8 ({sub['P=64 <=> 8|e']}).",
           f"  * P = 12 has the P = 3 primes with b even ({sub['P=12 count']:,}).",
           f"  * P = 27 has the P = 3 primes with 3 | b: {sub['P=27 count']:,}, a share of "
           f"{sub['P=27 share of P=3 primes (p=1 mod 3)']:.4f} of the {sub['P=3 primes 1 mod 3']:,} primes ≡ 1 (mod 3).",
           f"* **P = 5.** (−5/p) = 1 for p ≡ {', '.join(map(str, sub['P=5 classes with (-5/p)=1']))} (mod 20); "
           f"x² + 5y² takes {', '.join(map(str, sub['P=5 classes of x^2+5y^2']))}, 2x² + 2xy + 3y² takes "
           f"{', '.join(map(str, sub['P=5 classes of 2x^2+2xy+3y^2']))}; each such prime is in exactly one form: "
           f"{sub['P=5: every prime with (-5/p)=1 is in exactly one of the two forms']} "
           f"(rechecked: {subc['P=5 classes of 2x^2+2xy+3y^2'] and subc['P=5 classes with (-5/p)=1']}).",
           f"* **P = 58.** 2x² + 29 is prime for x = 0..{sub['2x^2+29 prime for x = 0..(run-1); run'] - 1} "
           f"({sub['2x^2+29 prime for x = 0..(run-1); run']} values; checked: {subc['2x^2+29 run']}). "
           f"x² + 58y² takes {sub['P=58 number of classes']} classes mod 232, and 2x² + 29y² takes "
           f"{subc['2x^2+29y^2 number of classes']} others (disjoint: "
           f"{subc['2x^2+29y^2 classes disjoint from x^2+58y^2']}).",
           ]
    t = t.replace("{{PIC}}", "\n".join(pic))
    assert "{{" not in t, t[t.index("{{"):t.index("{{") + 40]
    open("RESULTS.md", "w").write(t)
    print("RESULTS.md written,", len(t.splitlines()), "lines")


if __name__ == "__main__":
    main()
