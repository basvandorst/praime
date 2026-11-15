"""Format lines.json (checked by lines_check.json) and the question-3 JSON
files into markdown tables: lines_table.md, q3_table.md.  Formatting only."""
import json


def lines_tables():
    L = json.load(open("lines.json"))
    chk = json.load(open("lines_check.json"))
    title = {"1": "quadform P = 1", "2": "quadform P = 2", "3": "quadform P = 3",
             "eis": "eisenstein"}
    out = []
    for k in ("1", "2", "3", "eis"):
        r = L[k]
        okset = {x["line"]: x["ok"] for x in chk[k]["lines"]}
        out.append(f"#### {title[k]}  (background B = {r['B_theory']:.4f}, "
                   f"empirical {r['B_empirical']:.4f}; {r['primes_in_R']:,} primes on "
                   f"{r['points_in_R']:,} lattice points)\n")
        hdr = ("| line | n along the line | points | primes | Bateman-Horn | obs/pred | C | C/B | rechecked |\n"
               "|---|---|---:|---:|---:|---:|---:|---:|:---:|")
        for tag, name in (("top_by_count", "5 lines with the most primes (b, a, a - b or a + b fixed)"),
                          ("top_by_contrast", "5 brightest relative to the background (same directions, >= 100 points)"),
                          ("top_by_contrast_any_direction", "brightest in any direction (steps up to 3), >= 100 points"),
                          ("named", "reference lines")):
            items = list(r[tag].values()) if tag == "named" else r[tag]
            out.append(f"*{name}*\n")
            out.append(hdr)
            for d in items:
                rp = d["ratio_obs_pred"]
                out.append(f"| {d['line']} | {d['n(t)']} | {d['points']} | {d['primes']} | "
                           f"{d['predicted']:.1f} | {rp:.3f} | {d['C_BH']:.4f} | "
                           f"{d['contrast_C_over_B']:.2f} | {'yes' if okset.get(d['line']) else 'NO'} |")
            out.append("")
    # aggregate over the distinct lines reported (each line once per layout)
    seen, obs, pred, ratios = set(), 0, 0.0, []
    for k, r in L.items():
        for tag in ("top_by_count", "top_by_contrast", "top_by_contrast_any_direction", "named"):
            items = r[tag].values() if tag == "named" else r[tag]
            for d in items:
                if (k, d["line"]) in seen or not d["ratio_obs_pred"]:
                    continue
                seen.add((k, d["line"]))
                obs += d["primes"]
                pred += d["predicted"]
                ratios.append(d["ratio_obs_pred"])
    out.append(f"Over the {len(ratios)} distinct lines above: {obs:,} primes observed, "
               f"{pred:,.1f} predicted by Bateman-Horn (ratio {obs / pred:.3f}); "
               f"per line obs/pred between {min(ratios):.3f} and {max(ratios):.3f}.")
    open("lines_table.md", "w").write("\n".join(out) + "\n")


def q3_table():
    a = json.load(open("qf_extra.json"))
    b = {int(k): v for k, v in json.load(open("qf_extra_check.json")).items()}
    rows = ["| P | h(-4P) | criterion (method A, qf_extra.py) | primes checked | agree | exceptions | "
            "method B (qf_extra_check.py): class polynomial H_{-4P} | method B: simple polynomial |",
            "|---:|---:|---|---:|---:|---|---|---|"]
    spec = [
        (27, 3, "P=27 cubic residue", "p = 1 (mod 3) and 2^((p-1)/3) = 1 (mod p)"),
        (64, 4, "P=64 biquadratic residue", "p = 1 (mod 4) and 2^((p-1)/4) = 1 (mod p)"),
        (14, 4, "P=14 (x^2+1)^2=8", "(-14/p) = 1 and (x^2+1)^2 = 8 (mod p) solvable"),
        (11, 3, "P=11 x^3-x^2+x+1 three roots", "x^3 - x^2 + x + 1 has three roots mod p"),
        (23, 3, "P=23 x^3-x-1 three roots", "x^3 - x - 1 has three roots mod p"),
    ]
    for P, h, key, text in spec:
        r = a[key]
        c = b[P]
        sp = "; ".join(f"p={p}: drawn={'yes' if rep else 'no'}, criterion={'yes' if cr else 'no'}"
                       for p, rep, cr in r["special"])
        exc = (f"{r['disagree']}" + (f" (set aside: {sp})" if sp else ""))
        rows.append(
            f"| {P} | {h} | {text} | {r['checked']:,} | {r['agree']:,} | {exc} | "
            f"{c['classpoly_agree']:,} of {c['classpoly_checked']:,} agree "
            f"(skipped p dividing 4P*disc H: {c['classpoly_excluded_primes']}) | "
            f"{c['simple_criterion']}: {c['simple_agree']:,} of {c['simple_checked']:,} agree"
            f"{', exceptions ' + str(c['simple_exceptions']) if c['simple_exceptions'] else ''}; "
            f"same field as polredabs(H) = {c['polredabs(H)']}: {c['same_field_as_simple']} |")
    extra = [k for k in a if k not in [s[2] for s in spec]]
    rows.append("")
    rows.append("Further consistency rows of qf_extra.py:\n")
    rows.append("| statement | primes checked | agree | disagree |")
    rows.append("|---|---:|---:|---:|")
    for k in extra:
        r = a[k]
        if "checked" in r:
            rows.append(f"| {k} | {r['checked']:,} | {r['agree']:,} | {r['disagree']} |")
        else:
            rows.append(f"| {k} | {r} | | |")
    open("q3_table.md", "w").write("\n".join(rows) + "\n")


if __name__ == "__main__":
    lines_tables()
    q3_table()
    print(open("q3_table.md").read())
