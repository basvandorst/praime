"""Write table.md (questions 1-2) from qf_count.json, cross-referenced with
qf_check.json (Cornacchia/PARI), and print the summary numbers used in
RESULTS.md.  No new computation of primes here: only formatting and the
agreement test between the two JSON files.
"""
import json


def main():
    a = json.load(open("qf_count.json"))
    b = {r["P"]: r for r in json.load(open("qf_check.json"))}
    lines = ["| P | h(-4P) | primes drawn | / 41,538 | excl. p mid 4P | 1/(2h) | congruence-decided | idoneal |",
             "|---:|---:|---:|---:|---:|---:|:---:|:---:|"]
    lines[0] = lines[0].replace("excl. p mid 4P", "share among p not dividing 4P")
    agree = 0
    for r in a["rows"]:
        P = r["P"]
        c = b[P]
        ok = (c["drawn_all"] == r["drawn_all"] and c["h_pari"] == r["h"]
              and c["decided"] == r["decided"] and c["idoneal_genus"] == r["idoneal_list"])
        agree += ok
        lines.append(
            f"| {P} | {r['h']} | {r['drawn_all']:,} | {r['fraction_all']:.4f} | "
            f"{r['fraction_excl']:.4f} | {r['chebotarev']:.4f} | "
            f"{'yes' if r['decided'] else 'no'} | {'yes' if r['idoneal_list'] else 'no'} |")
    open("table.md", "w").write("\n".join(lines) + "\n")
    rows = a["rows"]
    print("rows agreeing with qf_check.json:", agree, "of", len(rows))
    print("decided == idoneal for all P:", all(r["decided"] == r["idoneal_list"] for r in rows))
    print("number decided:", sum(r["decided"] for r in rows))
    print("min primes in a coprime class (over all P):", min(r["min_primes_per_class"] for r in rows))
    print("empty classes anywhere:", sum(r["n_empty_classes"] for r in rows))
    nd = [r for r in rows if not r["decided"]]
    print("non-idoneal: mixed classes share min/max:",
          min(r["n_mixed_classes"] / r["n_classes"] for r in nd),
          max(r["n_mixed_classes"] / r["n_classes"] for r in nd))
    print("all reduced forms ambiguous == idoneal:",
          all(r["all_reduced_forms_ambiguous"] == r["idoneal_list"] for r in rows))
    dev = [(abs(r["fraction_excl"] * 2 * r["h"] - 1), r["P"]) for r in rows]
    print("max relative deviation |frac_excl*2h - 1|:", max(dev))
    print("special primes drawn (p | 4P):",
          {r["P"]: r["special_drawn"] for r in rows if r["special_drawn"]})
    for P in (1, 2, 3, 5, 7, 58):
        r = rows[P - 1]
        print(f"P={P} classes mod {4*P}:", r["rep_classes"])
    hs = {}
    for r in rows:
        hs.setdefault(r["h"], []).append(r["P"])
    print("P by h:", dict(sorted(hs.items())))


if __name__ == "__main__":
    main()
