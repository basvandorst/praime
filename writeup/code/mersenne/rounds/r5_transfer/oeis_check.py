#!/usr/bin/env python3
"""Cross-check data.jsonl against the local OEIS copy (names.gz, stripped.gz).
Finds sequences named 'Numbers k such that (b^k - 1)/(b-1) is prime' (and the +1 form),
plus A004023 (base 10) and A000978 (Wagstaff), and compares the odd terms <= N_b."""
import gzip, json, os, re

HERE = os.path.dirname(os.path.abspath(__file__))
OE = "/home/bas/praime/references/oeis"
pat = re.compile(r"^(A\d+) Numbers [kn] such that \((\d+)\^[kn] ?([-+]) ?1\)/(\d+) is (?:a )?prime\.?\s*$")
ids = {}
with gzip.open(os.path.join(OE, "names.gz"), "rt", errors="replace") as fh:
    for line in fh:
        m = pat.match(line.strip())
        if m:
            a, b, s, d = m.group(1), int(m.group(2)), m.group(3), int(m.group(4))
            kind = "m" if s == "-" else "p"
            if d == (b - 1 if kind == "m" else b + 1):
                ids.setdefault((kind, b), a)
ids[("m", 10)] = "A004023"; ids[("p", 2)] = "A000978"
want = set(ids.values())
terms = {}
with gzip.open(os.path.join(OE, "stripped.gz"), "rt") as fh:
    for line in fh:
        if line[:7] in want:
            terms[line[:7]] = [int(x) for x in line[8:].strip().strip(",").split(",") if x]
ok = bad = 0
for line in open(os.path.join(HERE, "data.jsonl")):
    d = json.loads(line)
    key = (d["kind"], d["b"]); N = d["N"]
    ours = sorted(r[0] for r in d["recs"] if r[1] == "PRP")
    if key not in ids:
        continue
    a = ids[key]
    t = terms.get(a, [])
    if not t:
        print(f"{d['kind']}{d['b']}: {a} no terms"); continue
    covered = max(t)
    theirs = sorted(x for x in t if x % 2 == 1 and x <= N and x > 2)
    cmp_ours = [x for x in ours if x <= covered]
    same = theirs == cmp_ours
    ok += same; bad += not same
    print(f"{d['kind']}{d['b']:>3} {a}: OEIS data to {covered}, compared n <= {min(N, covered)}: "
          f"{len(theirs)} terms, {'AGREE' if same else 'DISAGREE ours=' + str(cmp_ours) + ' oeis=' + str(theirs)}")
print(f"families compared {ok + bad}: agree {ok}, disagree {bad}")
