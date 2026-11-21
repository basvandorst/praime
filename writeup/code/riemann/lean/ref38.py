"""Independent recount of every number on prime38.html, written separately from graph_scan.py/page_data.py.

Reads the archive through `cat parts | ` a subprocess pipe, counts lines with splitlines(),
parses imports with a token walk, computes reachability by BFS and layers with Kahn's algorithm.
Usage: python3 -I ref38.py PARTS_DIR  -> ref38.json
"""
import sys, os, glob, subprocess, tarfile, json, collections, hashlib

parts = sorted(glob.glob(os.path.join(sys.argv[1], "lean.tar.gz.part*")))
assert len(parts) == 14
h = hashlib.sha256()
for p in parts:
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)

proc = subprocess.Popen(["cat"] + parts, stdout=subprocess.PIPE)


def strip_comments(t):
    out, i, n = [], 0, len(t)
    while i < n:
        a, b = t.find("/-", i), t.find("--", i)
        a = n if a < 0 else a
        b = n if b < 0 else b
        j = min(a, b)
        out.append(t[i:j])
        if j == n:
            break
        if j == b:
            k = t.find("\n", j)
            i = n if k < 0 else k
            continue
        depth, i = 1, j + 2
        while depth and i < n:
            o, c = t.find("/-", i), t.find("-/", i)
            if c < 0:
                i = n; break
            if 0 <= o < c:
                depth += 1; i = o + 2
            else:
                depth -= 1; i = c + 2
    return "".join(out)


def imports_of(code):
    res = []
    toks = code.split()
    i = 0
    while i < len(toks):
        t = toks[i]
        if t in ("module", "public", "private", "meta", "all"):
            i += 1; continue
        if t == "import":
            i += 1
            while i < len(toks) and toks[i] not in ("import", "public", "private", "meta") and \
                    (toks[i][0].isalpha() or toks[i][0] == "_"):
                if toks[i] == "all":
                    i += 1; continue
                if i + 1 < len(toks) and toks[i + 1] in (":", ":=", "where"):
                    break
                res.append(toks[i]); i += 1
                if res and i < len(toks) and toks[i] != "import" and toks[i] not in ("public", "private", "meta"):
                    break
            continue
        break
    return res


mods = {}
theorems = 0
with tarfile.open(fileobj=proc.stdout, mode="r|gz") as tf:
    for m in tf:
        if not m.isfile() or not m.name.endswith(".lean"):
            continue
        raw = tf.extractfile(m).read()
        txt = raw.decode("utf-8", "replace")
        name = m.name.split("lean/src/", 1)[1][:-5].replace("/", ".")
        code = strip_comments(txt)
        mods[name] = (len(raw), len(txt.splitlines()), imports_of(code))
        for ln in code.splitlines():
            w = ln.split()
            while w and (w[0].startswith("@[") or w[0] in ("private", "protected", "noncomputable")):
                if w[0].startswith("@[") and not w[0].endswith("]"):
                    while w and not w[0].endswith("]"):
                        w = w[1:]
                w = w[1:]
            if w and w[0] in ("theorem", "lemma"):
                theorems += 1
proc.wait()

listing = [l.split() for l in open("lean_list.txt")]
lean_rows = [r for r in listing if r[-1].endswith(".lean")]
names = set(mods)
indeg = collections.Counter()
for v in mods.values():
    for t in v[2]:
        if t in names:
            indeg[t] += 1

seen = {"Main"}
dq = collections.deque(["Main"])
while dq:
    x = dq.popleft()
    for t in mods[x][2]:
        if t in names and t not in seen:
            seen.add(t); dq.append(t)

deps = {m: {t for t in v[2] if t in names} for m, v in mods.items()}
users = collections.defaultdict(set)
for m, ds in deps.items():
    for d in ds:
        users[d].add(m)
left = {m: len(ds) for m, ds in deps.items()}
layer = {m: 0 for m, c in left.items() if c == 0}
q = collections.deque(layer)
done = 0
while q:
    x = q.popleft(); done += 1
    for u in users[x]:
        layer[u] = max(layer.get(u, 0), layer[x] + 1)
        left[u] -= 1
        if left[u] == 0:
            q.append(u)
hist = collections.defaultdict(lambda: [0, 0])
for m, d in layer.items():
    hist[d][0] += 1; hist[d][1] += mods[m][1]
top = collections.Counter()
for m, v in mods.items():
    top[m.split(".")[0]] += v[1]

R = {
    "sha256": h.hexdigest(),
    "files_listing": len(lean_rows), "bytes_listing": sum(int(r[2]) for r in lean_rows),
    "files": len(mods), "bytes": sum(v[0] for v in mods.values()), "lines": sum(v[1] for v in mods.values()),
    "imports": sum(len(v[2]) for v in mods.values()),
    "mathlib_all": sum(1 for v in mods.values() if "Mathlib" in v[2]),
    "unresolved_internal": sum(1 for v in mods.values() for t in v[2] if t not in names and not t.split(".")[0] in ("Mathlib", "Lean", "Batteries")),
    "reached_from_Main": len(seen), "acyclic": done == len(mods),
    "top_layer": layer["Main"], "layers": sorted([d, *hist[d]] for d in hist),
    "hubs": [[m, n, mods[m][1]] for m, n in indeg.most_common(12)],
    "top_lines": dict(top.most_common()),
    "theorems_and_lemmas": theorems,
    "median_lines": sorted(v[1] for v in mods.values())[len(mods) // 2],
}
json.dump(R, open("ref38.json", "w"), indent=0)
print({k: v for k, v in R.items() if k not in ("layers", "hubs")})
