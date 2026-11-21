"""Aggregate graph.json + cone.json into the small data block for the page (data38.json)."""
import json, collections

G = json.load(open("graph.json"))
K = json.load(open("cone.json"))
depth = K["depth"]
names = set(G)

# tree to depth 3 (folder levels); leaves below are summed
tree = {}
for m, v in G.items():
    p = m.split(".")
    key = p[:3] if len(p) > 3 else p[:-1] if len(p) > 1 else p
    node = tree
    for part in key:
        node = node.setdefault(part, {})
    node.setdefault("__", [0, 0])
    node["__"][0] += 1
    node["__"][1] += v["l"]


def fold(name, node):
    files = lines = 0
    kids = []
    for k, sub in node.items():
        if k == "__":
            files += sub[0]; lines += sub[1]
            continue
        c = fold(k, sub)
        kids.append(c)
        files += c["f"]; lines += c["l"]
    kids.sort(key=lambda c: -c["l"])
    out = {"n": name, "f": files, "l": lines}
    if kids:
        own = node.get("__")
        if own:
            kids.append({"n": "(files here)", "f": own[0], "l": own[1]})
        out["c"] = kids
    return out


root = fold("src", tree)

# per depth layer: files, lines
layers = collections.defaultdict(lambda: [0, 0])
for m, d in depth.items():
    layers[d][0] += 1
    layers[d][1] += G[m]["l"]
L = [[d, *layers[d]] for d in sorted(layers)]

# in-degree (internal), hubs
indeg = collections.Counter()
for m, v in G.items():
    for t in v["i"]:
        if t in names:
            indeg[t] += 1
hubs = [[m, n, G[m]["l"]] for m, n in indeg.most_common(12)]

mathlib = collections.Counter()
for v in G.values():
    for t in v["i"]:
        if t.startswith("Mathlib."):
            mathlib[".".join(t.split(".")[:2])] += 1
mtop = mathlib.most_common(12)

kinds = collections.Counter()
for v in G.values():
    for k, n in v["k"].items():
        kinds[k] += n
sizes = sorted(v["l"] for v in G.values())
n = len(sizes)
stats = {
    "files": n, "lines": sum(sizes), "bytes": sum(v["b"] for v in G.values()),
    "imports": sum(len(v["i"]) for v in G.values()),
    "mathlib_all": sum(1 for v in G.values() if "Mathlib" in v["i"]),
    "internal_imports": sum(1 for v in G.values() for t in v["i"] if t in names),
    "mathlib_imports": sum(1 for v in G.values() for t in v["i"] if t.startswith("Mathlib")),
    "cone": len(K["cone"]), "chain": len(K["chain"]) - 1,
    "median_lines": sizes[n // 2], "max_lines": sizes[-1],
    "max_file": max(G, key=lambda m: G[m]["l"]),
    "kinds": dict(kinds.most_common()),
}
data = {"stats": stats, "tree": root, "layers": L, "hubs": hubs, "mathlib": mtop, "chain": K["chain"]}
json.dump(data, open("data38.json", "w"), separators=(",", ":"))
print(json.dumps(stats, indent=1))
print("tree top:", [(c["n"], c["l"], len(c.get("c", []))) for c in root["c"]])
print("layers:", len(L), "max layer files", max(L, key=lambda r: r[1]))
print("hubs:", hubs[:6])
print("mathlib:", mtop[:8])
import os
print("data38.json bytes", os.path.getsize("data38.json"))
