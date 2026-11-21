"""Analyse graph.json: internal/external imports, cone of Main, depth, sizes per folder."""
import json, sys, collections

G = json.load(open("graph.json"))
names = set(G)
ext = collections.Counter()
missing = collections.Counter()
for m, v in G.items():
    for t in v["i"]:
        if t in names:
            continue
        root = t.split(".")[0]
        if root in ("Mathlib", "Lean", "Std", "Batteries", "Aesop", "Qq", "Init", "Plausible", "ProofWidgets", "ImportGraph", "LeanSearchClient", "Cli"):
            ext[root] += 1
        else:
            missing[t] += 1
print("external import roots:", dict(ext))
print("internal imports that name no file:", len(missing), list(missing.items())[:10])

# cone of Main (iterative DFS)
def cone(root):
    seen = {root}
    st = [root]
    while st:
        x = st.pop()
        for t in G[x]["i"]:
            if t in names and t not in seen:
                seen.add(t)
                st.append(t)
    return seen

C = cone("Main")
tot_l = sum(v["l"] for v in G.values())
c_l = sum(G[m]["l"] for m in C)
print(f"cone of Main: {len(C)} files of {len(G)}, {c_l} lines of {tot_l} ({100*c_l/tot_l:.1f} %)")

# topological order over the cone; detect cycles; longest import chain from Main
order = []
state = {}
for r in ["Main"]:
    st = [(r, 0)]
    while st:
        x, i = st.pop()
        if i == 0:
            if state.get(x) == 2:
                continue
            if state.get(x) == 1:
                continue
            state[x] = 1
        deps = [t for t in G[x]["i"] if t in names]
        if i < len(deps):
            st.append((x, i + 1))
            t = deps[i]
            if state.get(t) == 1:
                print("CYCLE at", x, "->", t)
            elif state.get(t) is None:
                st.append((t, 0))
        else:
            state[x] = 2
            order.append(x)
depth = {}
nxt = {}
for x in order:  # leaves first
    best, arg = 0, None
    for t in G[x]["i"]:
        if t in names and depth[t] + 1 > best:
            best, arg = depth[t] + 1, t
    depth[x] = best
    nxt[x] = arg
chain = ["Main"]
while nxt[chain[-1]]:
    chain.append(nxt[chain[-1]])
print("longest import chain from Main:", depth["Main"], "steps")
for x in chain[:12]:
    print("   ", x)
print("    ...", chain[-1])

# per top folder: files and lines, in cone / not in cone
agg = collections.defaultdict(lambda: [0, 0, 0, 0])
for m, v in G.items():
    k = m.split(".")[0]
    a = agg[k]
    a[0] += 1; a[1] += v["l"]
    if m in C:
        a[2] += 1; a[3] += v["l"]
print(f"{'folder':20s} {'files':>7s} {'lines':>10s} {'cone files':>11s} {'cone lines':>11s}")
for k, a in sorted(agg.items(), key=lambda kv: -kv[1][1]):
    print(f"{k:20s} {a[0]:7d} {a[1]:10d} {a[2]:11d} {a[3]:11d}")

kinds = collections.Counter()
ckinds = collections.Counter()
for m, v in G.items():
    for k, n in v["k"].items():
        kinds[k] += n
        if m in C:
            ckinds[k] += n
print("declarations all :", dict(kinds.most_common()))
print("declarations cone:", dict(ckinds.most_common()))
mathlib_files = sum(1 for m in G if any(t.startswith("Mathlib") for t in G[m]["i"]))
print("files importing Mathlib directly:", mathlib_files)
json.dump({"cone": sorted(C), "depth": depth, "chain": chain}, open("cone.json", "w"))
