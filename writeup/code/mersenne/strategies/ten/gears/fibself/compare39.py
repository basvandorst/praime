import json, sys
P = json.load(open(sys.argv[1])); R = json.load(open("ref39.json"))
bad = 0; out = []
def chk(name, ok, detail=""):
    global bad; bad += (not ok); out.append(f"{'OK ' if ok else 'BAD'} {name} {detail}")
g = {int(k): v for k, v in P["grid"].items()}
diff = [k for k in range(3, 1501) if g[k] != R["grid"][str(k)]]
chk("grid k = 3..1500, page (index rule) vs gmpy2 on F(k)", not diff, f"differences {diff[:10]}")
chk("Fibonacci pseudoprimes below 1e7", P["fpp"] == R["fpp1e7"], f"page {len(P['fpp'])} ref {len(R['fpp1e7'])}")
fer = sum(1 for n in R["fpp1e7"] if pow(2, n - 1, n) == 1)
chk("those that also fool the doubling wheel", P["fer"] == fer, f"page {P['fer']} python {fer}")
for name, v in R["examples"].items():
    pv = P["ex"][name]; chk(f"example {name} ({v['digits']} digits)", pv["home"] == v["home"] and pv["prime"] == v["prime"], f"home {pv['home']} prime {pv['prime']}")
out.append(f"{bad} mismatches"); open("compare39.txt", "w").write("\n".join(out) + "\n"); print("\n".join(out))
