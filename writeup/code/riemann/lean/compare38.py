import json, random, mpmath as mp
mp.mp.dps = 30
P = json.load(open("data38.json")); R = json.load(open("ref38.json")); Q = json.load(open("page38.json"))
out = []; bad = 0
def chk(name, a, b):
    global bad
    ok = a == b; bad += (not ok)
    out.append(f"{'OK ' if ok else 'BAD'} {name}: page {a} | reference {b}")
s = P["stats"]
for k in ["files", "bytes", "lines", "imports", "mathlib_all", "median_lines"]:
    chk(k, s[k], R[k])
chk("files in tar listing", s["files"], R["files_listing"])
chk("bytes in tar listing", s["bytes"], R["bytes_listing"])
chk("theorems+lemmas", s["kinds"]["theorem"] + s["kinds"]["lemma"], R["theorems_and_lemmas"])
chk("reached from Main", s["cone"], R["reached_from_Main"])
chk("top layer", s["chain"], R["top_layer"])
chk("chain length (files)", len(P["chain"]), R["top_layer"] + 1)
chk("layers", P["layers"], R["layers"])
chk("hubs", P["hubs"], R["hubs"])
chk("top folder lines", {c["n"]: c["l"] for c in P["tree"]["c"]}, R["top_lines"])
chk("acyclic, no unresolved", [True, 0], [R["acyclic"], R["unresolved_internal"]])
ref = []
k = 1
while True:
    z = mp.zetazero(k)
    if z.imag > 120: break
    ref.append(float(z.imag)); k += 1
chk("number of zeros up to t=120", len(Q["zeros"]), len(ref))
md = max(abs(a - b) for a, b in zip(Q["zeros"], ref))
ok = md < 1e-9; bad += not ok
out.append(f"{'OK ' if ok else 'BAD'} zeros up to 120: max |page - mpmath| = {md:.2e}")
me = 0
for (sr, si, zr, zi) in Q["zeta"]:
    w = mp.zeta(mp.mpc(sr, si)); me = max(me, abs(complex(w) - complex(zr, zi)) / max(1, abs(complex(w))))
ok = me < 1e-8; bad += not ok
out.append(f"{'OK ' if ok else 'BAD'} zeta at {len(Q['zeta'])} points in the drawn window: max relative error {me:.2e}")
mz = max(abs(float(mp.siegelz(t)) - v) for t, v in Q["Z"])
ok = mz < 0.01; bad += not ok
out.append(f"{'OK ' if ok else 'BAD'} Riemann-Siegel curve Z(t), {len(Q['Z'])} points 10..120: max |page - mpmath| = {mz:.2e} (drawing only)")
mn = min(v for _, _, v in Q["offline"])
out.append(f"info: |zeta(0.6 + i t)| at the zeros up to 120, smallest {mn:.4f}")
out.append(f"archive sha256 of part01..part14 concatenated: {R['sha256']}")
out.append(f"{bad} mismatches")
open("compare38.txt", "w").write("\n".join(out) + "\n"); print("\n".join(out))
