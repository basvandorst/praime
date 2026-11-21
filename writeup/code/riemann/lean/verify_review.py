import sys, re, tarfile, glob, collections
sys.path.insert(0, sys.argv[2])
from graph_scan_cat import Cat
parts = sorted(glob.glob(sys.argv[1] + "/lean.tar.gz.part*"))
P = {"kexactB": r"\bkexactB\b", "kc%": r"\bkc%", "kcopy": r"\bkcopy\b", "run_tac": r"\brun_tac\b",
     "hygiene false": r"set_option\s+hygiene\s+false", "allowUnsafeReducibility": r"allowUnsafeReducibility",
     "respectTransparency false": r"backward\.isDefEq\.respectTransparency\s+false"}
rx = {k: re.compile(v) for k, v in P.items()}
uses = collections.Counter(); files = collections.defaultdict(set)
DEFP = re.compile(r"^\s*(?:@\[[^\]]*\]\s*)?(?:(?:private|protected|noncomputable)\s+)*(?:def|abbrev)\s+\S+(.*?):\s*Prop\s*:=", re.M)
lit_noparam = lit_param = 0; cex_files = cex_lines = 0
EXACT = {"PublishedTheorems", "LiteratureResults", "PublishedInputs"}
with tarfile.open(fileobj=Cat(parts), mode="r|gz") as tf:
    for m in tf:
        if not m.isfile() or not m.name.endswith(".lean"): continue
        txt = tf.extractfile(m).read().decode("utf-8", "replace")
        code = re.sub(r"/-.*?-/", "", txt, flags=re.S); code = re.sub(r"--[^\n]*", "", code)
        for k, r in rx.items():
            n = len(r.findall(code))
            if n: uses[k] += n; files[k].add(m.name)
        dirs = m.name.split("/")[:-1]
        if EXACT & set(dirs):
            for mm in DEFP.finditer(code):
                if mm.group(1).strip(): lit_param += 1
                else: lit_noparam += 1
        if any("Counterexample" in d for d in dirs):
            cex_files += 1; cex_lines += txt.count("\n")
for k in P: print(f"{k:28s} uses {uses[k]:6d} files {len(files[k]):5d}")
kfiles = set().union(*(files[k] for k in ["kexactB", "kc%", "kcopy", "run_tac"]))
print("files with kexactB/kc%/kcopy/run_tac:", len(kfiles), "uses:", sum(uses[k] for k in ["kexactB", "kc%", "kcopy", "run_tac"]))
print("Prop defs in exact literature folders: no params", lit_noparam, "with params", lit_param)
print("Counterexample folders: files", cex_files, "lines", cex_lines)
