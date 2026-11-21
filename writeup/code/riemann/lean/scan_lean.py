"""Stream the lean tarball parts and scan every .lean file for trust-relevant constructs."""
import sys, re, tarfile, glob, collections, json

parts = sorted(glob.glob(sys.argv[1] + "/lean.tar.gz.part*"))

class Cat:
    def __init__(self, files):
        self.files = files; self.i = 0; self.f = open(files[0], "rb")
    def read(self, n=-1):
        out = b""
        while n < 0 or len(out) < n:
            chunk = self.f.read(-1 if n < 0 else n - len(out))
            if chunk:
                out += chunk; continue
            self.i += 1
            if self.i >= len(self.files):
                break
            self.f.close(); self.f = open(self.files[self.i], "rb")
        return out

PATS = {
    "axiom": r"^\s*(private\s+|protected\s+)?axiom\s",
    "sorry": r"\bsorry\b",
    "admit": r"\badmit\b",
    "native_decide": r"\bnative_decide\b|Lean\.ofReduceBool|ofReduceNat",
    "implemented_by/extern": r"implemented_by|@\[extern",
    "unsafe": r"\bunsafe\b",
    "opaque": r"^\s*(private\s+)?opaque\s",
    "run_cmd/elab-env": r"\brun_cmd\b|\brun_elab\b|\brun_meta\b|modifyEnv|addDecl|addAndCompile|Environment\.add|setEnv",
    "kernel-off": r"debug\.skipKernelTC|skipKernelTC|Kernel\.Environment|unsafeCast",
    "macro/syntax": r"^\s*(macro|syntax|elab|macro_rules|elab_rules)\b",
    "def RiemannHypothesis": r"\bRiemannHypothesis\b",
    "redef zeta/LFunction": r"(def|abbrev|notation|theorem)\s+(riemannZeta|DirichletCharacter\.LFunction)\b",
    "set_option": r"^\s*set_option\s",
}
rx = {k: re.compile(v, re.M) for k, v in PATS.items()}
hits = collections.defaultdict(list)
counts = collections.Counter()
nfiles = 0; nlines = 0
with tarfile.open(fileobj=Cat(parts), mode="r|gz") as tf:
    for m in tf:
        if not m.isfile() or not m.name.endswith(".lean"):
            continue
        txt = tf.extractfile(m).read().decode("utf-8", "replace")
        nfiles += 1; nlines += txt.count("\n")
        # strip comments for keyword checks (block and line comments, rough)
        code = re.sub(r"/-.*?-/", "", txt, flags=re.S)
        code = re.sub(r"--[^\n]*", "", code)
        for k, r in rx.items():
            found = r.findall(code)
            if found:
                counts[k] += len(found)
                if len(hits[k]) < 400:
                    for mm in r.finditer(code):
                        line = code[code.rfind("\n", 0, mm.start()) + 1: code.find("\n", mm.end())]
                        hits[k].append((m.name, line.strip()[:200]))
                        if len(hits[k]) >= 400:
                            break
print("files", nfiles, "lines", nlines)
for k in PATS:
    print(f"{k:28s} {counts[k]}")
json.dump(hits, open(sys.argv[2], "w"), indent=0)
