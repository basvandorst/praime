"""One streaming pass over the Lean archive parts.

For every .lean file: module name, bytes, lines, number of theorem/lemma/def-like
declarations, and the list of imported modules. Writes graph.json.
Usage: python3 -I graph_scan.py PARTS_DIR OUT.json
"""
import sys, re, tarfile, glob, json


class Cat:
    def __init__(self, files):
        self.files = files
        self.i = 0
        self.f = open(files[0], "rb")

    def read(self, n=-1):
        out = b""
        while n < 0 or len(out) < n:
            chunk = self.f.read(-1 if n < 0 else n - len(out))
            if chunk:
                out += chunk
                continue
            self.i += 1
            if self.i >= len(self.files):
                break
            self.f.close()
            self.f = open(self.files[self.i], "rb")
        return out


IMPORT = re.compile(r"^\s*(?:(?:public|private|meta)\s+)*import\s+(?:all\s+)?([^\n]+)", re.M)
DECL = re.compile(r"^\s*(?:@\[[^\]]*\]\s*)?(?:(?:private|protected|noncomputable|nonrec|partial)\s+)*"
                  r"(theorem|lemma|def|abbrev|instance|structure|class|inductive)\b", re.M)


def header(txt):
    """Text before the first non-import, non-comment, non-blank line."""
    code = re.sub(r"/-.*?-/", "", txt, flags=re.S)
    code = re.sub(r"--[^\n]*", "", code)
    lines = []
    for ln in code.split("\n"):
        s = ln.strip()
        if not s or s == "module":
            continue
        if re.match(r"(?:(?:public|private|meta)\s+)*import\s", s):
            lines.append(s)
            continue
        break
    return "\n".join(lines), code


def main():
    parts = sorted(glob.glob(sys.argv[1] + "/lean.tar.gz.part*"))
    mods = {}
    with tarfile.open(fileobj=Cat(parts), mode="r|gz") as tf:
        for m in tf:
            if not m.isfile() or not m.name.endswith(".lean"):
                continue
            assert m.name.startswith("lean/src/"), m.name
            mod = m.name[len("lean/src/"):-len(".lean")].replace("/", ".")
            raw = tf.extractfile(m).read()
            txt = raw.decode("utf-8", "replace")
            head, code = header(txt)
            imps = []
            for mm in IMPORT.finditer(head):
                imps += mm.group(1).split()
            kinds = {}
            for mm in DECL.finditer(code):
                kinds[mm.group(1)] = kinds.get(mm.group(1), 0) + 1
            mods[mod] = {"b": len(raw), "l": txt.count("\n") + (0 if txt.endswith("\n") or not txt else 1),
                         "i": imps, "k": kinds}
    json.dump(mods, open(sys.argv[2], "w"), separators=(",", ":"))
    print("modules", len(mods), "bytes", sum(v["b"] for v in mods.values()),
          "lines", sum(v["l"] for v in mods.values()),
          "imports", sum(len(v["i"]) for v in mods.values()))


main()
