import json, re
src = open("prime38.src.html").read()
data = open("data38.json").read()
out = src.replace("__DATA__", data)
assert "__DATA__" not in out
body = re.sub(r"https?://[^\s\"'<>]+", "", out)
assert not re.search(r"(?<![:\\])//", body.split("<script>")[1]), "comment-like // in script"
assert "/*" not in out and "<!--" not in out, "comment found"
open("prime38.html", "w").write(out)
print("prime38.html", len(out), "bytes")
