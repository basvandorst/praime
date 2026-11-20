import json, os, sys, math
D = os.path.expanduser("~/praime/writeup/code/mersenne/strategies/ten/pairclocks/")
def load(x):
    fn = D + "count_%d.json" % x
    if not os.path.exists(fn):
        return None
    js = json.load(open(fn))
    return {d: [js["tally"][str(d)][k][-1] for k in ("N", "H2", "L2")] for d in js["d"]}
e9, e10 = load(10**9), load(10**10)
note = sys.argv[1] if len(sys.argv) > 1 else ""
big = {"note": note, "e9": e9, "e10": e10 or {}}
src = open("prime35.src.html").read()
open("prime35.html", "w").write(src.replace("BIGDATA", json.dumps(big, separators=(",", ":"))))
print("written; e10:", bool(e10))
