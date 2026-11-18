"""The page's code for machine 19 (n2Test in prime3.html, run under node by n2_node.js) against norm2.py:
the start (a, c), both cells after n turns and the verdict, for every sampled list line, 64-bit primes, random odd numbers and all odd n to 20,001."""
import json, sys
from norm2 import start, power, home

rows = json.load(open(sys.argv[1]))
bad = 0
kinds = {'home': 0, 'not home': 0, 'no wheel': 0}
for s, a, c, x, y, ok, msg in rows:
    n = int(s)
    st = start(n)
    if st is None or st[0] == 'factor':
        good = x == '' and ok == 0
        kinds['no wheel'] += 1
    else:
        good = (a, c) == st and (int(x), int(y)) == power(a, c, n, n) and bool(ok) == bool(home(n))
        kinds['home' if ok else 'not home'] += 1
    bad += 0 if good else 1
t = "page code for machine 19 against norm2.py: %d numbers (%d came home, %d did not, %d had no wheel), different: %d" % (len(rows), kinds['home'], kinds['not home'], kinds['no wheel'], bad)
print(t)
open('norm2_js_check.txt', 'w').write(t + '\n')
assert bad == 0
