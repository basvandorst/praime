"""Independent reference for prime8.html: the turtle's end position, turns, left and right turns, farthest distance and
crossings of the start, for every deterministic rule, numbers up to X.   usage: ref8.py X out.json"""
import sys, json, math
X = int(sys.argv[1])
MEXP = {2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593}
L = X + 3
om = [0] * (L + 1); rem = list(range(L + 1))
for p in range(2, L + 1):
    if rem[p] == p and om[p] == 0:
        for m in range(p, L + 1, p):
            while rem[m] % p == 0: rem[m] //= p; om[m] += 1
DX = [1, 0, -1, 0]; DY = [0, 1, 0, -1]
out = {}
for rule in ('left', 'mod4', 'dig', 'twin', 'liou', 'mexp'):
    x = y = d = turns = Lt = Rt = ret = 0; far = 0.0
    for n in range(1, X + 1):
        ox = x; x += DX[d]; y += DY[d]
        if x == 0 and y == 0: ret += 1
        elif ox != 0 and y == 0 and (ox > 0) != (x > 0) and x != 0: ret += 1
        pr = om[n] == 1
        if rule == 'left': t = 1 if pr else 0
        elif rule == 'mod4': t = (1 if n % 4 == 1 else -1) if pr and n > 2 else 0
        elif rule == 'dig': t = (1 if n % 10 in (1, 9) else -1) if pr and n > 5 else 0
        elif rule == 'twin': t = (1 if om[n + 2] == 1 else -1) if pr and n > 2 else 0
        elif rule == 'liou': t = 1 if om[n] % 2 == 0 else -1
        else: t = 1 if n in MEXP else 0
        if t:
            d = (d + (1 if t > 0 else 3)) % 4; turns += 1
            if t > 0: Lt += 1
            else: Rt += 1
        far = max(far, math.hypot(x, y))
    out[rule] = dict(x=x, y=y, turns=turns, L=Lt, R=Rt, ret=ret, far=round(far, 6))
json.dump(out, open(sys.argv[2], 'w')); print(out)
