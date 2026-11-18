"""Arithmetic in Z[y], y^3 + y^2 - 2y - 1 = 0, for the three-cell wheel of the ring of 7 cells (u = 2 + y)."""
def mul(a, b, m=None):
    c = [0]*5
    for i in range(3):
        for j in range(3): c[i+j] += a[i]*b[j]
    for k in (4, 3):
        t = c[k]; c[k] = 0; c[k-1] -= t; c[k-2] += 2*t; c[k-3] += t
    return [x % m for x in c[:3]] if m else c[:3]
def norm(v):
    cols = [mul(v, e) for e in ([1, 0, 0], [0, 1, 0], [0, 0, 1])]
    a, b, c = cols
    return a[0]*(b[1]*c[2] - b[2]*c[1]) - b[0]*(a[1]*c[2] - a[2]*c[1]) + c[0]*(a[1]*b[2] - a[2]*b[1])
U = [2, 1, 0]
def upow(e, m):
    r = [1, 0, 0]; b = [x % m for x in U]
    while e:
        if e & 1: r = mul(r, b, m)
        b = mul(b, b, m); e >>= 1
    return r
def order(p, d):
    for k in sorted(k for k in range(1, d + 1) if d % k == 0):
        if upow(k, p) == [1, 0, 0]: return k
    return None
