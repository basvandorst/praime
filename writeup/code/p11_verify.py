from sympy import isprime, nextprime

def run(p, n):
    ps = [p]
    while len(ps) < n + 1:
        ps.append(nextprime(ps[-1]))
    assert all(isprime(q) for q in ps)
    return [b - a for a, b in zip(ps, ps[1:])]

distinct = {21: 69931991, 22: 203674907, 23: 1092101119, 24: 1363592621, 25: 1363592677,
            26: 2124140323, 27: 23024158649, 28: 30282104173, 29: 30282104173,
            30: 196948778371, 31: 263552821783, 32: 691340780027,
            33: 1625800359439, 34: 1625800359439, 35: 1625800359439, 36: 1625800359439}
for n, p in distinct.items():
    g = run(p, n)
    assert len(set(g)) == n, n                      # n pairwise distinct gaps

# The page's increasing / decreasing records, confirmed first below 1.7e12 by our pass.
g = run(1472840004017, 14); assert all(b > a for a, b in zip(g, g[1:]))
g = run(1273135176799, 13); assert all(b < a for a, b in zip(g, g[1:]))
print("all runs verified")
