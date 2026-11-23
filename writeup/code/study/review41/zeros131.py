import cypari2, json, sys
pari = cypari2.Pari(); pari.allocatemem(400*10**6)
pari.default("realprecision", 30)
out = {}
for q in [5,7,11,13,17,19,23,29,31]:
    G = pari.znstar(q,1)
    g = int(pari('(G)->lift(G.gen[1])')(G))
    n = q-1
    res = {'g': g, 'chars': {}}
    for k in range(1, n):
        L = pari.lfuncreate([G, [k]])
        z = [float(t) for t in pari.lfunzeros(L, 50)]
        res['chars'][k] = z
    out[q] = res
    print(q, g, sum(len(v) for v in res['chars'].values()), flush=True)
json.dump(out, open('zeros131.json','w'))
