import json, time, sys
import cypari2
pari = cypari2.Pari(); pari.allocatemem(1100 * 10 ** 6)
free = json.load(open('curves.json'))['free']
neg = [v for v in free if v['w'] == -1]
prox = [[v for v in neg if v['N'] >= 180][0], [v for v in neg if v['N'] >= 140][0]]
for v in prox:
    E = pari.ellinit(v['ainv']); L = pari.lfunsympow(E, 3)
    try:
        t = time.time(); Li = pari.lfuninit(L, [1, 73]); ti = time.time() - t
        t = time.time(); z = pari.lfunzeros(Li, [0.01, 71]); tz = time.time() - t
        print('Sym3 proxy N=%d %s w(E)=%d eps(Sym3)=%s: init to 73 %.1fs, zeros to 71: %d in %.1fs' % (v['N'], v['ainv'], v['w'], pari.lfunrootres(L)[2], ti, len(z), tz), flush=True)
    except Exception as e:
        print('N=%d failed: %s' % (v['N'], str(e)[:100]), flush=True)
json.dump(prox, open('proxies_sym3.json', 'w'))
