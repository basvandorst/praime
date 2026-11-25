"""timing proxies (excluded from the draws): Sym^3 at N ~ 215, T2 = 50; genus 2 with the largest N_L <= 3000, T2 = 150;
Rankin-Selberg L(E1 x E2) for 11a x 14a (built from local factors, conductor (N1 N2)^2), lfuncheckfeq for eps = +1 and -1, T2 = 60.
usage: c3_timing.py sym3|g2|rs"""
import json, time, sys
import cypari2
pari = cypari2.Pari(); pari.allocatemem(1100 * 10 ** 6)
what = sys.argv[1]
def run(L, T, label):
    t = time.time(); Li = pari.lfuninit(L, [1, T + 3]); ti = time.time() - t
    t = time.time(); z = pari.lfunzeros(Li, [0.01, T + 1]); tz = time.time() - t
    print('%s: init %.1fs, zeros to %d: %d in %.1fs' % (label, ti, T + 1, len(z), tz), flush=True)
if what == 'sym3':
    free = json.load(open('curves.json'))['free']
    v = min([u for u in free if u['w'] == -1], key=lambda u: abs(u['N'] - 215))
    E = pari.ellinit(v['ainv']); L = pari.lfunsympow(E, 3)
    json.dump(dict(sym3=v), open('timing_proxy_sym3.json', 'w'))
    run(L, 50, 'Sym3 N=%d %s T2=50' % (v['N'], v['ainv']))
elif what == 'g2':
    g = json.load(open('genus2.json'))['primitive'][-1]
    json.dump(dict(g2=g), open('timing_proxy_g2.json', 'w'))
    L = pari('lfungenus2([%s,%s])' % (g[1], g[2]))
    run(L, 150, 'genus 2 N=%d T2=150' % g[0])
else:
    pari('E1g = ellinit([0,-1,1,-10,-20]); E2g = ellinit([1,0,1,4,-6]); N1g = 11; N2g = 14')
    pari('rsloc(p) = my(a = ellap(E1g,p), b = ellap(E2g,p)); if(N1g%p==0, 1/(1 - a*b*X + a^2*p*X^2), if(N2g%p==0, 1/(1 - a*b*X + b^2*p*X^2), 1/(1 - a*b*X + (a^2*p + b^2*p - 2*p^2)*X^2 - a*b*p^2*X^3 + p^4*X^4)))')
    for eps in (1, -1):
        L = pari('lfuncreate([n -> direuler(p=2, n, rsloc(p)), 0, [-1,0,0,1], 3, (N1g*N2g)^2, %d])' % eps)
        print('RS 11a x 14a eps=%+d lfuncheckfeq %s' % (eps, pari.lfuncheckfeq(L)), flush=True)
    L = pari('lfuncreate([n -> direuler(p=2, n, rsloc(p)), 0, [-1,0,0,1], 3, (N1g*N2g)^2, 1])')
    run(L, 60, 'RS 11a x 14a T2=60')
