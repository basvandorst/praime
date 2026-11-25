"""seeded selection of the degree-4 test set (SEED4 = 20261126; numpy default_rng(SEED4 + k) permutations of sorted lists).
Sym^3 E (free semistable classes of curves.json, keyed on minimal models; the two timing proxies (N = 142 and 215, minimal models) excluded):
  S70: N <= 150, T2 = 70: every class with w(E) = -1 (11) and the same number with w(E) = +1 (draw k = 1);
  S50: 150 < N <= 230, T2 = 50: n50 = min(#w=-1, 19) classes with w = -1 (draw k = 2) and n50 with w = +1 (draw k = 3).
  eps(Sym^3 E) = w(E) for semistable E (local root numbers: w_inf(Sym^3) = w_inf(E) = -1 and w_p(Sym^3) = w_p(E) = -a_p at
  multiplicative p); checked against PARI (lfunsympow) for every chosen curve, and the opposite sign must fail lfuncheckfeq.
Genus 2: the primitive list of genus2.json without the timing proxy (largest N, 2981), draw k = 4: 40 test, 2 proxies.
Rankin-Selberg E1 x E2: all semistable classes with N <= 60 from the curve box (used curves allowed: the RS L-function is new),
  pairs with gcd(N1, N2) = 1 and N1 N2 <= 1500, not the timing pair 11a x 14a; draw k = 5: 20 test, 1 proxy; T2 = 60 (set after timing)."""
import json, math
import numpy as np, cypari2
pari = cypari2.Pari(); pari.allocatemem(600 * 10 ** 6)
SEED4 = 20261126
free = json.load(open('curves.json'))['free']
prox = [[int(x) for x in pari.ellminimalmodel(pari.ellinit(a))[:5]] for a in ([1, -1, 1, -12, 15], [0, 0, 1, -8, -12])]
free = [v for v in free if v['ainv'] not in prox]
def draw(lst, k, n):
    rng = np.random.default_rng(SEED4 + k); perm = rng.permutation(len(lst)); return [lst[i] for i in perm[:n]], [lst[i] for i in perm[n:n + 2]]
s70m = [v for v in free if v['N'] <= 150 and v['w'] == -1]
s70p_all = [v for v in free if v['N'] <= 150 and v['w'] == 1]
s70p, s70p_prox = draw(s70p_all, 1, len(s70m))
s50m_all = [v for v in free if 150 < v['N'] <= 230 and v['w'] == -1]
s50p_all = [v for v in free if 150 < v['N'] <= 230 and v['w'] == 1]
n50 = min(len(s50m_all), 19)
s50m, s50m_prox = draw(s50m_all, 2, n50)
s50p, s50p_prox = draw(s50p_all, 3, n50)
sym3 = [dict(v, stratum='S70', T2=70) for v in s70m + s70p] + [dict(v, stratum='S50', T2=50) for v in s50m + s50p]
chk = []
for v in sym3:
    E = pari.ellinit(v['ainv']); L = pari.lfunsympow(E, 3)
    eps = int(round(float(pari.real(pari.lfunrootres(L)[2]))))
    v['eps_sym3'] = eps; chk.append(eps == v['w'])
print('Sym3: S70 %d (w=-1 %d, w=+1 %d), S50 %d (w=-1 %d of %d available, w=+1 %d of %d); eps(Sym3) = w(E) for all: %s' % (
    len(s70m + s70p), len(s70m), len(s70p), len(s50m + s50p), len(s50m), len(s50m_all), len(s50p), len(s50p_all), all(chk)))
g2 = json.load(open('genus2.json'))['primitive'][:-1]
g2t, g2p = draw(g2, 4, 40)
# Rankin-Selberg pairs
PR = [int(p) for p in pari.primes([2, 97])]
classes = {}
for a1 in (0, 1):
    for a2 in (-1, 0, 1):
        for a3 in (0, 1):
            for a4 in range(-30, 31):
                for a6 in range(-60, 61):
                    try:
                        E = pari.ellinit([a1, a2, a3, a4, a6])
                        if len(E) == 0: continue
                        N = int(pari.ellglobalred(E)[0])
                    except Exception:
                        continue
                    if N > 60 or not pari.issquarefree(N): continue
                    E = pari.ellminimalmodel(E)
                    key = (N,) + tuple(int(pari.ellap(E, p)) for p in PR)
                    if key not in classes: classes[key] = [N, [int(x) for x in E[:5]]]
cl = sorted(classes.values())
pairs = [[x, y] for i, x in enumerate(cl) for y in cl[i + 1:] if math.gcd(x[0], y[0]) == 1 and x[0] * y[0] <= 1500
         and not (x[1] == [0, -1, 1, -10, -20] and y[1] == [1, 0, 1, 4, -6])]
rng = np.random.default_rng(SEED4 + 5); perm = rng.permutation(len(pairs))
rst = [pairs[i] for i in perm[:20]]; rsp = [pairs[i] for i in perm[20:21]]
print('genus 2: %d test of %d (N %d-%d), proxies %s' % (len(g2t), len(g2), min(r[0] for r in g2t), max(r[0] for r in g2t), [r[0] for r in g2p]))
print('RS: classes N <= 60: %d, coprime pairs with N1 N2 <= 1500: %d; test N1 x N2: %s' % (len(cl), len(pairs), [(p[0][0], p[1][0]) for p in rst]))
json.dump(dict(seed=SEED4, rule=__doc__, sym3=sym3, sym3_proxies=s70p_prox + s50m_prox + s50p_prox, g2=g2t, g2_proxies=g2p,
               rs=rst, rs_proxies=rsp, timing_proxies=dict(sym3=prox, g2=json.load(open('genus2.json'))['primitive'][-1], rs=[[11, [0, -1, 1, -10, -20]], [14, [1, 0, 1, 4, -6]]])),
          open('testset4.json', 'w'), indent=1)
