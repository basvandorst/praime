"""r77b: are the intercepts specific to x = W(y) = sum_{q>y} 1/q^2? Refit with other x(y) that also vanish as y -> oo."""
import json, math
import numpy as np
S = '/home/bas/praime/writeup/data/study/'
runs = {7: 'r77_y7.json', 11: 'r77_y11.json', 13: 'r77_y13.json', 17: 'r76.json'}
R = {y: json.load(open(S + f)) for y, f in runs.items()}
data = json.load(open(S + 'r73.json'))['diag']
r78 = {5: json.load(open(S + 'r78_y5.json')), 3: json.load(open(S + 'r78_y3.json'))}
pr = [q for q in range(2, 200000) if all(q % d for d in range(2, int(q ** .5) + 1))]
W = lambda y: sum(1 / q ** 2 for q in pr if q > y) + 1 / (200000 * math.log(200000))
forms = {'W=sum 1/q^2': W, '1/(y ln y)': lambda y: 1 / (y * math.log(y)), '1/y': lambda y: 1 / y,
         '1/y^2': lambda y: 1 / y ** 2, '1/ln y': lambda y: 1 / math.log(y)}
for name, f in forms.items():
    print(name)
    for q in (3, 5, 7, 11):
        ys = [y for y in runs if y >= q]
        x = np.array([f(y) for y in ys]); t = np.array([R[y]['mean_T_q%d' % q] for y in ys])
        B, A = np.polyfit(x, t, 1); res = np.max(np.abs(t - A - B * x))
        pred = {yp: A + B * f(yp) for yp in (5, 3) if yp >= q}
        obs = {yp: r78[yp]['mean_T_q%d' % q] for yp in pred}
        print('  q%-2d A=%.3f data=%.3f diff=%+.3f maxres=%.4f' % (q, A, data['mean_T_q%d' % q], A - data['mean_T_q%d' % q], res),
              ' '.join('y%d pred %.4f obs %.4f' % (yp, pred[yp], obs[yp]) for yp in pred))
