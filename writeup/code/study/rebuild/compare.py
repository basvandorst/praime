import json, math, os, sys
B = '/home/bas/scratch/rebuild-C/'
data = json.load(open(B + 'data.json'))
r81 = json.load(open('/home/bas/praime/writeup/data/study/r81.json'))
QQ = ['3', '5', '7', '11', '13']
def dstat(k, Dm):
    d = data[str(k)][str(Dm)]; n = d['n']
    return {q: (d['rep'][q] / n, math.sqrt(d['rep'][q] / n * (1 - d['rep'][q] / n) / n)) for q in QQ}
laws = {y: json.load(open(B + f'law_y{y}.json')) for y in (7, 11, 13) if os.path.exists(B + f'law_y{y}.json')}
for Dm in (200, 240):
    print(f'\n===== D <= {Dm} =====')
    for y, L in laws.items():
        for m in ('full', 'noR', 'lin2', 'S1', 'cramer'):
            zs = []
            for k in range(1, 5):
                ds = dstat(k, Dm)
                for q in QQ:
                    p, se = ds[q]; zs.append((L[m][str(Dm)]['rep'][str(k)][q] - p) / se)
            print(f'y={y:2d} {m:7s} chi2={sum(z*z for z in zs):10.1f}  max|z|={max(abs(z) for z in zs):6.2f}  '
                  f'k=1 max|z|={max(abs(z) for z in zs[:5]):6.2f}  k=2-4 max|z|={max(abs(z) for z in zs[5:]):6.2f}')
# detailed table for y=7 full at D<=200
for Dm in (200, 240):
    print(f'\nTable y=7 full, D<={Dm}: law / data / z   (+ noR z, S1 z)')
    L = laws[7]
    for k in range(1, 5):
        ds = dstat(k, Dm)
        row = []
        for q in QQ:
            p, se = ds[q]; t = L['full'][str(Dm)]['rep'][str(k)][q]
            zr = (L['noR'][str(Dm)]['rep'][str(k)][q] - p) / se
            zs1 = (L['S1'][str(Dm)]['rep'][str(k)][q] - p) / se
            row.append(f'q{q}: {t:.5f}/{p:.5f} z={(t-p)/se:+.2f} (noR {zr:+.1f}, S1 {zs1:+.1f})')
        print(f'k={k}', ' | '.join(row))
# against repo r81 (y=13, D<=240)
print('\nRepo r81 theory (y=13, D<=240) minus ours, in units of data se:')
for y, L in laws.items():
    mx = 0; mxa = 0
    for k in range(1, 5):
        ds = dstat(k, 240)
        for q in QQ:
            d = L['full']['240']['rep'][str(k)][q] - r81['theory'][str(k)][q]
            mx = max(mx, abs(d) / ds[q][1]); mxa = max(mxa, abs(d))
    print(f'y={y}: max |ours - r81| = {mxa:.2e} = {mx:.4f} sigma; mass ours {L["full"]["240"]["mass"]} vs r81 {list(r81["mass_within_DMAX"].values())}')
# data check
mx = max(abs(dstat(k, 240)[q][0] - r81['data'][str(k)][q]['p']) for k in range(1, 5) for q in QQ)
print('max |our data - r81 data| =', mx)
