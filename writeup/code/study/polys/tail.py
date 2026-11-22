"""Sum the run histograms in all6.err and fit the tail decay per extra length (runs >= 36). usage: tail.py"""
import re, math, collections
H = collections.Counter(); tested = 0; shapes = 0
for line in open('all6.err'):
    m = re.search(r'shapes (\d+) .* tested (\d+)', line)
    if not m: continue
    shapes += int(m.group(1)); tested += int(m.group(2))
    for a, b in re.findall(r'(\d+):(\d+)', line.split('|')[1]): H[int(a)] += int(b)
cum = {L: sum(v for k, v in H.items() if k >= L) for L in range(30, 60)}
print('shapes', shapes, 'tested', tested)
print('runs >= L:', {L: cum[L] for L in range(34, 50) if cum[L]})
pts = [(L, math.log(cum[L])) for L in range(34, 46) if cum[L] >= 3]
n = len(pts); mx = sum(p[0] for p in pts) / n; my = sum(p[1] for p in pts) / n
sl = sum((p[0] - mx) * (p[1] - my) for p in pts) / sum((p[0] - mx) ** 2 for p in pts)
print('decay per length %.3f; expected runs >= 56 so far %.3g; per million tested %.3g' % (math.exp(sl), math.exp(my + sl * (56 - mx)), math.exp(my + sl * (56 - mx)) / tested * 1e6))
