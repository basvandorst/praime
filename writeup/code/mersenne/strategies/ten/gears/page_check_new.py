"""Runs machines 20 to 23 of the page in headless Chrome with 5,000 wheels and writes tour/page_check_20_23.txt."""
import subprocess, json, re, os
OUT = '$("status").innerText.match(/Lamps lit so far: [\\d,]+\\. Mistakes so far: \\d+\\./)[0]'
spec = [('20', 'tick 11', 'v=sq&t=11&n=5000'), ('20', 'tick 12', 'v=sq&t=12&n=5000'), ('21', 'tick 18', 'v=sf&t=18&n=5000'), ('21', 'tick 19', 'v=sf&t=19&n=5000'),
        ('22', 'tick 5,000', 'v=cnt&t=5000&n=5000'), ('23', 'round 2', 'v=inh&t=2&n=5000'), ('23', 'round 3', 'v=inh&t=3&n=5000')]
runs = [dict(q=q, tag='pc%d' % i, wait=2000, out=OUT) for i, (m, t, q) in enumerate(spec)]
o = subprocess.run(['node', 't3.js', os.getcwd() + '/prime3.html', json.dumps(runs)], capture_output=True, text=True, timeout=400).stdout
outs = re.findall(r'  out: (.*)', o)
errs = re.findall(r'^pc\d+ .*? errors (\d+)', o, re.M)
assert len(outs) == len(spec) == len(errs), o[:500]
lines = ["The page's own code in headless Chrome, 5,000 wheels (status line of the machine at the tick named):"]
for (m, t, q), st, e in zip(spec, outs, errs):
    lines.append('machine %s at %s: %s script errors %s' % (m, t, st, e))
open('tour/page_check_20_23.txt', 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
