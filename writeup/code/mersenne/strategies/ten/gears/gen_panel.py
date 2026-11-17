import json,re
d=json.load(open('analysis/machines.json'))
pc=lambda x:"%.1f%%"%(100*x)
rows=[('div','Wheels, one tooth per tick (machines 1 to 6, 8 to 10 and 18)','p itself'),
('square','Square jumps (15)','p itself'),
('double','Doubling wheels (11)','the order of 2: 3 (2), 7 (3), 5 (4), 17 (8), 31 (5), 127 (7)'),
('fib','Fibonacci wheels (12)','the first Fibonacci number p divides: 2 (3), 3 (4), 5 (5), 13 (7), 7 (8)'),
('accel','Accelerating wheels (14)','p − 1 (3 for p = 2)'),
('flip','Back and forth (16)','2, 5, 11, 17 for p = 2, 3, 5, 7; then 36 for every p'),
('pclock','Prime clock (13)','the p-th prime: 2 (3), 3 (5), 5 (11), 7 (17)'),
('slip','Slipping wheels (17)','p + 1'),
('fast','One dial (7)','60, 40, 24 for p = 2, 3, 5; then 120')]
def rec(r,N):
    x=d[r]['N'][str(N)]
    return "—" if r=='fast' else "%d (tick %s)"%(x['rec_p'][0],format(x['rec_p'][1],','))
tb="".join("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>"%(nm,pc(d[r]['N']['30']['precision']),pc(d[r]['N']['100']['precision']),pc(d[r]['N']['250']['precision']) if d[r]['N']['250']['quiet_share']>0 else "no quiet tick",pc(d[r]['N']['250']['recall']),rec(r,30),rec(r,250),fh) for r,nm,fh in rows)
a=d['accel']['N']['250']['accel_pairs']
html='''  <section class="panel" id="compare">
    <h2>All machines, run and compared</h2>
    <p class="note">Every rule of motion was run for 100,000 ticks with 30, 100 and 250 gears. A quiet tick is a tick at which no gear other than gear 1 is home. Of all ticks, %s are prime.</p>
    <ul class="found">
      <li><b>Four machines sieve the primes.</b> The plain wheels and the square jumps do it alike: with 250 gears %s of the quiet ticks are prime, and %s of the prime ticks are quiet (the rest are the primes up to 250, whose own gear is home). The doubling and the Fibonacci wheels do it at about half strength (%s and %s): there a product of two large primes is a quiet tick too.</li>
      <li><b>The accelerating wheels find pairs of primes.</b> With 250 gears there are %s quiet ticks below 250² = 62,500, and every one is either t = 2p with p and 2p + 1 prime, or a prime t with (t + 1)/2 prime.</li>
      <li><b>Three machines say nothing about the tick.</b> On the prime clock and the slipping wheels a quiet tick is prime as often as any tick (%s and %s). On the back-and-forth machine it is twice as often (%s), only because its distance is odd exactly at the odd ticks.</li>
      <li><b>More gears: harder or easier for primes to be home alone.</b> On the plain wheels the best within 100,000 ticks is 5, 4 and 3 prime gears with 30, 100 and 250 gears. On the slipping, Fibonacci and doubling wheels it grows: 8, 7 and 6 prime gears with 250.</li>
      <li><b>The doubling wheels as a sieve for Mersenne numbers.</b> A wheel can only divide 2<sup>p</sup> − 1 when p is at most half its size, so 250 wheels reach the 30 prime ticks up to 124. They rule out five (11, 23, 29, 37 and 83, by the wheels 23, 47, 233, 223 and 167), 11 are Mersenne exponents, and 14 pass although 2<sup>p</sup> − 1 is not prime. The Lucas–Lehmer wheels (19) decide them all.</li>
      <li><b>The Mersenne-exponent gears behave like the other prime gears.</b> On the doubling and the Fibonacci wheels they come home neither sooner nor later for their size (first tick divided by p, median %.2f against %.2f, and %.2f against %.2f; both within chance). On the other machines the first tick is a fixed rule (last column).</li>
    </ul>
    <div class="boardwrap"><table>
      <thead><tr><th>Machine</th><th>Quiet ticks that are prime, 30 gears</th><th>100 gears</th><th>250 gears</th><th>Prime ticks that are quiet, 250 gears</th><th>Most prime gears home alone, 30 gears</th><th>250 gears</th><th>First tick prime gear p is home</th></tr></thead>
      <tbody>%s</tbody>
    </table></div>
    <p class="note">On the dial the hands that are home are the multiples of one number, so prime hands are never home alone there (the record box counts rays instead); its quiet ticks are the ticks that share little with 120, a sieve by 2, 3 and 5 only. All figures were recomputed by a second program that moves every wheel tick by tick.</p>
  </section>
'''%(pc(d['div']['N']['250']['base']),pc(d['div']['N']['250']['precision']),pc(d['div']['N']['250']['recall']),pc(d['double']['N']['250']['precision']),pc(d['fib']['N']['250']['precision']),format(a[0],','),pc(d['pclock']['N']['250']['precision']),pc(d['slip']['N']['250']['precision']),pc(d['flip']['N']['250']['precision']),d['double']['exp_med'],d['double']['oth_med'],d['fib']['exp_med'],d['fib']['oth_med'],tb)
assert a[0]==a[1]
s=open('prime3.html',encoding='utf-8').read()
if 'id="compare"' in s:
    i=s.index('  <section class="panel" id="compare">'); j=s.index('</section>\n',i)+len('</section>\n'); s=s[:i]+s[j:]
k=s.index('</main>'); s=s[:k]+html+s[k:]
css='''h2 { margin:0 0 4px; font-size:17px; }
.found { margin:8px 0 12px; padding-left:20px; color:#d6e2f2; font-size:13px; }
.found li { margin:6px 0; }
.boardwrap { overflow:auto; border:1px solid var(--edge); border-radius:8px; }
table { border-collapse:collapse; width:100%; font-size:12.5px; font-variant-numeric:tabular-nums; }
th, td { padding:5px 8px; text-align:right; border-bottom:1px solid #1e2d44; vertical-align:top; }
th:first-child, td:first-child, th:last-child, td:last-child { text-align:left; min-width:150px; }
th { background:#16243a; color:#b8c8dc; font-weight:600; }
'''
if '.boardwrap' not in s: s=s.replace('a { color:var(--accent); }','a { color:var(--accent); }\n'+css)
if 'href="#compare"' not in s:
    x='At tick t the gears at home are exactly the divisors of t, and t is prime when only gear 1 and gear t are home. Press Start and watch.</p>'
    assert s.count(x)==1; s=s.replace(x,'At tick t the gears at home are exactly the divisors of t, and t is prime when only gear 1 and gear t are home. Press Start and watch. What the machines find when they all run for 100,000 ticks: <a href="#compare">the comparison</a> at the bottom.</p>')
open('prime3.html','w',encoding='utf-8').write(s)
print(len(s))
