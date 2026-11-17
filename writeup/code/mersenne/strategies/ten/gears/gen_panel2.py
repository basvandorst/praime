import json
r={x[0]:x for x in json.load(open('analysis/best.json'))['rows']}
fd=json.load(open('analysis/best2.json'))['finder']
f=lambda n:format(n,',')
rows=[('10. Back and forth',r['Back and forth (16)']),('8. Fibonacci wheels',r['Fibonacci (12)']),('9. Doubling wheels',r['Doubling (11)']),('7. Pattern machine, set to quadruplets',r['Pattern: quadruplets (20)']),('7. Pattern machine, set to twin primes',r['Pattern: twin primes (20)']),('5. Accelerating wheels',r['Accelerating (14)']),('6. Slipping wheels',r['Slipping (17)']),('1 to 3. Wheels, gear train, tape',r['Wheels, one tooth per tick (1-6, 8-10, 18, 21)']),('4. Square jumps',r['Square jumps (15)']),('11. Mersenne finder (250 prime wheels)',fd)]
tb="".join("<tr><td>%s</td><td>%d (tick %s)</td><td>%d (tick %s)</td><td>%d (tick %s)</td></tr>"%(nm,x[1],f(x[2]),x[3],f(x[4]),x[5],f(x[6])) for nm,x in rows)
html='''  <section class="panel" id="compare">
    <h2>Which machine has the best results</h2>
    <p class="note">Every machine was run for 100,000 ticks with 250 wheels. The table gives the most prime wheels, Mersenne-exponent wheels and twin-prime wheels that were home at once, and the tick. Among the numbers 1 to 250 there are 53 primes, 12 Mersenne exponents and 33 twin primes.</p>
    <div class="boardwrap"><table>
      <thead><tr><th>Machine</th><th>Primes at home</th><th>Mersenne exponents at home</th><th>Twin primes at home</th></tr></thead>
      <tbody>%s</tbody>
    </table></div>
    <ul class="found">
      <li><b>Back and forth wins, for a simple reason.</b> At tick 36 its wheels are exactly back where they started, so every wheel is home: all 53 primes, all 12 exponents, all 33 twin primes.</li>
      <li><b>Fibonacci and doubling are the real winners.</b> On these a wheel that is home at tick a is home at every multiple of a, so at a tick with many divisors (27,720, 55,440, 83,160) half of the prime wheels are home together.</li>
      <li><b>On the plain wheels at most six primes are home at once.</b> The tick must be a multiple of every wheel that is home, and 2 × 3 × 5 × 7 × 11 × 13 = 30,030 is the last such product below 100,000.</li>
      <li><b>For finding primes the order turns round.</b> The machines on which few wheels are home say the most: a tick with no wheel home is a prime (machines 1 to 4), a pair or a group of primes (7), or a candidate for a Mersenne prime (9 and 11).</li>
    </ul>
  </section>
'''%tb
s=open('prime3.html',encoding='utf-8').read()
i=s.index('  <section class="panel" id="compare">'); j=s.index('</section>\n',i)+len('</section>\n'); s=s[:i]+html+s[j:]
a='What the machines find when they all run for 100,000 ticks: <a href="#compare">the comparison</a> at the bottom.'
assert s.count(a)==1; s=s.replace(a,'Which machine has the best results: <a href="#compare">the table</a> at the bottom.')
a='th:first-child, td:first-child, th:last-child, td:last-child { text-align:left; min-width:150px; }'
assert s.count(a)==1; s=s.replace(a,'th:first-child, td:first-child { text-align:left; min-width:150px; }')
open('prime3.html','w',encoding='utf-8').write(s); print(len(s))
