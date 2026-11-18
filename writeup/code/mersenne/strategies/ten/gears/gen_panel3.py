import json
r={x[0]:x for x in json.load(open('analysis/best.json'))['rows']}
fd=json.load(open('analysis/best2.json'))['finder']
T=json.load(open('tour/tournament.json'))
f=lambda n:format(n,',')
pas=['pascal']+T['row250']
duo=['duo']+T['duo_row250']
rows=[("10. Pascal's wheel",pas),("11. Two-cell wheel",duo),('8. Back and forth',r['Back and forth (16)']),('6. Fibonacci wheels',r['Fibonacci (12)']),('7. Doubling wheels',r['Doubling (11)']),('5. Pattern machine, set to quadruplets',r['Pattern: quadruplets (20)']),('5. Pattern machine, set to twin primes',r['Pattern: twin primes (20)']),('3. Accelerating wheels',r['Accelerating (14)']),('4. Slipping wheels',r['Slipping (17)']),('1 and 2. Wheels, gear train',r['Wheels, one tooth per tick (1-6, 8-10, 18, 21)']),('9. Mersenne finder (250 prime wheels)',fd)]
tb="".join("<tr><td>%s</td><td>%d (tick %s)</td><td>%d (tick %s)</td><td>%d (tick %s)</td></tr>"%(nm,x[1],f(x[2]),x[3],f(x[4]),x[5],f(x[6])) for nm,x in rows)
S=T['single']; B='2<sup>64</sup>'
def one(name,cells,first,found):
    return "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>"%(name,cells,first,found)
cen=T['census']; rings=sorted(int(k) for k in cen)
t2=one('Doubling wheel: one cell, doubled every turn',1,f(S['doubling']['first']),f(S['doubling']['below1e8'])+' below 10<sup>8</sup>; '+f(T['list_count'])+' below '+B)
t2+=one('Fibonacci wheel: two cells, each turn the sum of the two',2,f(S['fib']['first']),f(S['fib']['below1e8'])+' below 10<sup>8</sup>')
t2+=one("Perrin's wheel: three cells with a feedback",3,f(S['perrin']['first']),f(S['perrin']['below1e8'])+' below 10<sup>8</sup>')
t2+=one("Pascal's wheel, always 5 cells",5,f(S['ring5']['first']),f(S['ring5']['below1e8'])+' below 10<sup>8</sup>; '+f(S['ring5']['below1e10'])+' below 10<sup>10</sup>')
t2+=one('Doubling wheel and Fibonacci wheel together, for numbers that leave 2 or 3 after division by 5',3,'none','0 below '+B+' ('+f(T['pair']['tested'])+' candidates)')
t2+=one("Pascal's wheel on a ring the number does not fit (machine 10)",'%d to %d'%(T['misfit']['min_ring'],T['misfit']['max_ring']),'none','0 below '+B+' ('+f(T['misfit']['tested'])+' candidates)')
t2+=one("Two-cell wheel with a step the number does not fit (machine 11)",2,'none','0 among all odd numbers below '+S['duo']['bound_html']+'; with the doubling wheel 0 below '+B+' ('+f(T['duo_list']['tested'])+' candidates)')
ringtxt=(" We also tried every one of them on each ring of "+", ".join(str(x) for x in rings[:-1])+" and "+str(rings[-1])+" cells that it does not fit: none came home.") if len(rings)>1 and all(cen[str(r)]['liars']==0 for r in rings) else ""
html='''  <section class="panel" id="compare">
    <h2>Which machine has the best results</h2>
    <p class="note">Every machine was run for 100,000 ticks with 250 wheels. The table gives the most prime wheels, Mersenne-exponent wheels and twin-prime wheels that were home at once, and the tick. Among the numbers 1 to 250 there are 53 primes, 12 Mersenne exponents and 33 twin primes.</p>
    <div class="boardwrap"><table>
      <thead><tr><th>Machine</th><th>Primes at home</th><th>Mersenne exponents at home</th><th>Twin primes at home</th></tr></thead>
      <tbody>%s</tbody>
    </table></div>
    <ul class="found">
      <li><b>Pascal's wheel and the two-cell wheel win without a trick.</b> Each of their lamps is a decision: wheel n goes home at tick n only when n is prime, and stays there. At tick 241 all 53 prime lamps are on, and no other lamp.</li>
      <li><b>Back and forth also has them all, for a simple reason.</b> At tick 36 its wheels are exactly back where they started, so every wheel is home, prime or not.</li>
      <li><b>Fibonacci and doubling come next.</b> On these a wheel that is home at tick a is home at every multiple of a, so at a tick with many divisors (27,720, 55,440, 83,160) half of the prime wheels are home together.</li>
      <li><b>On the plain wheels at most six primes are home at once.</b> The tick must be a multiple of every wheel that is home, and 2 × 3 × 5 × 7 × 11 × 13 = 30,030 is the last such product below 100,000.</li>
      <li><b>Two machines were removed.</b> The tape and the square jumps gave the same results as machine 1.</li>
    </ul>
    <h2 id="decides">One wheel decides: which wheel can be trusted</h2>
    <p class="note">A harder contest. Wheel n must be home after n turns exactly when n is prime. On every machine below each prime comes home; that is proved. The mistakes are the numbers that are not prime and come home too. All machines were run against the odd numbers that are not prime.</p>
    <div class="boardwrap"><table class="dec">
      <thead><tr><th>Machine</th><th>Cells</th><th>First mistake</th><th>Mistakes found</th></tr></thead>
      <tbody>%s</tbody>
    </table></div>
    <ul class="found">
      <li><b>The rule that wins: give the number a wheel it does not fit.</b> A number n fits a ring of r cells when n × n leaves remainder 1 after division by r, and it fits a two-cell wheel with step c when c is a square in the arithmetic of remainders after division by n. Every mistake of Pascal's wheel that we found is a number that fits its ring: %s on rings of 2 to 40 cells below 10<sup>10</sup>, not one that does not fit.</li>
      <li><b>Checked as far as the lists go.</b> A number that fools Pascal's wheel must fool the doubling wheel first (add up the cells: the sum doubles every turn). All %s odd numbers below %s that fool the doubling wheel are known, and we put every one on the smallest ring it does not fit and on its two-cell wheel: none came home.%s</li>
      <li><b>Is this the gear that amazes mathematicians?</b> Not yet, and we say so plainly. The last three rows are open problems. Nobody has proved that they never make a mistake. For the pair of wheels there is a prize of 620 dollars, offered by Pomerance, Selfridge and Wagstaff; for Pascal's wheel the claim is Agrawal's conjecture; the two-cell wheel is a Frobenius test, for which no mistake is known either. There are reasons to expect a mistake at some very large number (Pomerance for the pair; Lenstra and Pomerance for Pascal's wheel). A proof for any of the three, or the first number that fools one, would be news.</li>
      <li><b>What they can do today.</b> Both wheels need a few cells and about as many steps as the number has binary digits, so they also run on numbers of a thousand digits: choose machine 10 or 11 and type one.</li>
    </ul>
  </section>
'''%(tb,t2,f(T['fit_liars']),f(T['list_count']),B,ringtxt)
s=open('prime3.html',encoding='utf-8').read()
i=s.index('  <section class="panel" id="compare">'); j=s.index('</section>\n',i)+len('</section>\n'); s=s[:i]+html+s[j:]
import re
s,k=re.subn(r'In our runs none of the [\d,]+ odd numbers below [\d,]+ that are not prime came home','In our runs none of the %s odd numbers below %s that are not prime came home'%(f(S['duo']['tested']),f(S['duo']['bound'])),s); assert k==1
open('prime3.html','w',encoding='utf-8').write(s); print(len(s))
