import json
r={x[0]:x for x in json.load(open('analysis/best.json'))['rows']}
fd=json.load(open('analysis/best2.json'))['finder']
T=json.load(open('tour/tournament.json'))
f=lambda n:format(n,',')
pas=['pascal']+T['row250']
duo=['duo']+T['duo_row250']
R12=json.load(open('tour/rows_12_13.json'))
rows=[("10. Pascal's wheel",pas),("11. Two-cell wheel",duo),("12. Three-cell wheel",['tri']+R12['three']),("13. Wilson's wheel",['wil']+R12['wilson']),("14. Lucas's wheel",['luc']+R12['lucas']),("15. Two cells with check turns",['np']+R12['nplus']),('8. Back and forth',r['Back and forth (16)']),('6. Fibonacci wheels',r['Fibonacci (12)']),('7. Doubling wheels',r['Doubling (11)']),('5. Pattern machine, set to quadruplets',r['Pattern: quadruplets (20)']),('5. Pattern machine, set to twin primes',r['Pattern: twin primes (20)']),('3. Accelerating wheels',r['Accelerating (14)']),('4. Slipping wheels',r['Slipping (17)']),('1 and 2. Wheels, gear train',r['Wheels, one tooth per tick (1-6, 8-10, 18, 21)']),('9. Mersenne finder (250 prime wheels)',fd)]
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
ringtxt=(" We also tried every one of them on each ring of up to "+str(T['all_r_upto'])+" cells that it does not fit: none came home.") if T['all_r_upto']>=16 else ""
import os,re
FA=json.load(open('tour/folded_alone_1e8.json' if os.path.exists('tour/folded_alone_1e8.json') else 'tour/folded_alone_1e7.json'))
fab='10<sup>'+FA['bound'][2:]+'</sup>'; g5=[x for x in FA['folds'] if x['r']==5][0]; fa_tot=sum(x['passed'] for x in FA['folds'])
TC=json.load(open('tour/three_check.json')); al=sorted([k for k in TC if k.startswith('alone_')],key=lambda x:float(x.split('_')[1]))[-1]
tl=('list_0' in TC and 'list_1' in TC)
t2+=one("Golden wheel: Pascal's 5 cells folded onto 2 (the Fibonacci wheel in other clothes)",2,f(g5['first'][0]),f(g5['passed'])+' below '+fab)
t2+=one("Three-cell wheel: Pascal's ring folded onto 3 cells that the number does not fit (machine 12)",3,'none','0 among all odd numbers below 10<sup>'+al.split('_')[1][2:]+'</sup>'+(('; with the doubling wheel 0 below '+B+' ('+f(TC['list_0']['tested']+TC['list_1']['tested'])+' candidates)') if tl and TC['list_0']['home']+TC['list_1']['home']==0 else ''))
t2+=one("Wilson's wheel: one cell, multiplied by 1, 2, 3, … turn after turn (machine 13)",1,'none, proved','never, but no turn can be skipped: a number of 20 digits needs more turns than can ever be made')
t2+=one("Lucas's wheel: one cell, home after n − 1 turns and at no check turn (machine 14)",1,'none, proved','never, and it is fast; but it needs the primes that divide n − 1')
t2+=one("Two cells with check turns: home after n turns, and the round trip not finished at any check turn (machine 15)",2,'none, proved','never, and it is fast; but it needs the primes that divide n + 1. For 2<sup>p</sup> − 1 it always works')
DM=open('tour/duo_mersenne.txt').read(); bm=re.search(r'p = 4423 \(1332 digits\): Lucas-Lehmer True ([\d.]+) \| two-cell wheel True ([\d.]+) \| Pascal on (\d+) cells True ([\d.]+)',DM).groups()
PW=open('tour/proved_wheels.txt').read(); w3=re.search(r'\(log2 n\)\^2 = ([\d,]+) cells',PW).group(1)
fib=re.search(r'liars of the 5-cell ring below 1e10: (\d+)',open('tour/fold_units.txt').read()).group(1)
round2="""      <li><b>Second round: fold the wheel.</b> A ring has symmetries: renumber the cells and it is the same ring. Multiply the wheel by its mirror images and the product lives on a few cells: 5 cells fold onto 2, 7 cells onto 3. A number that fools the big ring must fool the folded wheel; that is proved.%s</li>
      <li><b>Pascal's wheel on 5 cells contains the Fibonacci wheel.</b> Folded onto 2 cells it turns by the square of the golden ratio. So every mistake of the 5 cells is a mistake of machine 6 and of machine 7 too (all %s below 10<sup>10</sup> are), and a number that does not fit 5 cells and still came home would win the 620 dollars.</li>
      <li><b>The rule holds one level down.</b> Folded wheels alone made %s mistakes below %s, every one by a number that moves the cells round with order 2: it fits. Where a number must move three cells round, nothing came home. That is machine 12, the three-cell wheel.</li>
      <li><b>On Mersenne numbers the two-cell wheel is proved.</b> For 2<sup>p</sup> − 1 its step is always 3, and home is exactly the Lucas–Lehmer test.</li>
      <li><b>Proved or fast, and both only for special numbers.</b> Wilson's wheel (machine 13) never errs and that is proved, but it needs n − 1 turns and none can be skipped. Lucas's wheel (machine 14) is proved and fast, but only when the primes that divide n − 1 are known; that is how numbers such as 3 × 2<sup>2208</sup> + 1 are proved. Machine 15 does the same with n + 1 and proves every Mersenne prime: type 2^4423-1 there. The proved ring test of Agrawal, Kayal and Saxena needs more than %s cells for 2<sup>4423</sup> − 1, where machine 10 uses 5. In one timing run on that number Lucas–Lehmer took %s seconds, the two-cell wheel %s and Pascal's wheel %s. For every number there is a method that is proved and works in practice, but it is not one wheel: elliptic curves give wheels whose size can be chosen, and one searches for a wheel whose size comes apart; the largest such proof we found reported has 109,297 digits. One small wheel with few turns, for every number and with a proof, is what is missing.</li>
"""%(((" No number below "+B+" that fools the doubling wheel passed a folded wheel it does not fit (rings of "+", ".join(sorted(T['folded'],key=int))+" cells), so a large ring can be checked by turning a small wheel.") if T['folded'] and T['fold_pass_total']==0 else ""),fib,f(fa_tot),fab,w3,bm[0],bm[1],bm[3])
html='''  <section class="panel" id="compare">
    <h2>Which machine has the best results</h2>
    <p class="note">Every machine was run for 100,000 ticks with 250 wheels. The table gives the most prime wheels, Mersenne-exponent wheels and twin-prime wheels that were home at once, and the tick. Among the numbers 1 to 250 there are 53 primes, 12 Mersenne exponents and 33 twin primes.</p>
    <div class="boardwrap"><table>
      <thead><tr><th>Machine</th><th>Primes at home</th><th>Mersenne exponents at home</th><th>Twin primes at home</th></tr></thead>
      <tbody>%s</tbody>
    </table></div>
    <ul class="found">
      <li><b>The wheels that decide win without a trick.</b> On machines 10 to 15 each lamp is a decision: wheel n goes home at tick n only when n is prime, and stays there. At tick 241 all 53 prime lamps are on, and no other lamp.</li>
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
      <li><b>Is this the gear that amazes mathematicians?</b> Not yet, and we say so plainly. The last three rows are open problems. Nobody has proved that they never make a mistake. For the pair of wheels there is a prize of 620 dollars, offered by Pomerance, Selfridge and Wagstaff; for Pascal's wheel the claim is a sibling of Agrawal's conjecture (his has x − 1 where the wheel has x + 1; on a ring with an even number of cells the two are the same); the two-cell wheel is a Frobenius test, for which no mistake is known either. There are reasons to expect a mistake at some very large number (Pomerance for the pair; Lenstra and Pomerance for Pascal's wheel). A proof for any of the three, or the first number that fools one, would be news.</li>
%s      <li><b>What they can do today.</b> Both wheels need a few cells and about as many steps as the number has binary digits, so they also run on numbers of a thousand digits: choose machine 10 or 11 and type one.</li>
    </ul>
  </section>
'''%(tb,t2,f(T['fit_liars']),f(T['list_count']),B,ringtxt,round2)
s=open('prime3.html',encoding='utf-8').read()
i=s.index('  <section class="panel" id="compare">'); j=s.index('</section>\n',i)+len('</section>\n'); s=s[:i]+html+s[j:]
import re
s,k=re.subn(r'In our runs none of the [\d,]+ odd numbers below [\d,]+ that are not prime came home','In our runs none of the %s odd numbers below %s that are not prime came home'%(f(S['duo']['tested']),f(S['duo']['bound'])),s); assert k==1
tri="In our runs of the three cells, %s odd numbers below %s that are not prime were tried and none came home"%(f(TC[al]['tested']),f(int(float(al.split('_')[1]))))
if tl and TC['list_0']['home']+TC['list_1']['home']==0: tri+="; nor did any of the %s numbers below 2^64 that fool the doubling wheel"%f(TC['list_0']['tested']+TC['list_1']['tested'])
assert TC[al]['home']==0
s,k=re.subn(r'In our runs of the three cells, [\d,]+ odd numbers below [\d,]+ that are not prime were tried and none came home(; nor did any of the [\d,]+ numbers below 2\^64 that fool the doubling wheel)?',tri,s); assert k==1
open('prime3.html','w',encoding='utf-8').write(s); print(len(s))
