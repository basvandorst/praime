import json
r={x[0]:x for x in json.load(open('analysis/best.json'))['rows']}
fd=json.load(open('analysis/best2.json'))['finder']
T=json.load(open('tour/tournament.json'))
f=lambda n:format(n,',')
pas=['pascal']+T['row250']
duo=['duo']+T['duo_row250']
R12=json.load(open('tour/rows_12_13.json'))
CH=json.load(open('tour/chain.json'))
RC=json.load(open('tour/receipt.json'))
N2=json.load(open('tour/norm2.json'))
ES=json.load(open('tour/early_sum.json')); IS=json.load(open('tour/inherit_sum.json'))
rows=[("10. Pascal's wheel",pas),("11. Two-cell wheel",duo),("12. Three-cell wheel",['tri']+R12['three']),("13. Wilson's wheel",['wil']+R12['wilson']),("14. Lucas's wheel",['luc']+R12['lucas']),("15. Two cells with check turns",['np']+R12['nplus']),("16. Curve wheel",['ec']+json.load(open('tour/curve_wheel.json'))['row250']),("17. Chain of wheels",['ch']+CH['row250']),("19. Two cells from a start of norm 2",['n2']+N2['row250']),("23. Wheels that light each other",['inh']+IS['row250']),("21. Gears from its own lamps",['sf']+ES['row250_sf']),("20. Squaring wheels",['sq']+ES['row250_sq']),("22. Counting wheels",['cnt']+ES['row250_cnt']),("18. Two cells with a receipt (wheel p decides 2<sup>p</sup> − 1)",['rc']+RC['row250']),('8. Back and forth',r['Back and forth (16)']),('6. Fibonacci wheels',r['Fibonacci (12)']),('7. Doubling wheels',r['Doubling (11)']),('5. Pattern machine, set to quadruplets',r['Pattern: quadruplets (20)']),('5. Pattern machine, set to twin primes',r['Pattern: twin primes (20)']),('3. Accelerating wheels',r['Accelerating (14)']),('4. Slipping wheels',r['Slipping (17)']),('1 and 2. Wheels, gear train',r['Wheels, one tooth per tick (1-6, 8-10, 18, 21)']),('9. Mersenne finder (250 prime wheels)',fd)]
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
t2+=one("Curve wheel: the points of a curve form a wheel whose size can be chosen (machine 16)",'1 point','none, proved','never, and for every number; on this page up to 3,000,000 only, because it counts the positions one by one')
t2+=one("Chain of wheels: one cell, two cells or a curve wheel, whichever size comes apart, and then the same for its large prime part (machine 17)",'1, 2 or 1 point','none, proved','never, it is fast, and it proved every number we gave it up to %d digits; but it is a search for the right wheels, not one wheel'%CH['js_maxdigits'])
big=lambda d,mode:[x for x in CH['big'] if x['digits']==d and x['mode'].startswith(mode)][0]
j100=[x for x in CH['js'] if x.get('digits')==100][0]; pm60=big(60,'n - 1 and n + 1 only'); f9=big(100,'all sizes of wheel, 9 curve')
round3='''      <li><b>Third round: a chain of wheels has all three.</b> A number n has more wheels than one: one cell with n − 1 positions, two cells with a round trip of n + 1 turns, and curve wheels. For %d families of curves the number of positions follows from n at once, without counting. Machine 17 takes a size that comes apart, checks that wheel, and goes on with the large prime part: a chain of wheels down to a small number. The page's code proved all %d random primes of 20 to %d digits that we gave it, in about %.1f seconds for 100 digits, and a separate program accepted every chain. With n − 1 and n + 1 alone only %d of %d primes of 60 digits were proved, and with 9 families of curves %d of %d primes of 100 digits. It built no chain for any of the %s numbers that are not prime and pass the doubling wheel which we tried. This is not our idea: it is elliptic curve primality proving (Atkin and Morain), cut down to fit a page. And it is a search, with no proof that it always ends soon.</li>
'''%(CH['families'],CH['js_proved'],CH['js_maxdigits'],int(j100['median_ms'])/1000,pm60['proved'],pm60['primes'],f9['proved'],f9['primes'],f(CH['adv']['spsp']['count']+sum(sum(CH['adv'][k].values()) for k in CH['adv'] if k.startswith('pq'))))
assert CH['js_proved']==CH['js_primes']
t2+=one("Two cells from a start of norm 2: the start a and 1 with a × a two more or two less than the step (machine 19)",2,'none','0 below '+B+', alone: a number that fools it must fool the doubling wheel, and none of the '+f(N2['list']['lines'])+' that do came home')
bp=N2['bench']['primes']
cen_txt=''
if os.path.exists('tour/norm2_census.json'):
    CJ=json.load(open('tour/norm2_census.json')); ct=CJ['total']; assert ct['misfit_home']==0
    cen_txt=" And it does not hang on which start is taken: all %d starts with a up to 8 were tried on every number of the list. In %s tests on a start the number does not fit, none came home; in %s tests on a start it fits, %s behaved like a prime. A number cheats only on a wheel it fits."%(len(CJ['rows']),f(ct['misfit']),f(ct['fit']),f(ct['fit_same']))
round5='''      <li><b>Fifth round: two cells, alone, settled below 2<sup>64</sup>.</b> Pascal's wheel is settled alone below 2<sup>64</sup>, because the sum of its cells is the doubling wheel; but it needs 5 to 47 cells. Machine 11 has two cells, but was checked up to 2<sup>64</sup> only with the doubling wheel at its side. Machine 19 has both. It changes the start of the two cells to a and 1, with a × a two more or two less than the step. Then the two cells, multiplied into one number, give 2 or −2, and that number is doubled at every turn. So whatever fools this wheel fools the doubling wheel too (proved). All %s numbers below 2<sup>64</sup> that fool the doubling wheel were put on it: none came home. As a check without the list, every odd number below 10<sup>10</sup> went on the wheel: %s came home, and a sieve counts %s odd primes there. On primes of 19 digits it took %s billionths of a second, against %s for the doubling wheel alone and %s for the seven tests that are known to settle every number below 2<sup>64</sup>. Khashin reported in 2018 that his version of the two-cell test makes no mistake below 2<sup>64</sup>, by a long computation; the short reason here is the 2. Above 2<sup>64</sup> nothing is proved.%s</li>
'''%(f(N2['list']['lines']),f(N2['range']['home']),f(N2['range']['primes']),f(bp['norm-2 wheel']['ns']),f(bp['doubling wheel']['ns']),f(bp['seven strong tests']['ns']),cen_txt)
assert N2['list']['home']==0 and N2['range']['home']==N2['range']['primes']
two=[x for x in RC['timing'] if x['group']=='two cells']; t0=two[0]; t9=two[-1]; fo={(x['group'],x['factor'][:2]):x for x in RC['forgery']}
f1=fo[('one','-1')]; f4=[x for x in RC['forgery'] if x['group']=='two' and 'order 4' in x['factor']][0]; tp=RC['tamper']['two']; sf=RC['safe']
round4='''      <li><b>Fourth round: a test with a receipt.</b> A different way to win on the largest numbers: not a faster test, but a test that need not be done twice. Machine 18 runs the Lucas–Lehmer test in its two-cell form and writes down a few numbers, the receipt. With it a checker confirms the result by turning only a small piece of the strip of turns. In our runs the check was %.0f times faster than the test at p = %s and %.0f times at p = %s (table below). Every wrong result (%d of %d), every changed number in a receipt (%d of %d) and every bit flipped during a test (%d of %d) was rejected.</li>
      <li><b>A receipt is not a proof of the result.</b> We tried to cheat. A tester who has made the turns of a prime 2<sup>p</sup> − 1 can multiply the true result by a number that returns to 1 after a few multiplications, and the check still passes. On the one-cell form of the test the factor −1 was accepted in %d of %d attempts; on two cells a factor of that kind was accepted in %d of %d. Each time the verdict changed from prime to not prime. So a receipt shows that the work was done, not that the reported result is the true one. This weakness is known (Hoffmann, Hubáček, Kamath and Pietrzak, 2023, who also repair it). A simple guard worked in our runs: compare a power of the result that wipes out such factors. All %d forged results were then sent to a full test, and none of %d numbers that are not prime was sent without need.</li>
'''%(t0['ratio'],f(t0['p']),t9['ratio'],f(t9['p']),tp['wrong_y_tried']-tp['wrong_y_accepted'],tp['wrong_y_tried'],tp['wrong_mu_tried']-tp['wrong_mu_accepted'],tp['wrong_mu_tried'],tp['flips']-tp['flips_accepted'],tp['flips'],f1['accepted'],f1['tried'],f4['accepted'],f4['tried'],sf['one_forged']+sf['two_forged'],sf['composites'])
rct="".join("<tr><td>%s</td><td>%s</td><td>%.2f s</td><td>%.2f s</td><td>%s and %d</td><td>%.3f s</td><td>%.0f</td></tr>"%(f(x['p']),f(x['squarings']),x['test_s'],x['receipt_s'],f(x['check_squarings']),x.get('small_powers',2*x['k']),x['check_s'],x['ratio']) for x in two)
DM=open('tour/duo_mersenne.txt').read(); bm=re.search(r'p = 4423 \(1332 digits\): Lucas-Lehmer True ([\d.]+) \| two-cell wheel True ([\d.]+) \| Pascal on (\d+) cells True ([\d.]+)',DM).groups()
PW=open('tour/proved_wheels.txt').read(); w3=re.search(r'\(log2 n\)\^2 = ([\d,]+) cells',PW).group(1)
fib=re.search(r'liars of the 5-cell ring below 1e10: (\d+)',open('tour/fold_units.txt').read()).group(1)
round2="""      <li><b>Second round: fold the wheel.</b> A ring has symmetries: renumber the cells and it is the same ring. Multiply the wheel by its mirror images and the product lives on a few cells: 5 cells fold onto 2, 7 cells onto 3. A number that fools the big ring must fool the folded wheel; that is proved.%s</li>
      <li><b>Pascal's wheel on 5 cells contains the Fibonacci wheel.</b> Folded onto 2 cells it turns by the square of the golden ratio. So every mistake of the 5 cells is a mistake of machine 6 and of machine 7 too (all %s below 10<sup>10</sup> are), and a number that does not fit 5 cells and still came home would win the 620 dollars.</li>
      <li><b>The rule holds one level down.</b> Folded wheels alone made %s mistakes below %s, every one by a number that moves the cells round with order 2: it fits. Where a number must move three cells round, nothing came home. That is machine 12, the three-cell wheel.</li>
      <li><b>On Mersenne numbers the two-cell wheel is proved.</b> For 2<sup>p</sup> − 1 its step is always 3, and home is exactly the Lucas–Lehmer test.</li>
      <li><b>Proved or fast, and both only for special numbers.</b> Wilson's wheel (machine 13) never errs and that is proved, but it needs n − 1 turns and none can be skipped. Lucas's wheel (machine 14) is proved and fast, but only when the primes that divide n − 1 are known; that is how numbers such as 3 × 2<sup>2208</sup> + 1 are proved. Machine 15 does the same with n + 1 and proves every Mersenne prime: type 2^4423-1 there. The proved ring test of Agrawal, Kayal and Saxena needs more than %s cells for 2<sup>4423</sup> − 1, where machine 10 uses 5. In one timing run on that number Lucas–Lehmer took %s seconds, the two-cell wheel %s and Pascal's wheel %s. For every number there is a method that is proved and works in practice, but it is not one wheel: a curve gives wheels whose size can be chosen, and one searches for a wheel whose size comes apart. Machine 16 is a slow version of it and machine 17 a fast one; the largest such proof we found reported has 109,297 digits. One small wheel with few turns, for every number and with a proof, is what is missing.</li>
"""%(((" Below "+B+" the numbers that fool the doubling wheel were put on folded wheels "+f(sum(q['tested'] for d in T['folded'].values() for q in d['folds'].values()))+" times (rings of "+", ".join(sorted(T['folded'],key=int))+" cells): "+("none passed." if T['fold_pass_total']==0 else str(T['fold_pass_total'])+" passed, on a fold onto 2 cells that it fits one level down, and the full ring turned it away.")+" So a large ring can be checked by turning a small wheel.") if T['folded'] else ""),fib,f(fa_tot),fab,w3,bm[0],bm[1],bm[3])
html='''  <section class="panel" id="compare">
    <h2>Which machine has the best results</h2>
    <p class="note">Every machine was run for 100,000 ticks with 250 wheels. The table gives the most prime wheels, Mersenne-exponent wheels and twin-prime wheels that were home at once, and the tick. Among the numbers 1 to 250 there are 53 primes, 12 Mersenne exponents and 33 twin primes.</p>
    <div class="boardwrap"><table>
      <thead><tr><th>Machine</th><th>Primes at home</th><th>Mersenne exponents at home</th><th>Twin primes at home</th></tr></thead>
      <tbody>%s</tbody>
    </table></div>
    <ul class="found">
      <li><b>The wheels that decide win without a trick.</b> On machines 10 to 17 and 19 each lamp is a decision: wheel n goes home at tick n only when n is prime, and stays there. At tick 241 all 53 prime lamps are on, and no other lamp. Machine 18 asks another question (is 2<sup>p</sup> − 1 prime?), so only the 12 Mersenne exponents light up.</li>
@@EARLY_LI@@      <li><b>Back and forth also has them all, for a simple reason.</b> At tick 36 its wheels are exactly back where they started, so every wheel is home, prime or not.</li>
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
      <li><b>The rule that wins: give the number a wheel it does not fit.</b> A number n fits a ring of r cells when n × n leaves remainder 1 after division by r, and it fits a two-cell wheel with step c when c is a square in the arithmetic of remainders after division by n. Every mistake of Pascal's wheel that we found is a number that fits its ring: %s on rings of 2 to 40 cells below 10<sup>10</sup>, not one that does not fit. For rings the rule is not ours: it is the exception that Popovych wrote into his version of Agrawal's conjecture in 2009. Here it is carried over to the wheels of two and three cells.</li>
      <li><b>Checked as far as the lists go.</b> A number that fools Pascal's wheel must fool the doubling wheel first (add up the cells: the sum doubles every turn). All %s odd numbers below %s that fool the doubling wheel are known, and we put every one on the smallest ring it does not fit and on its two-cell wheel: none came home.%s</li>
      <li><b>Is this the gear that amazes mathematicians?</b> Not yet, and we say so plainly. The rows that say “none” without “proved” are open problems: the pair of wheels and machines 10, 11, 12 and 19. Nobody has proved that they never make a mistake. For the pair of wheels there is a prize of 620 dollars, offered by Pomerance, Selfridge and Wagstaff; for Pascal's wheel the claim is a sibling of Agrawal's conjecture (his has x − 1 where the wheel has x + 1; on a ring with an even number of cells the two are the same); the two-cell wheel is a Frobenius test, for which no mistake is known either (Khashin reports that his own version of it makes none below 2<sup>64</sup>). There are reasons to expect a mistake at some very large number (Pomerance for the pair; Lenstra and Pomerance for Pascal's wheel). A proof for any of them, or the first number that fools one, would be news.</li>
%s%s%s%s      <li><b>What they can do today.</b> Both wheels need a few cells and about as many steps as the number has binary digits, so they also run on numbers of a thousand digits: choose machine 10 or 11 and type one.</li>
    </ul>
@@EARLY@@    <h2 id="winners">Who wins where</h2>
    <p class="note">After six rounds, the best machine on this page for each question, and what it cannot do.</p>
    <div class="boardwrap"><table class="dec">
      <thead><tr><th>Question</th><th>Winner here</th><th>Why</th><th>Limit</th></tr></thead>
      <tbody>
        <tr><td>Is n prime, for n below 2<sup>64</sup>?</td><td>Machine 19: one wheel of two cells</td><td>Nothing below 2<sup>64</sup> fools it. Every candidate is on a published list, and none came home</td><td>Above 2<sup>64</sup> nothing is proved</td></tr>
        <tr><td>Is n prime, with a proof, for n of any size?</td><td>Machine 17: a chain of wheels</td><td>Every link is a proved test. It proved all %d random primes of 20 to %d digits we gave it</td><td>It is a search, and this page stops at 150 digits</td></tr>
        <tr><td>Is 2<sup>p</sup> − 1 prime?</td><td>The Lucas–Lehmer test: machines 11, 15 and 18 are that test in gear form</td><td>Proved, and one squaring per turn. No wheel here needs fewer turns</td><td>The p − 2 turns cannot be skipped by anything we know</td></tr>
        <tr><td>Was a long test done right?</td><td>Machine 18: the receipt</td><td>The check was %.0f times faster than the test at p = %s</td><td>It shows that the work was done, not that the reported result is the true one</td></tr>
@@WIN_ROWS@@      </tbody>
    </table></div>
    <h2 id="proved">Proved, fast, every number: which wheel has what</h2>
    <p class="note">Fast means that the wheel can jump: a few hundred steps do the work of n turns, so a number of a thousand digits is no problem. No single wheel on this page has all three for numbers of every size. The chain of machine 17 has them, and below 2<sup>64</sup> the wheels of machines 10 and 19 have them.</p>
    <div class="boardwrap"><table class="dec">
      <thead><tr><th>Machine</th><th>Proved</th><th>Fast</th><th>Every number</th></tr></thead>
      <tbody>
        <tr><td>7. Doubling wheel</td><td>no: 341 fools it</td><td>yes</td><td>yes</td></tr>
        <tr><td>10. Pascal's wheel</td><td>below 2<sup>64</sup> yes, given the published list; above that no, and no mistake is known</td><td>yes</td><td>yes</td></tr>
        <tr><td>11. Two-cell wheel</td><td>no; no mistake found (below 2<sup>64</sup>: none together with the doubling wheel). Proved for 2<sup>p</sup> − 1</td><td>yes</td><td>yes</td></tr>
        <tr><td>12. Three-cell wheel</td><td>no; no mistake found (below 2<sup>64</sup>: none together with the doubling wheel)</td><td>yes</td><td>yes</td></tr>
        <tr><td>13. Wilson's wheel</td><td>yes</td><td>no: every turn is needed</td><td>yes</td></tr>
        <tr><td>14. Lucas's wheel</td><td>yes</td><td>yes</td><td>no: only when n − 1 comes apart</td></tr>
        <tr><td>15. Two cells with check turns</td><td>yes</td><td>yes</td><td>no: only when n + 1 comes apart (every 2<sup>p</sup> − 1 does)</td></tr>
        <tr><td>16. Curve wheel</td><td>yes</td><td>no: it counts its positions one by one</td><td>yes</td></tr>
        <tr><td>17. Chain of wheels</td><td>yes</td><td>yes</td><td>yes for every number we tried, up to %d digits here; it is a search, so there is no guarantee</td></tr>
        <tr><td>19. Two cells from a start of norm 2</td><td>below 2<sup>64</sup> yes, given the published list; above that no, and no mistake is known</td><td>yes</td><td>yes</td></tr>
        <tr><td>18. Two cells with a receipt</td><td>yes (Lucas–Lehmer), and the work can be checked without doing it again</td><td>yes</td><td>no: 2<sup>p</sup> − 1 only</td></tr>
@@PROVED_ROWS@@      </tbody>
    </table></div>
    <h2 id="receipt">A test with a receipt: what the check costs</h2>
    <p class="note">Machine 18 on Mersenne primes 2<sup>p</sup> − 1, outside the page (Python, the smaller time of two runs). The receipt is 6 or 8 numbers, whichever made the check faster. A small power is a power with a number of 64 binary digits.</p>
    <div class="boardwrap"><table class="dec">
      <thead><tr><th>p</th><th>Turns of the test</th><th>Test</th><th>Writing the receipt</th><th>Turns and small powers of the check</th><th>Check</th><th>Times faster</th></tr></thead>
      <tbody>%s</tbody>
    </table></div>
  </section>
'''%(tb,t2,f(T['fit_liars']),f(T['list_count']),B,ringtxt,round2,round3,round4,round5,CH['js_proved'],CH['js_maxdigits'],t9['ratio'],f(t9['p']),CH['js_maxdigits'],rct)
ic=IS['c']; firsts=[x[2] for x in ic['cnt']]; icm=IS['c_minus']; law=IS['law']; sfx=ES['early']['self']; wk=ES['work']; dy={t[0]:t[2] for t in ES['dyadic']}
def sf_front(T):
    return sfx['rows'][T-1][3]
def sf_lamps(T):
    return sfx['rows'][T-1][5]
big60=[b for b in IS['big'] if b['digits']==60][0]; big8=[b for b in IS['big'] if b['digits']==8][0]
assert firsts[1:]==[3,19,173,21767,1688173,2223990533] and ic['cnt'][-1][0]==6
cross=ES['early']['cross']; box=ES['box']
early_li='      <li><b>Four machines do not wait for tick n.</b> Machines 20, 21 and 23 have all 53 prime lamps on at tick %d, %d and %d, and machine 22 counts ahead of its lamps. They are compared under <a href="#early">primes before their tick</a>.</li>\n'%(ES['row250_sq'][1],ES['row250_sf'][1],IS['row250'][1])
erow=lambda *c:"<tr>"+"".join("<td>%s</td>"%x for x in c)+"</tr>"
etab=erow('10 to 19: wheel n makes one turn per tick','4','6','11','31','4,999','one turn of each wheel')
etab+=erow('21: a new prime gear, taken from the lamps',f(sf_front(3)),f(sf_front(5)),f(sf_front(10)),f(sf_front(30)),str(sfx['all_5000'][0]),'it looks at every number below the frontier: %.0f moves per prime with 5,000 wheels'%((wk['strikes']+wk['looked'])/wk['primes']))
etab+=erow('20: every wheel squares its cells','16','64',f(2**11),f(2**31),'12','every wheel works at every tick: %.0f ticks of a wheel per prime with 5,000 wheels'%(wk['ladder_ticks']/wk['primes']))
etab+=erow('23: a round in which every dark wheel asks the lit ones',f(firsts[4]),f(firsts[6]),'not known','not known',str(IS['upto5000'][-1][0]),'a round is heavy: every dark wheel tries in every round')
etab+=erow('22: a gear changes counters (it counts, it does not find)','25','49','121','961','4,999 for the lamps','%.2f moves per prime counted with 5,000 gears'%(wk['count_moves']/wk['count_primes']))
assert sfx['all_5000'][0]==19 and (5000).bit_length()-1==12 and IS['upto5000'][-1][0]==3
early_html='''    <h2 id="early">Primes before their tick</h2>
    <p class="note">On machines 10 to 19 wheel n decides at tick n, so 30 ticks give the 10 primes up to 30. Four machines break that rule, each in its own way. The table says below which number every prime is on (for machine 22: counted) after 3, 5, 10 and 30 ticks, and at which tick all 669 primes below 5,000 are on.</p>
    <div class="boardwrap"><table class="dec">
      <thead><tr><th>Machine, and what one tick is</th><th>After tick 3: every prime below</th><th>After tick 5</th><th>After tick 10</th><th>After tick 30</th><th>All 669 primes below 5,000 at tick</th><th>The price</th></tr></thead>
      <tbody>%s</tbody>
    </table></div>
    <ul class="found">
      <li><b>Sixth round: do not wait for the tick.</b> Machine 21 looks ahead: with the gears it has, every number below the square of the next gear is already decided. Machine 20 turns faster: a wheel that squares its cells is done when the binary digits of its number are used up. Machine 23 lets the wheels help each other. Machine 22 does not find the primes at all: it counts them. With 5,000 wheels the lamps are all on at tick %d, %d and %d instead of 4,999.</li>
      <li><b>Squaring is as early as one wheel can be.</b> Let a tick be one multiplication of values the wheel already holds and at most one plain turn. Then the turns made can at most double and grow by one, so a wheel that has to make n turns cannot be done before tick ⌊log<sub>2</sub> n⌋ (proved). Machine 20 is done exactly then: tick T decides every wheel below 2<sup>T+1</sup>. Its verdicts are those of machine 19, so with the published list there is no mistake below 2<sup>64</sup>, and at tick 63 all 425,656,284,035,217,743 primes below 2<sup>64</sup> would be on. Checked without the list: every odd number below 2<sup>32</sup> went on the wheel, and with the prime 2 that makes %s primes, the number machine 22 counts there.</li>
      <li><b>A sieve that makes its own gears, and an open question inside it.</b> Machine 21 is the sieve of Eratosthenes counted in gears: 10 gears decide every number below %s (%s primes), 30 gears every number below %s (%s primes). Does every tick after the first light at least four new lamps? That is Brocard's conjecture of 1904, and nobody has proved it. In our run of %d ticks the fewest were %d, at tick %d.</li>
      <li><b>More primes than moves, if counting is enough.</b> Machine 22 with the gears up to 5,000 ends at %s primes up to %s after %s changes of a counter: %.1f primes per move. From %d wheels on it always counts more primes than it makes moves (tried up to 2,999 wheels). On the page it counted the %s primes up to 10<sup>10</sup> in %s moves, %d primes per move, without looking at any number above 100,000. The price is plain: it knows how many primes there are, not which. The rule is Legendre's of 1808, in the form known as Lucy's recursion.</li>
      <li><b>A proof does not have to go deep when a number may ask both neighbours.</b> On machine 23 a wheel is proved by lit wheels that divide n − 1 or n + 1, and each of those was proved the same way, down to 2. Up to %s no prime needs more than %d rounds. The smallest primes that need 1 to 6 rounds are %s. Random primes of %d digits needed %.1f rounds on average, of 60 digits %.1f (never more than %d). All our sizes from 1,000 up to 60 digits fit one rule of thumb: rounds = log<sub>4</sub>(ln n) + %.2f, within %.2f. The reason we give for it is not a proof: when both neighbours help, each needs to supply only a quarter of the digits, so one more round raises the reachable size to about its fourth power. With n − 1 alone the rule becomes log<sub>2</sub>(ln n), and more rounds are needed: %s against %s on average for the primes between 2<sup>28</sup> and 2<sup>29</sup>, and up to %d below %s.@@UNSEEN@@</li>
      <li><b>What is ours here and what is not.</b> Squaring, the sieve and the counting are old and well tested; we only put them next to each other on one clock. The two proof steps of machine 23 are old too (Pocklington 1914; Morrison, and Brillhart, Lehmer and Selfridge, 1975). The count of rounds with both neighbours, its list of first primes and the rule log<sub>4</sub>(ln n) we did not find in a search; the list is not in the OEIS. Close relatives are known: the classes of Erdős and Selfridge and the height of the Pratt tree, which grows like 2.7 × ln ln n in the model of Ford, Konyagin and Luca.</li>
    </ul>
'''%(etab,12,sfx['all_5000'][0],IS['upto5000'][-1][0],f(dy[31]),f(sf_front(10)),f(sf_lamps(10)),f(sf_front(30)),f(sf_lamps(30)),sfx['ticks'],sfx['min_new'][0],sfx['min_new'][1],f(wk['count_primes']),f(wk['N']**2),f(wk['count_moves']),wk['count_primes']/wk['count_moves'],cross,f(box[1]),f(box[2]),round(box[1]/box[2]),f(ic['X']),ic['cnt'][-1][0],', '.join(f(x) for x in firsts[1:]),big8['digits'],big8['mean'],big60['mean'],big60['worst'][0],law['both'][2],max(law['both'][2]-law['both'][0],law['both'][1]-law['both'][2]),'%.2f'%dict((a,c) for a,b,c in icm['mean'])[28],'%.2f'%dict((a,c) for a,b,c in ic['mean'])[28],icm['cnt'][-1][0],f(icm['X']))
win_rows='''        <tr><td>How many primes are there up to x?</td><td>Machine 22: counting wheels</td><td>It counted the %s primes up to 10<sup>10</sup> in %s moves: %d primes per move</td><td>It knows how many, not which</td></tr>
        <tr><td>Which machine has its lamps on first?</td><td>Machine 23, then 20 and 21</td><td>All 669 primes below 5,000 are on after %d rounds, 12 ticks and %d ticks, instead of 4,999 ticks</td><td>A tick is heavier than one turn, and the total work is not less</td></tr>
'''%(f(box[1]),f(box[2]),round(box[1]/box[2]),IS['upto5000'][-1][0],sfx['all_5000'][0])
proved_rows='''        <tr><td>20. Squaring wheels</td><td>as machine 19: below 2<sup>64</sup> yes, given the published list</td><td>yes: done at tick log<sub>2</sub> n</td><td>yes</td></tr>
        <tr><td>21. Gears from its own lamps</td><td>yes (the sieve of Eratosthenes)</td><td>no: it looks at every number below the frontier</td><td>yes</td></tr>
        <tr><td>22. Counting wheels</td><td>yes (Legendre's rule), for the count</td><td>yes for counting: far fewer moves than primes</td><td>it does not say which numbers are prime</td></tr>
        <tr><td>23. Wheels that light each other</td><td>yes</td><td>yes when n − 1 and n + 1 can be taken apart</td><td>every prime up to %s in at most %d rounds; for a large number taking the neighbours apart is the hard part</td></tr>
'''%(f(ic['X']),ic['cnt'][-1][0])
_un=IS.get('unseen',[]); _r0=IS['deep']['stage2']['rows'][0]
unseen_txt=((" The rule was then tried on sizes it was not made from: "+"; ".join("%d random primes of %d digits needed %.3f rounds on average, where the rule says %.3f"%(u['count'],u['digits'],u['mean'],u['pred']) for u in _un)+".") if _un else "")+" Primes that need 7 rounds exist: the smallest we found is %s, by building on primes of 6 rounds; a smaller one may exist%s, and it is not proved that the number of rounds has no upper limit."%(f(_r0[0]),(", but not below %s"%f(IS['second8']['no7'])) if IS.get('second8') else "")
_src={d['k']:d for d in IS.get('sources',[])}
if 3 in _src: unseen_txt+=" The reason was tested as well: it says that a third source of helpers, the number n × n + 1, should change the rule to log<sub>6</sub>(ln n). Counted afterwards, the rounds grew with slope %.2f where the reason says %.2f (with one source %.2f against %.2f, with two %.2f against %.2f)."%(_src[3]['slope'],_src[3]['predicted_slope'],_src[1]['slope'],_src[1]['predicted_slope'],_src[2]['slope'],_src[2]['predicted_slope'])
if 5 in _src: unseen_txt+=" With five sources a straight line did not fit (slope %.2f where the reason says %.2f): the curve is a step, flat at about 2 rounds up to two million and rising after that. Over the whole range of sizes the mean rose by %.2f rounds, where the reason says %.2f."%(_src[5]['slope'],_src[5]['predicted_slope'],_src[5]['rise'],_src[5]['pred_rise'])
early_html=early_html.replace('@@UNSEEN@@',unseen_txt)
html=html.replace('@@EARLY_LI@@',early_li).replace('@@EARLY@@',early_html).replace('@@WIN_ROWS@@',win_rows).replace('@@PROVED_ROWS@@',proved_rows); assert '@@' not in html
s=open('prime3.html',encoding='utf-8').read()
i=s.index('  <section class="panel" id="compare">'); j=s.index('</section>\n',i)+len('</section>\n'); s=s[:i]+html+s[j:]
import re
s,k=re.subn(r'In our runs none of the [\d,]+ odd numbers below [\d,]+ that are not prime came home','In our runs none of the %s odd numbers below %s that are not prime came home'%(f(S['duo']['tested']),f(S['duo']['bound'])),s); assert k==1
tri="In our runs of the three cells, %s odd numbers below %s that are not prime were tried and none came home"%(f(TC[al]['tested']),f(int(float(al.split('_')[1]))))
if tl and TC['list_0']['home']+TC['list_1']['home']==0: tri+="; nor did any of the %s numbers below 2^64 that fool the doubling wheel"%f(TC['list_0']['tested']+TC['list_1']['tested'])
assert TC[al]['home']==0
s,k=re.subn(r'In our runs of the three cells, [\d,]+ odd numbers below [\d,]+ that are not prime were tried and none came home(; nor did any of the [\d,]+ numbers below 2\^64 that fool the doubling wheel)?',tri,s); assert k==1
open('prime3.html','w',encoding='utf-8').write(s); print(len(s))
