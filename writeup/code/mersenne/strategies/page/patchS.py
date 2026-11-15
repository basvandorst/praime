"""Insert the new strategy layouts (ls.js, lsnotes.js) into prime2.html. usage: patchS.py in.html out.html"""
import sys
src, dst = sys.argv[1], sys.argv[2]
t = open(src, encoding='utf-8').read()
ls = open('ls.js', encoding='utf-8').read().strip()
notes = open('lsnotes.js', encoding='utf-8').read().strip()
for s in (ls, notes):
    assert '//' not in s.replace('://', '') and '/*' not in s, 'no comments allowed'
def rep(old, new, count=1):
    global t
    assert t.count(old) == count, (old[:60], t.count(old))
    t = t.replace(old, new)
OPTS = [('quadform', 'Sums a² + P·b² (Fermat, Gauss)'), ('eisenstein', 'Sums a² + ab + b² (hexagonal)'), ('mersieve', 'Mersenne factor sieve (row k: 2kp + 1)'),
        ('gaptowers', 'Gap towers (a column per prime)'), ('polyspiral', 'Polygon spiral (P sides)'), ('digitrows', 'Rows by digit sum (base P)'),
        ('chebrace', 'Race of squares mod P (walk)'), ('totient', 'Totient: height φ(n)/n'), ('goldbach', "Goldbach's comet"), ('dickman', 'Largest prime factor (Dickman)'), ('primroot', 'Least primitive root (Artin)'), ('pili', 'π(n) − li(n) (walk)'), ('ord2', 'Order of 2 mod q (Mersenne factors)')]
rep('<option value="twinwave">Wave that grows with the twin primes</option></select>',
    '<option value="twinwave">Wave that grows with the twin primes</option>' + ''.join(f'<option value="{v}">{l}</option>' for v, l in OPTS) + '</select>')
rep('function sqW(){', ls + '\n' + notes + '\nfunction sqW(){')
rep('function sqW(){return Math.ceil(Math.sqrt(S.N*(LWR.includes(S.layout)?4:1)));}', 'function sqW(){return Math.ceil(Math.sqrt(S.N*((LWR.includes(S.layout)||LSW.includes(S.layout))?4:1)));}')
rep('$(id).disabled=!["rows","snake","wheel","honeycomb",...LWR].includes(S.layout));$("period").disabled=!LWN.includes(S.layout)||S.layout==="wave210";}',
    '$(id).disabled=!["rows","snake","wheel","honeycomb",...LWR,...LSW].includes(S.layout));$("period").disabled=!(LWN.includes(S.layout)||LSP.includes(S.layout))||S.layout==="wave210";}')
rep('function tipW(n){const L=S.layout,P=S.P;if(!LWN.includes(L))return "";', 'function tipW(n){const L=S.layout,P=S.P;if(LSN.includes(L))return S.l20&&S.l20.tip?S.l20.tip(n):"";if(!LWN.includes(L))return "";')
rep(' else if(S.layout==="hilbert"){', ' else if(LSN.includes(S.layout)){const o=LS(S.layout,N,S.w,S.P,pos,S.sieve);S.cols=o.cols;S.rows=o.rows;S.dw=o.dw;S.fx=o.fx;S.fy=o.fy;S.l20=o;S.par=o.parity||null;}\n else if(S.layout==="hilbert"){')
rep(' $("why").textContent=W[S.layout];', ' Object.assign(W,lsNotes());$("why").textContent=W[S.layout];')
rep('<input id="period" type="number" min="2" max="10000" value="12">', '<input id="period" type="number" min="1" max="10000" value="12">')
rep('let v=Math.round(Number(e.target.value))||12;v=Math.max(2,Math.min(10000,v));', 'let v=Math.round(Number(e.target.value))||12;v=Math.max(S.layout==="quadform"?1:2,Math.min(10000,v));')
rep('const pp=Number(q.get("period"));if(pp>=2&&pp<=10000)S.P=Math.round(pp);', 'const pp=Number(q.get("period"));if(pp>=(L==="quadform"?1:2)&&pp<=10000)S.P=Math.round(pp);')
rep('function setCtl(){', 'function setCtl(){if(S.P<2&&S.layout!=="quadform"){S.P=2;$("period").value=2;}')
open(dst, 'w', encoding='utf-8').write(t)
print('ok', len(t))
