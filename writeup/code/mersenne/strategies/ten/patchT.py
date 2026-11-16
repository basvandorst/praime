"""Insert the ten layouts (lx.js, opts.txt) into the 70-layout page.
Usage: put that page next to this file as base.html, then python3 patchT.py -> prime2.html.
Test with ../page/test.js."""
import re,sys
b=open('base.html').read(); lx=open('lx.js').read(); opts=open('opts.txt').read().strip().split('\n')
def rep(old,new,cnt=1):
    global b
    assert b.count(old)==cnt,(old[:60],b.count(old)); b=b.replace(old,new)
o=''.join('<option value="%s">%s</option>'%tuple(l.split('|',1)) for l in opts if l)
rep('<option value="llres">Lucas–Lehmer residues (100·P primes)</option>','<option value="llres">Lucas–Lehmer residues (100·P primes)</option>'+o)
rep('const LS_IDO=[',lx.strip()+'\nconst LS_IDO=[')
rep('Object.assign(W,lsNotes());','Object.assign(W,lsNotes(),lxNotes());')
rep(' else if(S.layout==="hilbert"){',' else if(LXN.includes(S.layout)){const o=LX(S.layout,N,S.w,S.P,pos,S.sieve);S.cols=o.cols;S.rows=o.rows;S.dw=o.dw;S.fx=o.fx;S.fy=o.fy;S.l20=o;S.par=o.parity||null;}\n else if(S.layout==="hilbert"){')
rep('S.pos=pos;S.vis=S.exps.filter(k=>pos[2*k]>=0);','S.pos=pos;S.vis=S.l20&&S.l20.noexp?[]:S.exps.filter(k=>pos[2*k]>=0);')
rep('...LWR,...LSW].includes(S.layout));$("period")','...LWR,...LSW,...LXW].includes(S.layout));$("period")')
rep('LSP.includes(S.layout))||S.layout==="wave210"','LSP.includes(S.layout)||LXP.includes(S.layout))||S.layout==="wave210"')
rep('if(S.P<2&&S.layout!=="quadform")','if(S.P<2&&S.layout!=="quadform"&&S.layout!=="spectrum")')
rep('if(LSN.includes(L))return','if(LSN.includes(L)||LXN.includes(L))return')
rep('else if(S.par)c=S.par[n]?[111,143,184]:[38,54,76];','else if(S.par)c=S.par[n]===2?[52,211,153]:S.par[n]?[111,143,184]:[38,54,76];')
rep('tw=!S.par&&S.twins&&isTwin(n),pr=S.par?S.par[n]===1:tw||(S.showPrimes&&!S.sieve[n])','tw=S.par?S.par[n]===2:S.twins&&isTwin(n),pr=S.par?S.par[n]>=1:tw||(S.showPrimes&&!S.sieve[n])',2)
rep('ctx.fillText(n,px,py)','ctx.fillText(LB(n),px,py)',2)
rep('setFont(n);','setFont(LB(n));',2)
rep('tip.hidden=false;tip.innerHTML=','tip.hidden=false;tip.innerHTML=S.l20&&S.l20.tipFull?S.l20.tipFull(n):')
rep('if(pp>=(L==="quadform"?1:2)','if(pp>=(L==="quadform"||L==="spectrum"?1:2)')
rep('Math.max(S.layout==="quadform"?1:2,','Math.max(S.layout==="quadform"||S.layout==="spectrum"?1:2,')
rep('sqW(){return Math.ceil(Math.sqrt(S.N*((LWR.includes(S.layout)||LSW.includes(S.layout))?4:1)));}','sqW(){return Math.ceil(Math.sqrt(S.N*((LWR.includes(S.layout)||LSW.includes(S.layout))?4:1)));}')
assert '//' not in lx.replace('https://','') and '/*' not in lx
open('prime2.html','w').write(b); print('built',len(b))
