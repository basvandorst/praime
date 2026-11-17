"""Insert the layouts (lx.js, opts.txt) into the 70-layout page.
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
rep('Object.assign(W,lsNotes());','Object.assign(W,lsNotes(),lxNotes());if(W[S.layout])W[S.layout]+=lxScoreText(S.layout);')
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
rep('sqW(){return Math.ceil(Math.sqrt(S.N*((LWR.includes(S.layout)||LSW.includes(S.layout))?4:1)));}','sqW(){return Math.ceil(Math.sqrt(S.N*((LWR.includes(S.layout)||LSW.includes(S.layout)||S.layout==="parity")?4:1)));}')
rep('.note { color:var(--muted); font-size:12.5px; margin:8px 0 0; }','.note { color:var(--muted); font-size:12.5px; margin:8px 0 0; }\n.score { display:flex; flex-wrap:wrap; gap:8px 10px; align-items:center; font-size:12.5px; }\n.score .t { background:#16243a; border:1px solid #2a3b55; border-radius:8px; padding:3px 10px; color:#c9d8ec; }\n.score .t b { color:#fff; font-size:15px; margin-right:5px; }\n.score .ok { color:#34d399; font-weight:600; }\n.score .warn { color:#fde68a; font-weight:600; }\nh2 { margin:0 0 4px; font-size:17px; }\n.found { margin:8px 0 12px; padding-left:20px; color:#d6e2f2; font-size:13px; }\n.found li { margin:6px 0; }\n.boardwrap { max-height:440px; overflow:auto; border:1px solid var(--edge); border-radius:8px; }\ntable { border-collapse:collapse; width:100%; font-size:12.5px; font-variant-numeric:tabular-nums; }\nth, td { padding:4px 8px; text-align:right; border-bottom:1px solid #1e2d44; }\nth:nth-child(2), td:nth-child(2) { text-align:left; min-width:170px; }\nth { position:sticky; top:0; background:#16243a; color:#b8c8dc; font-weight:600; }\nth[data-c] { cursor:pointer; }\ntr.cur td { background:#1b2a40; }\ntd.hi { color:#fde68a; font-weight:600; }')
rep('  <section class="panel">\n    <div class="stats" id="stats"></div>','  <section class="panel score" id="scoreBox" aria-live="polite"></section>\n  <section class="panel">\n    <div class="stats" id="stats"></div>')
rep('  </section>\n</main>','  </section>\n  <section class="panel" id="results">\n    <h2>What the backtest found</h2>\n    <p class="note">Every layout was scored at N = 500,000 against decoys: numbers that pass the sieve up to N<sup>0.4027</sup> and are then thinned by a coin. Decoys have the density of the primes and everything the small factors do, and nothing else.</p>\n    <ul class="found" id="found"></ul>\n    <div class="boardwrap"><table id="board"></table></div>\n    <p class="note">Pattern strength: how much a dot being prime tells about its 24 nearest neighbours, in thousandths of a bit. Decoys: the same for the decoys (mean of 8 sets). Beyond the decoys: how far the pair counts of the primes are from those of the decoys, in units of the scatter among the decoys; a decoy scored as if it were the primes gives 0.50 on average, 0.73 at the 95% point and 0.85 at most. Improved model: the products of two primes in simple ratios are struck before the coin is thrown. Select a column head to sort, a layout to open it. More in <a href="prime2-patterns.html#ten-found">the pattern report</a>.</p>\n  </section>\n</main>')
rep('(3.1 expected by chance).</p>','(3.1 expected by chance). Every layout is scored against decoys: the score of the layout on screen is under the picture, <a href="#results">all results</a> are at the bottom.</p>')
rep('$("why").textContent=W[S.layout];','$("why").textContent=W[S.layout];lxStrip();')
rep('fromURL();\n','fromURL();lxBoard();\n')
rep('with screenshots of all layouts: <a href="prime2-patterns.html">the pattern report</a>.</p>','with screenshots of all layouts: <a href="prime2-patterns.html">the pattern report</a>. The first 30 numbers as turning gears: <a href="prime3.html">prime3.html</a>.</p>')
assert '//' not in lx.replace('https://','') and '/*' not in lx
open('prime2.html','w').write(b); print('built',len(b))
