"""prime2.html: the gaps between Mersenne exponents as numbers on the links, and a toggle that highlights the exponents
(19 November 2026). Appends replacements to gen.py, which applies them to base.html at every build."""
g = open('gen.py', encoding='utf-8').read()
assert 'edit33' not in g
add = r'''
# edit33: gaps on the links, highlight toggle
rep('<option value="chain">In order (2 → 3 → 5 → …)</option>','<option value="chain">In order (2 → 3 → 5 → …)</option><option value="gaps">In order, with the gaps as numbers</option>')
rep('<label class="check"><input id="labels" type="checkbox" checked> Numbers at the pink dots</label>','<label class="check"><input id="labels" type="checkbox" checked> Numbers at the pink dots</label>\n    <label class="check"><input id="hiexp" type="checkbox"> Highlight the exponents</label>')
rep('$("labels").addEventListener("change",e=>{S.labels=e.target.checked;draw();});','$("labels").addEventListener("change",e=>{S.labels=e.target.checked;draw();});\n$("hiexp").addEventListener("change",e=>{S.hiexp=e.target.checked;draw();});')
rep('if(k&&["all","chain","none"].includes(k)){$("links").value=k;S.links=k;}','if(k&&["all","chain","gaps","none"].includes(k)){$("links").value=k;S.links=k;}if(q.get("hi")==="1"){$("hiexp").checked=true;S.hiexp=true;}')
rep('const E=S.vis||S.exps,lw=Math.max(1,0.9*S.dpr);','const E=S.vis||S.exps,lw=Math.max(1,0.9*S.dpr);if(S.hiexp&&E.length){ctx.setTransform(1,0,0,1,0,0);ctx.fillStyle="rgba(14,24,40,0.62)";ctx.fillRect(0,0,cv.width,cv.height);}')
rep('lk=S.links==="all"?pairs:S.links==="chain"?Math.max(0,m-1):0;','lk=S.links==="all"?pairs:(S.links==="chain"||S.links==="gaps")?Math.max(0,m-1):0;')
rep('(S.links==="all"?` (every pair: ${m}·${m-1}/2)`:S.links==="chain"?" (each to the next)":"")','(S.links==="all"?` (every pair: ${m}·${m-1}/2)`:S.links==="chain"?" (each to the next)":S.links==="gaps"?" (each to the next; the number on a link is the gap)":"")+(S.links==="gaps"&&m>1?`. The gaps: ${E.slice(1).map((k,i)=>(k-E[i]).toLocaleString("en-US")).join(", ")}`:"")')
rep('const R=Math.max(3.2*S.dpr,sd*0.42);','const R=S.hiexp?Math.max(6.5*S.dpr,sd*0.6):Math.max(3.2*S.dpr,sd*0.42);\n if(S.links==="gaps"&&E.length>1){ctx.save();ctx.font=`600 ${Math.round(11*S.dpr)}px system-ui,sans-serif`;ctx.textAlign="center";ctx.textBaseline="middle";\n  for(let i=1;i<E.length;i++){const [ax,ay]=center(E[i-1]),[bx,by]=center(E[i]),t=(E[i]-E[i-1]).toLocaleString("en-US"),mx=(ax+bx)/2,my=(ay+by)/2,tw=ctx.measureText(t).width+8*S.dpr,th=15*S.dpr;\n   ctx.fillStyle="rgba(11,18,32,0.85)";ctx.strokeStyle="rgba(250,204,21,0.9)";ctx.lineWidth=1*S.dpr;ctx.beginPath();ctx.rect(mx-tw/2,my-th/2,tw,th);ctx.fill();ctx.stroke();ctx.fillStyle="#fde68a";ctx.fillText(t,mx,my+0.5*S.dpr);}ctx.restore();}')
rep('ctx.fillStyle="#f472b6";ctx.strokeStyle="#fde7f3";ctx.lineWidth=1.2*S.dpr;ctx.beginPath();ctx.arc(px,py,R,0,6.2832);ctx.fill();ctx.stroke();','if(S.hiexp){ctx.strokeStyle="rgba(250,204,21,0.95)";ctx.lineWidth=2.4*S.dpr;ctx.beginPath();ctx.arc(px,py,R+3.5*S.dpr,0,6.2832);ctx.stroke();}ctx.fillStyle="#f472b6";ctx.strokeStyle="#fde7f3";ctx.lineWidth=1.2*S.dpr;ctx.beginPath();ctx.arc(px,py,R,0,6.2832);ctx.fill();ctx.stroke();')
'''
a = "assert '//' not in lx.replace('https://','') and '/*' not in lx"
assert g.count(a) == 1
g = g.replace(a, add.strip('\n') + '\n' + a)
open('gen.py', 'w', encoding='utf-8').write(g)
t = open('test.js', encoding='utf-8').read()
a = "+ (r.twins ? '&twins=1' : '');"
assert t.count(a) == 1
open('test.js', 'w', encoding='utf-8').write(t.replace(a, "+ (r.twins ? '&twins=1' : '') + (r.hi ? '&hi=1' : '');"))
print('ok')
