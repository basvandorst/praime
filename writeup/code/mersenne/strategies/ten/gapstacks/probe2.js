(()=>{const cfg=[[12,40,'rand',0],[12,40,'mod',0],[30,60,'rand',0],[7,40,'mod',0],[12,40,'low',300000]],res=[];
for(const [W,H,rule,lim] of cfg){$('bw').value=W;$('bh').value=H;$('rule').value=rule;reset();paused=true;let guard=0;
 while(!S.over&&(!lim||S.n<lim)&&guard<2000000){const b=nextBar();if(!b)break;place(b.i,b.c);guard++;}
 const o=S.over;res.push(o?{W,H,rule,placed:S.n,lines:S.lines,over_gap:o.g,over_after:o.end?null:PR[o.i],col:o.c}:{W,H,rule,placed:S.n,lines:S.lines,end:true});}
draw();return JSON.stringify(res);})()
