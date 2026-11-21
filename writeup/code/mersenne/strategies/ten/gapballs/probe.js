(()=>{const t0=performance.now();let steps=0;while(!S.full&&steps<400000){step();steps++;if(steps%200===0&&S.n>5&&topOfPile()<0)S.full=true;}
for(let i=0;i<3000;i++)step();
let maxOv=0,out=0;for(let i=0;i<S.n;i++){if(S.x[i]<S.r[i]-1e-6||S.x[i]>S.W-S.r[i]+1e-6||S.y[i]>100-S.r[i]+1e-6)out++;for(let j=i+1;j<S.n;j++){const d=Math.hypot(S.x[i]-S.x[j],S.y[i]-S.y[j]),o=(S.r[i]+S.r[j]-d)/Math.min(S.r[i],S.r[j]);if(o>maxOv)maxOv=o;}}
draw();return JSON.stringify({W:S.W,n:S.n,full:S.full,blocked:S.blocked,steps,ms:Math.round(performance.now()-t0),maxRelOverlap:+maxOv.toFixed(3),out,gaps30:GAPS.slice(0,30).join(","),primes:PR.length});})()
