const { spawn } = require('child_process'); const fs = require('fs'), path = require('path');
const file = path.resolve(process.argv[2]), X = Number(process.argv[3]), outp = process.argv[4];
const B = fs.readdirSync(process.env.HOME + '/.cache/ms-playwright').filter(d => d.startsWith('chromium_headless_shell')).map(d => `${process.env.HOME}/.cache/ms-playwright/${d}/chrome-headless-shell-linux64/chrome-headless-shell`).pop();
const port = 9400 + Math.floor(Math.random() * 400);
const ch = spawn(B, ['--no-sandbox', '--disable-gpu', `--remote-debugging-port=${port}`, 'about:blank'], { stdio: 'ignore' });
const sleep = ms => new Promise(r => setTimeout(r, ms));
(async () => {
  let L; for (let i = 0; i < 50; i++) { try { L = await (await fetch(`http://127.0.0.1:${port}/json/list`)).json(); break; } catch (e) { await sleep(200); } }
  const ws = new WebSocket(L.find(t => t.type === 'page').webSocketDebuggerUrl); let id = 0; const pend = new Map(); const errs = [];
  ws.onmessage = m => { const d = JSON.parse(m.data); if (d.id && pend.has(d.id)) { pend.get(d.id)(d); pend.delete(d.id); } if (d.method === 'Runtime.exceptionThrown') errs.push(1); };
  await new Promise(r => ws.onopen = r);
  const send = (method, params = {}) => new Promise(r => { const i = ++id; pend.set(i, r); ws.send(JSON.stringify({ id: i, method, params })); });
  await send('Page.enable'); await send('Runtime.enable');
  await send('Page.navigate', { url: 'file://' + file }); await sleep(1500);
  const js = `(()=>{const R=["alt","light","mod4","mod3","mod8","lastdig","repeat","binsum","decsum","twincous","twinsexy","sophie","gap24","gap26","updown","liou","mob","omega"],M=["alt","light","mod4","mod3","mod8","lastdig","binsum","decsum"],out={};
   for(const [src,rules] of [["primes",R],["mexp",M]])for(const r of rules)for(const wt of ["one","val","log","recip"]){S.rule=r;S.wt=wt;S.src=src;S.maxN=${X};reset();advance(1e9);out[src+"|"+r+"|"+wt]={k:S.k,L:S.L,R:S.R,eq:S.eq,aL:S.aheadL,aR:S.aheadR,ch:S.changes,bL:S.bestL[0],bR:S.bestR[0],fL:S.firstL,fR:S.firstR,eqsteps:S.eqLog.slice(0,200).map(e=>e[0])};}
   return JSON.stringify(out);})()`;
  const r = await send('Runtime.evaluate', { expression: js, returnByValue: true });
  fs.writeFileSync(outp, r.result.result.value); console.log('page runs written, script errors', errs.length);
  ws.close(); ch.kill();
})().catch(e => { console.error(e); ch.kill(); process.exit(1); });
