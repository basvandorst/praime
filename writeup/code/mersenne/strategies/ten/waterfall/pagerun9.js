const { spawn } = require('child_process'); const fs = require('fs'), path = require('path');
const file = path.resolve(process.argv[2]), outp = process.argv[3];
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
  const js = `(()=>{const out={};for(const X of [1000,5000,100000,1000000]){const s=seqOf("primes",X,0),r=analyse(s);out["primes_"+X]={depth:r.depth,holds:r.holds,firsts:r.firsts.slice(0,12),n:s.length};}const m=seqOf("mexp",1000000,0),r=analyse(m);out.mexp_1000000={depth:r.depth,holds:r.holds,firsts:r.firsts.slice(0,12),n:m.length};return JSON.stringify(out);})()`;
  const r = await send('Runtime.evaluate', { expression: js, returnByValue: true });
  fs.writeFileSync(outp, r.result.result.value); console.log('page results written, script errors', errs.length);
  ws.close(); ch.kill();
})().catch(e => { console.error(e); ch.kill(); process.exit(1); });
