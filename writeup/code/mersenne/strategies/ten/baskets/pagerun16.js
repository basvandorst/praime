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
  await send('Page.navigate', { url: 'file://' + file }); await sleep(4000);
  const js = `(()=>{S.src="primes";S.N=${X};prepare();const out={};for(const h of [10,30,100,300,1000,3000,10000]){S.h=h;exact();out[h]={ratio:Math.round(S.ex.ratio*1e6)/1e6,wv:Math.round(S.ex.wv*1e4)/1e4,ms:Math.round(S.ex.ms*1e4)/1e4,cr:Math.round(S.ex.cr*1e4)/1e4};}return JSON.stringify(out);})()`;
  const r = await send('Runtime.evaluate', { expression: js, returnByValue: true });
  fs.writeFileSync(outp, r.result.result.value); console.log('page results written, script errors', errs.length);
  ws.close(); ch.kill();
})().catch(e => { console.error(e); ch.kill(); process.exit(1); });
