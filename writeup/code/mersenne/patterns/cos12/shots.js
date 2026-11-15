const { spawn } = require('child_process'); const fs = require('fs'), path = require('path');
const file = path.resolve('prime2.html');
const B = fs.readdirSync(process.env.HOME + '/.cache/ms-playwright').filter(d => d.startsWith('chromium_headless_shell')).map(d => `${process.env.HOME}/.cache/ms-playwright/${d}/chrome-headless-shell-linux64/chrome-headless-shell`).pop();
const port = 9300 + Math.floor(Math.random() * 500);
const ch = spawn(B, ['--no-sandbox', '--disable-gpu', '--hide-scrollbars', `--remote-debugging-port=${port}`, '--window-size=900,2400', 'about:blank'], { stdio: 'ignore' });
const sleep = ms => new Promise(r => setTimeout(r, ms));
const JOBS = JSON.parse(process.argv[2]);
(async () => {
  let L; for (let i = 0; i < 50; i++) { try { L = await (await fetch(`http://127.0.0.1:${port}/json/list`)).json(); break; } catch (e) { await sleep(200); } }
  const ws = new WebSocket(L.find(t => t.type === 'page').webSocketDebuggerUrl); let id = 0; const pend = new Map();
  ws.onmessage = m => { const d = JSON.parse(m.data); if (d.id && pend.has(d.id)) { pend.get(d.id)(d); pend.delete(d.id); } };
  await new Promise(r => ws.onopen = r);
  const send = (method, params = {}) => new Promise(r => { const i = ++id; pend.set(i, r); ws.send(JSON.stringify({ id: i, method, params })); });
  const ev = async e => { const r = await send('Runtime.evaluate', { expression: e, returnByValue: true }); return r.result?.result?.value; };
  await send('Page.enable');
  for (const j of JOBS) {
    await send('Page.navigate', { url: 'file://' + file + '?' + j.q });
    for (let i = 0; i < 100; i++) { await sleep(300); if (await ev('typeof S!=="undefined"&&!!S.grid')) break; }
    await ev(`(()=>{const z=document.querySelector(".zoom");if(z)z.style.display="none";${j.js || ''}})()`);
    if (j.zoom) await ev(`(()=>{const v=S.view,f=${j.zoom};v.x=v.x*f;v.y=v.y*f;v.s*=f;draw();})()`);
    const clip = JSON.parse(await ev('(()=>{const b=document.getElementById("stage").getBoundingClientRect();return JSON.stringify({x:b.x,y:b.y+window.scrollY,width:b.width,height:b.height,scale:1});})()'));
    const s = await send('Page.captureScreenshot', { format: 'jpeg', quality: 70, captureBeyondViewport: false, clip });
    fs.writeFileSync(`${j.name}.png`, Buffer.from(s.result.data, 'base64'));
    console.log(j.name, await ev('JSON.stringify({w:S.w,N:S.N,layout:S.layout,P:S.P,primes:S.showPrimes,twins:S.twins})'));
  }
  ws.close(); ch.kill();
})().catch(e => { console.error(e); ch.kill(); process.exit(1); });
