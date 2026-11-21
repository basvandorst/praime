const { spawn } = require('child_process'); const fs = require('fs'), path = require('path');
const file = path.resolve(process.argv[2]); const runs = JSON.parse(process.argv[3]);
const B = fs.readdirSync(process.env.HOME + '/.cache/ms-playwright').filter(d => d.startsWith('chromium_headless_shell')).map(d => `${process.env.HOME}/.cache/ms-playwright/${d}/chrome-headless-shell-linux64/chrome-headless-shell`).pop();
const port = 9400 + Math.floor(Math.random() * 400);
const ch = spawn(B, ['--no-sandbox', '--disable-gpu', '--hide-scrollbars', `--remote-debugging-port=${port}`, '--window-size=' + (process.env.WW || 1100) + ',' + (process.env.WH || 1300), 'about:blank'], { stdio: 'ignore' });
const sleep = ms => new Promise(r => setTimeout(r, ms));
(async () => {
  let L; for (let i = 0; i < 50; i++) { try { L = await (await fetch(`http://127.0.0.1:${port}/json/list`)).json(); break; } catch (e) { await sleep(200); } }
  const ws = new WebSocket(L.find(t => t.type === 'page').webSocketDebuggerUrl); let id = 0; const pend = new Map(); const errs = [];
  ws.onmessage = m => { const d = JSON.parse(m.data); if (d.id && pend.has(d.id)) { pend.get(d.id)(d); pend.delete(d.id); } if (d.method === 'Runtime.exceptionThrown') errs.push(JSON.stringify(d.params.exceptionDetails).slice(0, 300)); };
  await new Promise(r => ws.onopen = r);
  const send = (method, params = {}) => new Promise(r => { const i = ++id; pend.set(i, r); ws.send(JSON.stringify({ id: i, method, params })); });
  const ev = async e => { const r = await send('Runtime.evaluate', { expression: e, returnByValue: true }); return r.result?.result?.value; };
  await send('Page.enable'); await send('Runtime.enable');
  for (const r of runs) {
    errs.length = 0;
    await send('Page.navigate', { url: 'file://' + file + '?' + (r.q || '') }); await sleep(900);
    for (const d of (r.drags || [])) {
      const b = JSON.parse(await ev(`(()=>{const b=cv.getBoundingClientRect();return JSON.stringify({x:b.left+S.x[${d.n}]/S.dpr,y:b.top+S.y[${d.n}]/S.dpr,w:b.width,h:b.height,l:b.left,t:b.top});})()`));
      const tx = b.l + d.fx * b.w, ty = b.t + d.fy * b.h;
      await send('Input.dispatchMouseEvent', { type: 'mousePressed', x: b.x, y: b.y, button: 'left', clickCount: 1 });
      const steps = 25; for (let i = 1; i <= steps; i++) { await send('Input.dispatchMouseEvent', { type: 'mouseMoved', x: b.x + (tx - b.x) * i / steps, y: b.y + (ty - b.y) * i / steps, button: 'left', buttons: 1 }); await sleep(25); }
      await sleep(250); await send('Input.dispatchMouseEvent', { type: 'mouseReleased', x: tx, y: ty, button: 'left', clickCount: 1 }); await sleep(d.wait || 400);
    }
    if (r.js) { await ev(r.js); }
    await sleep(r.wait || 1500);
    const info = JSON.parse(await ev('JSON.stringify({m:[...document.querySelectorAll(".status")].map(e=>e.innerText).join(" || "),sw:document.documentElement.scrollWidth+"/"+document.documentElement.clientWidth})')); info.o = r.out ? await ev(r.out) : '';
    const clip = JSON.parse(await ev(`(()=>{const e=document.getElementById("${r.clip||""}");if(!e)return JSON.stringify({x:0,y:0,width:document.documentElement.scrollWidth,height:document.documentElement.scrollHeight,scale:1});const b=e.getBoundingClientRect();return JSON.stringify({x:b.x,y:b.y+window.scrollY,width:b.width,height:b.height,scale:1});})()`));
    const s = await send('Page.captureScreenshot', { format: 'jpeg', quality: 75, captureBeyondViewport: true, clip });
    fs.writeFileSync(`shots/${r.tag}.jpg`, Buffer.from(s.result.data, 'base64'));
    console.log(r.tag, '| sw', info.sw, '| errors', errs.length, errs[0] || '', '\n  ', info.m, '\n   out:', info.o);
  }
  ws.close(); ch.kill();
})().catch(e => { console.error(e); ch.kill(); process.exit(1); });
