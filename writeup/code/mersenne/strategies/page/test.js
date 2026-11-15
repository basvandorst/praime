const { spawn } = require('child_process'); const fs = require('fs'), path = require('path');
const file = path.resolve(process.argv[2] || 'page/prime2.html');
const runs = JSON.parse(process.argv[3] || '[]');
const B = fs.readdirSync(process.env.HOME + '/.cache/ms-playwright').filter(d => d.startsWith('chromium_headless_shell')).map(d => `${process.env.HOME}/.cache/ms-playwright/${d}/chrome-headless-shell-linux64/chrome-headless-shell`).pop();
const port = 9400 + Math.floor(Math.random() * 400);
const ch = spawn(B, ['--no-sandbox', '--disable-gpu', '--hide-scrollbars', `--remote-debugging-port=${port}`, '--window-size=900,2400', 'about:blank'], { stdio: 'ignore' });
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
    const tag = r.tag || r.layout + (r.period ? '_P' + r.period : '') + (r.width ? '_w' + r.width : '');
    const q = `layout=${r.layout}&links=${r.links || 'none'}` + (r.period ? `&period=${r.period}` : '') + (r.width ? `&width=${r.width}` : '') + (r.n ? `&n=${r.n}` : '') + (r.twins ? '&twins=1' : '');
    errs.length = 0; const t0 = Date.now();
    await send('Page.navigate', { url: 'file://' + file + '?' + q });
    let ok = false; for (let i = 0; i < 200; i++) { await sleep(150); if (await ev(`typeof S!=="undefined"&&S.layout===${JSON.stringify(r.layout)}&&!!S.grid&&!!$("why").textContent`)) { ok = true; break; } }
    const dt = (Date.now() - t0) / 1000;
    const info = JSON.parse(await ev(`JSON.stringify({why:$("why").textContent,cols:S.cols,rows:S.rows,hidden:(()=>{let h=0;for(let n=1;n<=S.N;n++)if(S.pos[2*n]<0)h++;return h;})(),vis:S.vis.length,P:S.P,w:S.w})`) || '{}');
    await ev('(()=>{const z=document.querySelector(".zoom");if(z)z.style.display="none";})()');
    const clip = JSON.parse(await ev('(()=>{const b=document.getElementById("stage").getBoundingClientRect();return JSON.stringify({x:b.x,y:b.y+window.scrollY,width:b.width,height:b.height,scale:1});})()'));
    let s = await send('Page.captureScreenshot', { format: 'jpeg', quality: 70, captureBeyondViewport: false, clip });
    fs.writeFileSync(`shots/${tag}.jpg`, Buffer.from(s.result.data, 'base64'));
    if (r.zoom) {
      await ev(`(()=>{const v=S.view,f=${r.zoom},cx=cv.width*${r.zx ?? 0.5},cy=cv.height*${r.zy ?? 0.5};v.x=cx-(cx-v.x)*f;v.y=cy-(cy-v.y)*f;v.s*=f;draw();})()`);
      s = await send('Page.captureScreenshot', { format: 'jpeg', quality: 70, captureBeyondViewport: false, clip });
      fs.writeFileSync(`shots/${tag}_z.jpg`, Buffer.from(s.result.data, 'base64'));
    }
    fs.writeFileSync(`shots/${tag}.txt`, JSON.stringify({ q, ok, dt, errs: [...errs], ...info }, null, 1));
    console.log(`${tag}: ok ${ok} ${dt.toFixed(1)} s, ${info.cols}x${info.rows}, hidden ${info.hidden}, exps shown ${info.vis}, errors ${errs.length}${errs.length ? ' ' + errs[0] : ''}`);
  }
  ws.close(); ch.kill();
})().catch(e => { console.error(e); ch.kill(); process.exit(1); });
