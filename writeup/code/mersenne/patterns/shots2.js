// Report screenshots: links off; on grid layouts the brightest line of lines_out/<layout>.json drawn by the page's own
// line test. node shots2.js [layout ...] -> report/img/<layout>.jpg (full view) and report/img/<layout>_z.jpg (zoom x4)
const { spawn } = require('child_process'); const fs = require('fs'), path = require('path');
const file = path.resolve('prime2.html');
const B = fs.readdirSync(process.env.HOME + '/.cache/ms-playwright').filter(d => d.startsWith('chromium_headless_shell')).map(d => `${process.env.HOME}/.cache/ms-playwright/${d}/chrome-headless-shell-linux64/chrome-headless-shell`).pop();
const ALL = ['rows','spiral','tree','klauber','hex','sacks','hilbert','snake','diag','shells','morton','peano','sunflower','polar','wheel','binring','doubling','ulam41','floyd','legendre','crt','zigzag','moore','gray','bitrev','btree','decrings','ternrings','theodorus','omega','sieverows','honeycomb','diamond','turtle','mertens','race4','liouville','bintri','binblocks','ruler','binsquares','bands512','sierpinski'];
const CENTER = new Set(['spiral','ulam41','hex','sacks','sunflower','polar','wheel','theodorus','diamond','binring','decrings','ternrings','doubling','binsquares']);
const list = process.argv.length > 2 ? process.argv.slice(2) : ALL;
fs.mkdirSync('report/img', { recursive: true });
const port = 9300 + Math.floor(Math.random() * 500);
const ch = spawn(B, ['--no-sandbox', '--disable-gpu', '--hide-scrollbars', `--remote-debugging-port=${port}`, '--window-size=900,2400', 'about:blank'], { stdio: 'ignore' });
const sleep = ms => new Promise(r => setTimeout(r, ms));
(async () => {
  let L; for (let i = 0; i < 50; i++) { try { L = await (await fetch(`http://127.0.0.1:${port}/json/list`)).json(); break; } catch (e) { await sleep(200); } }
  const ws = new WebSocket(L.find(t => t.type === 'page').webSocketDebuggerUrl); let id = 0; const pend = new Map();
  ws.onmessage = m => { const d = JSON.parse(m.data); if (d.id && pend.has(d.id)) { pend.get(d.id)(d); pend.delete(d.id); } };
  await new Promise(r => ws.onopen = r);
  const send = (method, params = {}) => new Promise(r => { const i = ++id; pend.set(i, r); ws.send(JSON.stringify({ id: i, method, params })); });
  const ev = async e => { const r = await send('Runtime.evaluate', { expression: e, returnByValue: true }); return r.result?.result?.value; };
  await send('Page.enable');
  for (const lay of list) {
    await send('Page.navigate', { url: 'file://' + file + '?layout=' + lay + '&links=none' });
    for (let i = 0; i < 100; i++) { await sleep(300); if (await ev(`typeof S!=="undefined"&&S.layout===${JSON.stringify(lay)}&&!!S.grid`)) break; }
    let a = 0, b = 0;
    const jf = `lines_out/${lay === 'sierpinski' ? 'tree' : lay}.json`;
    if (fs.existsSync(jf)) { const j = JSON.parse(fs.readFileSync(jf, 'utf8')); if (j.top && j.top.length) { a = j.top[0].n0; b = j.top[0].n0 + j.top[0].d1; } }
    if (a && b) await ev(`(()=>{S.test.on=true;testLine(${a},${b});})()`);
    await ev('(()=>{const z=document.querySelector(".zoom");if(z)z.style.display="none";})()');
    const clip = JSON.parse(await ev('(()=>{const b=document.getElementById("stage").getBoundingClientRect();return JSON.stringify({x:b.x,y:b.y+window.scrollY,width:b.width,height:b.height,scale:1});})()'));
    let s = await send('Page.captureScreenshot', { format: 'jpeg', quality: 62, captureBeyondViewport: false, clip });
    fs.writeFileSync(`report/img/${lay}.jpg`, Buffer.from(s.result.data, 'base64'));
    await ev(CENTER.has(lay) ? `(()=>{const v=S.view,f=4,cx=cv.width/2,cy=cv.height/2;v.x=cx-(cx-v.x)*f;v.y=cy-(cy-v.y)*f;v.s*=f;draw();})()` : `(()=>{const v=S.view,f=4;v.x=v.x*f;v.y=v.y*f;v.s*=f;draw();})()`);
    s = await send('Page.captureScreenshot', { format: 'jpeg', quality: 62, captureBeyondViewport: false, clip });
    fs.writeFileSync(`report/img/${lay}_z.jpg`, Buffer.from(s.result.data, 'base64'));
    const txt = await ev('document.getElementById("lineTest").textContent');
    console.log(lay, a ? `line ${a},${b}: ${(txt || '').slice(0, 140)}` : '(no line)');
  }
  ws.close(); ch.kill();
})().catch(e => { console.error(e); ch.kill(); process.exit(1); });
