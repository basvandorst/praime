// For every layout of prime2.html: the exact dot positions (as the page computes them) and two screenshots.
// node extract.js [layout ...]   -> data/<layout>.json (+ .pos.bin uint16, .fx.bin/.fy.bin float32), shots/<layout>_full.jpg, _zoom.jpg
const { spawn } = require('child_process'); const fs = require('fs'), path = require('path');
const file = path.resolve('prime2.html');
const B = fs.readdirSync(process.env.HOME + '/.cache/ms-playwright').filter(d => d.startsWith('chromium_headless_shell')).map(d => `${process.env.HOME}/.cache/ms-playwright/${d}/chrome-headless-shell-linux64/chrome-headless-shell`).pop();
const ALL = ['rows','spiral','tree','klauber','hex','sacks','hilbert','snake','diag','shells','morton','peano','sunflower','polar','wheel','binring','doubling','ulam41','floyd','legendre','crt','zigzag','moore','gray','bitrev','btree','decrings','ternrings','theodorus','omega','sieverows','honeycomb','diamond','turtle','mertens','race4','liouville','bintri','binblocks','ruler','binsquares','bands512','sierpinski'];
const CENTER = new Set(['spiral','ulam41','hex','sacks','sunflower','polar','wheel','theodorus','diamond','binring','decrings','ternrings','doubling','binsquares']);
const list = process.argv.length > 2 ? process.argv.slice(2) : ALL;
fs.mkdirSync('data', { recursive: true }); fs.mkdirSync('shots', { recursive: true });
const port = 9300 + Math.floor(Math.random() * 500);
const ch = spawn(B, ['--no-sandbox', '--disable-gpu', '--hide-scrollbars', `--remote-debugging-port=${port}`, '--window-size=1100,1300', 'about:blank'], { stdio: 'ignore' });
const sleep = ms => new Promise(r => setTimeout(r, ms));
(async () => {
  let L; for (let i = 0; i < 50; i++) { try { L = await (await fetch(`http://127.0.0.1:${port}/json/list`)).json(); break; } catch (e) { await sleep(200); } }
  const ws = new WebSocket(L.find(t => t.type === 'page').webSocketDebuggerUrl); let id = 0; const pend = new Map();
  ws.onmessage = m => { const d = JSON.parse(m.data); if (d.id && pend.has(d.id)) { pend.get(d.id)(d); pend.delete(d.id); } };
  await new Promise(r => ws.onopen = r);
  const send = (method, params = {}) => new Promise(r => { const i = ++id; pend.set(i, r); ws.send(JSON.stringify({ id: i, method, params })); });
  const ev = async e => { const r = await send('Runtime.evaluate', { expression: e, returnByValue: true, awaitPromise: true }); if (r.result?.exceptionDetails) throw new Error(JSON.stringify(r.result.exceptionDetails).slice(0, 300)); return r.result?.result?.value; };
  await send('Page.enable');
  for (const lay of list) {
    const t0 = Date.now();
    await send('Page.navigate', { url: 'file://' + file + '?layout=' + lay });
    for (let i = 0; i < 100; i++) { await sleep(300); const ok = await ev(`typeof S!=="undefined"&&S.layout===${JSON.stringify(lay)}&&!!S.pos&&!!S.grid`).catch(() => false); if (ok) break; }
    const meta = JSON.parse(await ev(`(()=>{const b64=a=>{const u=new Uint8Array(a.buffer,a.byteOffset,a.byteLength);let s="";for(let i=0;i<u.length;i+=32768)s+=String.fromCharCode.apply(null,u.subarray(i,i+32768));return btoa(s);};
      let mx=0;for(let i=2;i<S.pos.length;i++)if(S.pos[i]>mx)mx=S.pos[i];
      const p16=mx<65536?Uint16Array.from(S.pos):null;
      return JSON.stringify({layout:S.layout,N:S.N,cols:S.cols,rows:S.rows,dw:S.dw,w:S.w,maxpos:mx,hidden:Array.from(S.pos).filter((v,i)=>i>=2&&i%2===0&&v<0).length,
        pos:p16?b64(p16):b64(Int32Array.from(S.pos)),pos16:!!p16,fx:S.fx?b64(Float32Array.from(S.fx)):null,fy:S.fy?b64(Float32Array.from(S.fy)):null,
        extra:S.l20?{ab:S.l20.ab||null,rowsCount:S.l20.rowsCount||null,walk:S.l20.walk||null}:null});})()`));
    fs.writeFileSync(`data/${lay}.pos.bin`, Buffer.from(meta.pos, 'base64'));
    if (meta.fx) { fs.writeFileSync(`data/${lay}.fx.bin`, Buffer.from(meta.fx, 'base64')); fs.writeFileSync(`data/${lay}.fy.bin`, Buffer.from(meta.fy, 'base64')); }
    delete meta.pos; delete meta.fx; delete meta.fy; meta.float = fs.existsSync(`data/${lay}.fx.bin`);
    fs.writeFileSync(`data/${lay}.json`, JSON.stringify(meta));
    const clip = JSON.parse(await ev('(()=>{const b=document.getElementById("stage").getBoundingClientRect();return JSON.stringify({x:b.x,y:b.y+window.scrollY,width:b.width,height:b.height,scale:1});})()'));
    let s = await send('Page.captureScreenshot', { format: 'jpeg', quality: 72, captureBeyondViewport: true, clip });
    fs.writeFileSync(`shots/${lay}_full.jpg`, Buffer.from(s.result.data, 'base64'));
    await ev(CENTER.has(lay) ? `(()=>{const v=S.view,f=4,cx=cv.width/2,cy=cv.height/2;v.x=cx-(cx-v.x)*f;v.y=cy-(cy-v.y)*f;v.s*=f;draw();})()` : `(()=>{const v=S.view,f=4;v.x=v.x*f;v.y=v.y*f;v.s*=f;draw();})()`);
    s = await send('Page.captureScreenshot', { format: 'jpeg', quality: 72, captureBeyondViewport: true, clip });
    fs.writeFileSync(`shots/${lay}_zoom.jpg`, Buffer.from(s.result.data, 'base64'));
    console.log(`${lay}: ${meta.cols}x${meta.rows} dw ${meta.dw} float ${meta.float} hidden ${meta.hidden} (${((Date.now() - t0) / 1000).toFixed(1)} s)`);
  }
  ws.close(); ch.kill();
})().catch(e => { console.error(e); ch.kill(); process.exit(1); });
