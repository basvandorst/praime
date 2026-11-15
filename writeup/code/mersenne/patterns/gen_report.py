"""Build report/prime2-patterns.html from lines_out/*.json (grid layouts) and texts.json (all layouts: seen, explained,
exponents for the round layouts and walks). Images: prime2-patterns/<layout>.jpg and _z.jpg."""
import json, os, html
ORDER = ['rows','snake','honeycomb','crt','bands512','binblocks','spiral','ulam41','hex','diamond','shells','klauber','legendre','tree','floyd','diag','zigzag','sierpinski',
         'hilbert','moore','peano','morton','gray','bitrev','sacks','sunflower','theodorus','polar','wheel','binring','decrings','ternrings','binsquares','doubling','bintri','btree',
         'omega','sieverows','ruler','turtle','mertens','race4','liouville']
NAMES = {'rows':'Rows (width 708)','spiral':'Ulam spiral','tree':'Tree (1, 2, 3, … per row)','klauber':'Klauber triangle','hex':'Hexagonal spiral','sacks':'Sacks spiral','hilbert':'Hilbert curve',
 'snake':'Snake rows','diag':'Diagonals (Cantor)','shells':'Square shells','morton':'Z-order (Morton)','peano':'Peano curve','sunflower':'Sunflower','polar':'Polar: angle n, distance n','wheel':'Wheel (708 spokes)',
 'binring':'Binary rings','doubling':'Doubling spiral','ulam41':'Ulam spiral from 41','floyd':"Floyd's triangle",'legendre':'Between squares (Legendre)','crt':'Chinese remainder grid','zigzag':'Zigzag diagonals',
 'moore':'Moore curve','gray':'Gray-code order','bitrev':'Bit-reversed order','btree':'Binary tree (heap)','decrings':'Decimal rings','ternrings':'Ternary rings','theodorus':'Spiral of Theodorus',
 'omega':'Rows by number of prime factors','sieverows':'Sieve of Eratosthenes (rows)','honeycomb':'Honeycomb rows','diamond':'Diamond spiral','turtle':'Prime turtle','mertens':'Mertens walk',
 'race4':'Prime race mod 4','liouville':'Liouville walk','bintri':'Binary triangle','binblocks':'Binary blocks','ruler':'Powers-of-2 pyramid','binsquares':'Binary squares','bands512':'Bands of 512',
 'sierpinski':"Pascal's triangle mod 2"}
T = json.load(open('texts.json', encoding='utf-8'))

def canon(r):
    A = r['d2'] / 2; B = r['d1'] - r['d2'] / 2; C = r['n0']
    if A == 0: return f"{int(B)}t + {C}"
    k = round(B / (2 * A)); B2 = B - 2 * A * k; C2 = A * k * k - B * k + C
    f = lambda x: (f"{x:g}" if x != int(x) else f"{int(x)}")
    s = (f"{f(A)}s²" if A != 1 else "s²") + (f" + {f(B2)}s" if B2 > 0 else (f" − {f(-B2)}s" if B2 < 0 else "")) + (f" + {f(C2)}" if C2 > 0 else (f" − {f(-C2)}" if C2 < 0 else ""))
    return s

def grid_text(lay):
    f = f"lines_out/{'tree' if lay == 'sierpinski' else lay}.json"
    if not os.path.exists(f): return None
    r = json.load(open(f)); e = r['exponents']
    if not r['top']:
        lines = "No straight lines: not one stretch of 30 numbers on a line follows a polynomial, in any of 24 directions."
    else:
        t = r['top'][0]; d = r['bottom'][0]; s = r['sample']
        lines = (f"Brightest line: <b>{html.escape(canon(t))}</b> ({t['m']} numbers): <b>{t['obs']}</b> primes, where random numbers of that size give {t['naive']:.1f} "
                 f"and Bateman–Horn predicts <b>{t['pred']:.1f}</b> (constant {t['C']:.2f}). Darkest: {html.escape(canon(d))}, constant {d['C']:.2f}, {d['obs']} primes. "
                 f"{r['segments']:,} polynomial stretches in all; a random 150 of them differ from the prediction by {s['mean']:+.2f} ± {s['sd']:.2f} standard deviations.")
    ex = (f"Exponents: at most {e['max_on_line']} on one line ({', '.join(f'{k:,}' for k in e['on'])}); random sets with the same growth give {e['null_median']:.0f} "
          f"(95%: {e['null_95']:.0f}), p = {e['p']:.2f}.")
    return lines, ex

cards = []; rows = []
for lay in ORDER:
    t = T.get(lay, {})
    g = grid_text(lay)
    num = g[0] if g else t.get('numbers', '')
    ex = g[1] if g else t.get('exponents', '')
    cards.append(f'''<section class="card" id="{lay}"><h2>{html.escape(NAMES[lay])} <a class="open" href="prime2.html?layout={lay}&amp;links=none">open</a></h2>
<div class="imgs"><img src="prime2-patterns/{lay}.jpg" alt="{html.escape(NAMES[lay])}, whole view" loading="lazy"><img src="prime2-patterns/{lay}_z.jpg" alt="{html.escape(NAMES[lay])}, zoomed four times" loading="lazy"></div>
<p><b>What you see.</b> {t.get('seen', '')}</p><p><b>What it is.</b> {t.get('explained', '')}</p><p class="num">{num}</p><p class="ex">{ex}</p></section>''')
    rows.append(f'<tr><td><a href="#{lay}">{html.escape(NAMES[lay])}</a></td><td>{t.get("short", "")}</td><td>{t.get("verdict", "explained")}</td></tr>')

page = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Prime Layout Patterns</title>
<style>
:root {{ color-scheme: dark; --bg:#0b1220; --panel:#121e30; --muted:#9fb0c6; --edge:#26364c; --pink:#f472b6; --accent:#38bdf8; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:#edf4ff; font:15px/1.5 system-ui,-apple-system,Segoe UI,sans-serif; }}
header {{ padding:26px 16px 8px; text-align:center; }}
h1 {{ margin:0 0 8px; font-size:clamp(22px,4vw,30px); }}
header p {{ margin:6px auto; max-width:820px; color:var(--muted); }}
main {{ max-width:1100px; margin:auto; padding:10px 16px 40px; }}
.panel, .card {{ background:var(--panel); border:1px solid var(--edge); border-radius:12px; padding:14px 16px; margin:14px 0; }}
table {{ width:100%; border-collapse:collapse; font-size:14px; }}
td, th {{ text-align:left; padding:6px 8px; border-top:1px solid #1f2e44; vertical-align:top; }}
th {{ color:#c9d8ec; font-weight:600; }}
a {{ color:var(--accent); }}
.card h2 {{ margin:0 0 10px; font-size:18px; display:flex; justify-content:space-between; gap:10px; }}
.card h2 .open {{ font-size:13px; font-weight:400; }}
.imgs {{ display:grid; grid-template-columns:1fr 1fr; gap:8px; }}
.imgs img {{ width:100%; border-radius:8px; border:1px solid var(--edge); background:#0e1828; }}
.card p {{ margin:8px 0 0; color:#d6e2f2; }}
.card p.num {{ color:#fde68a; }}
.card p.ex {{ color:#fbcfe8; }}
.wrap {{ overflow-x:auto; }}
@media (max-width:640px) {{ .imgs {{ grid-template-columns:1fr; }} }}
</style>
</head>
<body>
<header>
<h1>Patterns in the prime layouts</h1>
<p>Every layout of <a href="prime2.html">prime2.html</a>, photographed whole and zoomed in four times, with the exponent links off so that the pattern of the primes shows; on the layouts with straight lines the brightest line found is drawn in yellow by the page itself.</p>
<p>{T['_intro']['method']}</p>
</header>
<main>
<section class="panel"><h2 style="margin:0 0 8px;font-size:18px">The answer in short</h2><p style="margin:0 0 10px;color:#d6e2f2">{T['_intro']['answer']}</p>
<div class="wrap"><table><thead><tr><th>Layout</th><th>The pattern</th><th>Verdict</th></tr></thead><tbody>{''.join(rows)}</tbody></table></div></section>
{''.join(cards)}
</main>
</body>
</html>
'''
assert '<!--' not in page and '/*' not in page
open('report/prime2-patterns.html', 'w', encoding='utf-8').write(page)
print('ok', len(page))
