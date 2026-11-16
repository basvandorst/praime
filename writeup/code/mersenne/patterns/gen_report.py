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

def big_img(f, alt, cls=''):
    return f'<a href="prime2-patterns/{f}"><img class="{cls}" src="prime2-patterns/{f}" alt="{html.escape(alt)}" loading="lazy"></a>'
BIG = f'''<section class="card" id="cos12big"><h2>The big picture: cosine rows, period 12, zoomed out <a class="open" href="prime2.html?layout=cosrows&amp;period=12&amp;n=100000&amp;width=240&amp;links=none">open</a></h2>
<p><b>Zoomed out naively, there is nothing to see.</b> At 500,000 numbers in rows of 1,416 there is less than one pixel per number and the picture is an even speckle of blue and orange: the two lanes are equally full, and inside the places a prime may sit the primes look random. To see the structure, the row width has to fit the wave. Click any picture for full size.</p>
<div class="imgs one">{big_img('cos12_big_zoom.png', 'Cosine rows, period 12, width 240, first 24 rows at full size')}</div>
<p><b>Close up (the first 24 rows of width 240).</b> Each row holds 20 waves; the faint dots are the composites tracing the wave, the primes sit on the crest (blue, p ≡ ±1 mod 12) or in the trough (orange, p ≡ ±5). Because 240 is a multiple of 12, every column keeps its height from row to row.</p>
<div class="imgs one">{big_img('cos12_big_stripes.png', 'Cosine rows, period 12, width 240, numbers 1 to 100,000, primes coloured by lane')}</div>
<p><b>Zoomed out (1 to 100,000, 417 rows of 240).</b> The waves line up into stripes. Only 64 of the 240 columns can hold primes (the rest share a factor 2, 3 or 5 with 240), and each of them is a single colour: the pattern of blue and orange stripes repeats every 60 columns, four times across. All 64 stripes are equally full, 137 to 164 primes each, mean 149.8, as Dirichlet’s theorem says; inside a stripe the primes come at random.</p>
<div class="imgs one">{big_img('cos12_big_twins.png', 'Cosine rows, period 12, width 240, twin primes')}</div>
<p><b>Twin primes (same view; other primes dimmed).</b> Twins run on 24 double tracks, pairs of columns two apart that can both hold primes: 12 blue and 12 orange, because a twin pair never changes lane. All 1,222 twin pairs from (11, 13) to 100,000 sit on these tracks (the track in column 239 continues in column 1 of the next row); only (5, 7) in the first row is off them, since 5 divides 240.</p>
<div class="imgs">{big_img('cos12_big_diag.png', 'Cosine rows, period 12, width 277, primes')}{big_img('cos12_big_braid.png', 'Cosine rows, period 12, width 100, primes')}</div>
<p><b>Other widths.</b> Left, width 277 (≡ 1 mod 12): each row starts one step later in the wave, so the lanes slide one column per row and the stripes tilt into short diagonals. Right, width 100 (≡ 4 mod 12, 40,000 numbers): each row moves the wave by four places, so the 40 columns that can hold primes change lane in a cycle of three rows (column 1: 1 up, 101 down, 201 a multiple of 3) and the colours mix within each column.</p>
<div class="chart">{open('cos12/race.svg').read()}</div>
<p><b>Zoomed out to 10 million: the lane race.</b> The only difference between the lanes at large scale is a small lead for the down lane: 332,352 primes against 332,225. Its size is predicted. Every prime square p² (p ≥ 5) is ≡ 1 (mod 12), an up-lane remainder, and in the count of primes and prime powers the four classes are even on average; so the primes alone in class 1 fall short by about half the number of squares, ½ Li(√x) ≈ √x / ln x (Chebyshev’s bias, made precise by Rubinstein and Sarnak in 1994). Measured: averaged over ln x from 1,000 to 10 million the lead is 0.986 times ½ Li(√x) (dashed). The zeta zeros make the wobble around it: the up lane is ahead in 102 short stretches (blue), the first at 61,463, the longest from 1,172,807 to 1,179,403 (481 primes), the last ending at 2,233,787; in all after 1,476 of the 664,577 primes (0.22%).</p>
<div class="chart">{open('cos12/rails.svg').read()}</div>
<p class="ex"><b>The Mersenne exponents on the two lanes.</b> From 5 to 136,279,841, 16 are up and 34 down (16 of the 49 to no. 51); the heuristic expects about half. Hover over a dot for its value. As said above, after all the splits tested this is not significant; each new exponent is a fresh test.</p></section>'''

NT = json.load(open('texts_new.json', encoding='utf-8')) if os.path.exists('texts_new.json') else {'_intro': '', 'cards': []}
def new_card(c):
    imgs = ''.join(f'<a href="prime2-patterns/{f}"><img src="prime2-patterns/{f}" alt="{html.escape(c["name"])}" loading="lazy"></a>' for f in c['imgs'])
    return (f'<section class="card" id="new-{c["id"]}"><h2>{html.escape(c["name"])} <a class="open" href="prime2.html?{c["q"]}">open</a></h2>'
            f'<div class="imgs">{imgs}</div><p><b>What you see.</b> {c["seen"]}</p><p><b>What it is.</b> {c["explained"]}</p><p class="num">{c["checked"]}</p>'
            + (f'<p class="ex">{c["exponents"]}</p>' if c.get('exponents') else '') + '</section>')
NEWS = ('' if not NT['cards'] else f'''<section class="panel" id="new"><h2 style="margin:0 0 8px;font-size:18px">New strategies ({len(NT['cards'])} layouts, 16 November)</h2><p style="margin:0;color:#d6e2f2">{NT['_intro']}</p>
<div class="wrap"><table><thead><tr><th>Layout</th><th>The pattern</th><th>Checked</th></tr></thead><tbody>{''.join(f'<tr><td><a href="#new-{c["id"]}">{html.escape(c["name"])}</a></td><td>{c["short"]}</td><td>{c["verdict"]}</td></tr>' for c in NT['cards'])}</tbody></table></div></section>'''
        + ''.join(new_card(c) for c in NT['cards']))

TN = json.load(open('texts_ten.json', encoding='utf-8')) if os.path.exists('texts_ten.json') else None
TEN = '' if not TN else (f'''<section class="panel" id="ten"><h2 style="margin:0 0 8px;font-size:18px">{TN['_title']}</h2><p style="margin:0 0 8px;color:#d6e2f2">{TN['_intro']}</p><p style="margin:0 0 10px;color:#d6e2f2">{TN['_reveal']}</p>
<div class="wrap"><table><thead><tr><th>Layout</th><th>The pattern</th><th>Checked</th></tr></thead><tbody>{''.join(f'<tr><td><a href="#new-{c["id"]}">{html.escape(c["name"])}</a></td><td>{c["short"]}</td><td>{c["verdict"]}</td></tr>' for c in TN['cards'])}</tbody></table></div></section>
<section class="card" id="ten-found"><h2>Found along the way</h2>''' + ''.join(f'<p><b>{f[0]}</b> ({f[1]}). {f[2]}</p>' + (f'<p class="num">{f[3]}</p>' if f[3] else '') for f in TN['found']) + '</section>'
       + ''.join(new_card(c) for c in TN['cards']))
TT = json.load(open('texts_theory.json', encoding='utf-8')) if os.path.exists('texts_theory.json') else {'_intro': '', 'cards': []}
def th_card(c):
    op = f' <a class="open" href="prime2.html?{c["open"]}">open</a>' if c.get('open') else ''
    ch = f'<div class="chart">{open("theory/" + c["chart"]).read()}</div>' if c.get('chart') else ''
    return (f'<section class="card" id="th-{c["id"]}"><h2>{html.escape(c["name"])}{op}</h2><p class="lab">{c["label"]}</p>{ch}'
            f'<p><b>The idea.</b> {c["idea"]}</p><p class="num"><b>The test.</b> {c["test"]}</p><p class="ex"><b>Known?</b> {c["known"]}</p></section>')
THEORY = ('' if not TT['cards'] else f'''<section class="panel" id="theories"><h2 style="margin:0 0 8px;font-size:18px">Theories of our own ({len(TT['cards'])}, 16 November)</h2><p style="margin:0;color:#d6e2f2">{TT['_intro']}</p>
<div class="wrap"><table><thead><tr><th>Theory</th><th>The claim</th><th>The test</th></tr></thead><tbody>{''.join(f'<tr><td><a href="#th-{c["id"]}">{html.escape(c["name"])}</a><br><span class="lab">{c["label"]}</span></td><td>{c["short"]}</td><td>{c["verdict"]}</td></tr>' for c in TT['cards'])}</tbody></table></div></section>'''
        + ''.join(th_card(c) for c in TT['cards']))

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
.imgs.one {{ grid-template-columns:1fr; margin-top:10px; }}
.lab {{ color:#a5b4fc; font-size:13px; }}
.card p.lab {{ margin:0 0 6px; color:#a5b4fc; font-size:13px; }}
.chart {{ margin-top:12px; background:#0e1828; border:1px solid var(--edge); border-radius:8px; padding:8px; }}
@media (max-width:640px) {{ .imgs {{ grid-template-columns:1fr; }} }}
</style>
</head>
<body>
<header>
<h1>Patterns in the prime layouts</h1>
<p>Every layout of <a href="prime2.html">prime2.html</a>, photographed whole and zoomed in four times, with the exponent links off so that the pattern of the primes shows; on the layouts with straight lines the brightest line found is drawn in yellow by the page itself.</p>
<p>{T['_intro']['method']}</p>
<p><a href="#theories">Theories of our own ({len(TT['cards'])}, 16 November)</a> · <a href="#new">{len(NT['cards'])} number-theory strategies (16 November)</a> · <a href="#cos12big">the cosine rows zoomed out</a></p>
</header>
<main>
<section class="panel"><h2 style="margin:0 0 8px;font-size:18px">The answer in short</h2><p style="margin:0 0 10px;color:#d6e2f2">{T['_intro']['answer']}</p>
<div class="wrap"><table><thead><tr><th>Layout</th><th>The pattern</th><th>Verdict</th></tr></thead><tbody>{''.join(rows)}</tbody></table></div></section>
<section class="card" id="cos12"><h2>Deep dive: cosine rows, period 12 <a class="open" href="prime2.html?layout=cosrows&amp;period=12&amp;n=10000&amp;width=100&amp;links=none">open</a></h2>
<div class="imgs"><img src="prime2-patterns/cos12_w100_primes.jpg" alt="Cosine rows, period 12, width 100, primes" loading="lazy"><img src="prime2-patterns/cos12_w60_primes.jpg" alt="Cosine rows, period 12, width 60, primes" loading="lazy"><img src="prime2-patterns/cos12_twins_zoom.jpg" alt="Cosine rows, period 12, twin primes, zoomed" loading="lazy"><img src="prime2-patterns/cos12_w60_twins.jpg" alt="Cosine rows, period 12, width 60, twin primes" loading="lazy"></div>
<p><b>What you see.</b> Every prime from 5 on rides one of two lanes: up (+0.87) when p ≡ ±1 (mod 12), that is when 3 is a square mod p, and down (−0.87) when p ≡ ±5. With a width that is a multiple of 12 (60, lower left) every column keeps one height and the primes fill only the 16 columns coprime to 60; with width 100 (≡ 4 mod 12, top) each row shifts the wave by four places and the lanes braid. Twin primes (right) are green dashes at one height.</p>
<p><b>The lanes are equally full</b>, but not quite: up to 500,000 the down lane holds 20,806 primes and the up lane 20,730, because 1 is the only square mod 12 and the class 1 lags (Chebyshev’s bias, 73 behind the other three classes on average). The up lane first leads at p = 61,463 and leads at only 695 of 41,536 primes.</p>
<p><b>Consecutive primes zigzag.</b> The gap to the next prime decides the lane exactly: a gap ≡ 2, 10 or 12 (mod 12) keeps it, 4, 6 or 8 switches it. Gaps ≡ 6 are the commonest (26.8%), so consecutive primes switch lanes 54.6% of the time up to 500,000 (53.7% up to 10 million), not 50%. Deeper down, consecutive primes repeat their remainder mod 12 only 57% as often as independence would give (1,464 times for remainder 1 against 2,569), the bias found by Lemke Oliver and Soundararajan in 2016; it fades slowly (66% up to 10 million), as their theory says.</p>
<p><b>Twin primes share a lane, always</b>: a pair (p, p + 2) has p ≡ 11 (both up, like 11 and 13) or p ≡ 5 (both down, like 5 and 7). Up and down pairs are equally common, as Hardy–Littlewood predicts: 2,318 against 2,246 up to 500,000, and 29,497 against 29,482 up to 10 million.</p>
<p class="ex"><b>The exponents.</b> Of the 49 Mersenne exponents from 5 to 82,589,933 only 16 sit in the up lane, where the heuristic expects about 23 (P = 0.023 one-sided, 0.046 two-sided); both halves lean the same way (p = 0.06 below 100,000, 0.15 above), and no. 52, 136,279,841, is down too. This is one split of the classes mod 12 among the many tested since the first round (all classes mod 12 together: p = 0.095), so it is not significant, and no mechanism is known. A clean test lies ahead: if the lean were real, about two of every three future exponents would be ≡ ±5 (mod 12); the heuristic says about half.</p></section>
{THEORY}
{BIG}
{TEN}
{NEWS}
{''.join(cards)}
</main>
</body>
</html>
'''
assert '<!--' not in page and '/*' not in page
open('report/prime2-patterns.html', 'w', encoding='utf-8').write(page)
print('ok', len(page))
