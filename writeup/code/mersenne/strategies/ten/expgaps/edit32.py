"""Two layouts for the gaps between Mersenne exponents (19 November 2026): expgaprows and expgapmod.
The fixed numbers in the notes come from ~/scratch/expgaps/mexp_gaps.json."""
import json, math, os
J = json.load(open(os.path.expanduser('~/scratch/expgaps/mexp_gaps.json')))
lx = open('lx.js', encoding='utf-8').read()


def rep(a, b, c=1):
    global lx
    assert lx.count(a) == c, (a[:70], lx.count(a))
    lx = lx.replace(a, b)


s5 = J['sizes_from_5']
mods = {r['m']: r for r in J['mods']['rows']}
pairs = J['mods']['pairs']
assert pairs == 42 and abs(s5['mean'] - 0.4998) < 0.01
rep('LXM=["nocarry",', 'LXM=["expgaprows","expgapmod","nocarry",')
rep('const LXN=["nocarry",', 'const LXN=["expgaprows","expgapmod","nocarry",')
rep('LXP=["nocarry",', 'LXP=["expgapmod","nocarry",')
rep('LXSC={"nocarry":0,', 'LXSC={"expgaprows":0,"expgapmod":0,"nocarry":0,')
code = r''' else if(name==="expgaprows"){const E=MEXP.filter(k=>k<=N),C=48,RH=4,W=Math.ceil(2.2*C),row=new Int32Array(N+1).fill(-1),xs=new Float64Array(N+1);let j=0;
  for(let n=2;n<=N;n++){while(j+1<E.length&&E[j+1]<=n)j++;if(!isP(n))continue;row[n]=j;xs[n]=Math.min(W-1,Math.log2(n/E[j])*C);}
  hid(W,E.length*RH,n=>row[n]<0?null:[xs[n],row[n]*RH+1]);out.parity=null;out.noexp=false;out.dw=2;
  const r=[];for(let i=0;i+1<E.length;i++)r.push([Math.log2(E[i+1]/E[i]),E[i],E[i+1]]);const mean=r.reduce((a,x)=>a+x[0],0)/Math.max(1,r.length);let lo=r[0],hi=r[0];for(const x of r){if(x[0]<lo[0])lo=x;if(x[0]>hi[0])hi=x;}
  out.tipFull=n=>row[n]<0?"":`<b>${f(n)}</b> (prime) lies in row ${row[n]+1}, between the exponents ${f(E[row[n]])} and ${row[n]+1<E.length?f(E[row[n]+1]):"the next one, above N"}: ${(xs[n]/C).toFixed(2)} doublings from the start of the row`;
  Object.assign(sh,{ne:E.length,nr:r.length,mean,lo,hi,C});}
 else if(name==="expgapmod"){const M=Math.max(2,Math.min(1000,Math.round(P))),E=MEXP.filter(k=>k<=N),idx=new Map();E.forEach((e,i)=>idx.set(e,i));
  hid(M*3,E.length*3,n=>idx.has(n)?[(n%M)*3+1,idx.get(n)*3+1]:null);out.parity=null;out.noexp=false;out.dw=3;
  let phi=M,t=M;for(let q=2;q*q<=t;q++)if(t%q===0){while(t%q===0)t/=q;phi-=phi/q;}if(t>1)phi-=phi/t;
  const A=MEXP.slice(0,51).filter(p=>p>=61),n=A.length-1;let hit=0;for(let i=0;i<n;i++)if((A[i+1]-A[i])%M===0)hit++;
  const q=1/phi,lg=[0];for(let i=1;i<=n;i++)lg.push(lg[i-1]+Math.log(i));const pm=k=>Math.exp(lg[n]-lg[k]-lg[n-k]+k*Math.log(q)+(n-k)*Math.log(1-q));let up=0,dn=0;for(let k=0;k<=n;k++){if(k>=hit)up+=pm(k);if(k<=hit)dn+=pm(k);}
  let inN=0;for(let i=0;i+1<E.length;i++)if((E[i+1]-E[i])%M===0)inN++;
  out.tipFull=n=>{if(!idx.has(n))return "";const i=idx.get(n),g=i>0?n-E[i-1]:0;return `<b>${f(n)}</b>, exponent number ${i+1}: ${f(n)} ≡ ${n%M} (mod ${M})`+(i>0?`; gap from ${f(E[i-1])}: ${f(g)}, ${g%M===0?"a multiple of "+M:"remainder "+(g%M)}`:"");};
  Object.assign(sh,{M,phi,n,hit,ex:n/phi,pv:Math.min(1,2*Math.min(up,dn)),inN,ne:E.length});}
 else if(name==="bases"){'''
rep(' else if(name==="bases"){', code)
pv = {}
for m in (4, 5):
    from math import comb
    from sympy import totient
    q = 1 / int(totient(m))
    k = mods[m]['real']
    up = sum(comb(pairs, i) * q ** i * (1 - q) ** (pairs - i) for i in range(k, pairs + 1))
    dn = sum(comb(pairs, i) * q ** i * (1 - q) ** (pairs - i) for i in range(0, k + 1))
    pv[m] = min(1, 2 * min(up, dn))
minB = min(r['pB'] for r in J['mods']['rows'])
mA = [m for m in mods]
notes = r'''
  expgaprows:()=>`Only the primes are drawn. Each row runs from one Mersenne exponent to the next (row 1 from 2 to 3, row 2 from 3 to 5, …), on a doubling scale with ${h.C} cells per doubling, so the length of a row is log₂ of the ratio of two consecutive exponents. Choose Links: In order to join the exponents. Up to ${NN} there are ${h.ne} exponents and ${h.nr} whole rows: mean length ${h.mean.toFixed(2)} doublings, the longest ${h.hi[0].toFixed(2)} (from ${f(h.hi[1])} to ${f(h.hi[2])}), the shortest ${h.lo[0].toFixed(3)} (from ${f(h.lo[1])} to ${f(h.lo[2])}). The model of Lenstra, Pomerance and Wagstaff (a conjecture: log₂ of the exponents behaves like a random scatter with e^γ ≈ 1.78 points per doubling) gives rows of mean length __MODEL__ doublings, spread like an exponential law. For all 51 exponents that are certainly consecutive (GIMPS has checked everything up to 82,589,933 twice), from 5 on: mean __MEAN5__ (p = __PMEAN__), the spread fits the law (Kolmogorov–Smirnov p = __PKS__), and a long row is not followed by a short one more often than chance (correlation __RHO__, p = __PRHO__). The counts per doubling are, if anything, a little more even than random scatter (variance/mean __DISP__). The rows look like random lengths, and that is what the test finds.`,
  expgapmod:()=>`Only the Mersenne exponents up to N are drawn, one row each in order (2 at the top), in the column of their remainder mod P (now ${h.M}). Choose Links: In order: a vertical link means that the gap between two consecutive exponents is a multiple of ${h.M}. Up to ${NN}: ${h.inN} of ${Math.max(0,h.ne-1)} links are vertical. Over all 51 exponents that are certainly consecutive, from 61 on (${h.n} gaps), ${h.hit} gaps are multiples of ${h.M}; random primes of the same sizes would give ${h.n}/φ(${h.M}) = ${h.ex.toFixed(1)} (two-sided p = ${h.pv.toFixed(2)}). Try P = 4: __R4__ where __E4__ are expected (p = __P4__); the exponents are more often ≡ 1 (mod 4) than random primes (__C4__ of __N4__ from 61 on), a known tilt, because a factor 2p + 1 of 2^p − 1 is only possible for p ≡ 3 (mod 4). Try P = 5: __R5__ where __E5__ are expected (p = __P5__). Eleven values were tested in advance (P = 3, 4, 5, 7, 8, 9, 11, 12, 13, 16, 24), also against random primes that survive trial division of 2^p − 1 by 2kp + 1 for k ≤ 100: the smallest p was __MINB__. With eleven tries that is what chance gives. The gaps between the exponents show no remainder rule.`,
'''
repl = {'__MODEL__': '%.2f' % s5['model'], '__MEAN5__': '%.2f' % s5['mean'], '__PMEAN__': '%.2f' % s5['pMean'], '__PKS__': '%.2f' % s5['pD'],
        '__RHO__': ('%.2f' % s5['rho']).replace('-', '−'), '__PRHO__': '%.2f' % s5['pRho'], '__DISP__': '%.2f' % (J['doublings']['var'] / J['doublings']['mean']),
        '__R4__': str(mods[4]['real']), '__E4__': '%.1f' % (pairs / 2), '__P4__': '%.2f' % pv[4], '__C4__': str(J['class4']['real']), '__N4__': str(J['class4']['n']),
        '__R5__': str(mods[5]['real']), '__E5__': '%.1f' % (pairs / 4), '__P5__': '%.2f' % pv[5], '__MINB__': '%.2f' % minB}
for k, v in repl.items():
    assert notes.count(k) == 1, k
    notes = notes.replace(k, v)
rep('\n  depthrows:()=>`', '\n' + notes.strip('\n') + '\n  depthrows:()=>`')
open('lx.js', 'w', encoding='utf-8').write(lx)
o = open('opts.txt', encoding='utf-8').read().rstrip('\n')
assert 'expgaprows|' not in o
o += '\nexpgaprows|Gaps between Mersenne exponents: one row per gap, on a doubling scale\nexpgapmod|Gaps between Mersenne exponents: remainders mod P (a vertical link = a multiple of P)\n'
open('opts.txt', 'w', encoding='utf-8').write(o)
print('ok', repl)
