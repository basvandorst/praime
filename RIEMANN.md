# Primes from the zeros of zeta, and zeros from the primes

> **In short:** Riemann's explicit formula reads Λ(n) from the zeros of ζ, but the plain step in ψ needs zeros up to about 2πn. A decoder that subtracts the echo of decided numbers, models every undecided number as a coin with the known prime density, fades out the top seventh of the zeros and searches the next 20 candidates exactly needs 0.075·πn on odd numbers and 0.045·πn on numbers coprime to 210. Below 10⁶ that is 3.6 and 2.6 zeros per prime instead of 74. Backwards, the prime powers up to 1.6·√(T/2π) fix every zero up to height T to 10⁻³; the 25 primes below 100 give all 10,142 zeros below 10,000, and those give back all 6,935 primes below 70,000.

**Page:** [prime.html](writeup/code/riemann/prime.html), an interactive page (transformations 125–140, section "Riemann's explicit formula")

**Code:** [code/riemann/](writeup/code/riemann/): the page itself ([prime.html](writeup/code/riemann/prime.html), the version of 11 November 2026) and the scripts that check every number below ([zeros.js](writeup/code/riemann/zeros.js), [zeros_check.py](writeup/code/riemann/zeros_check.py), [validate.js](writeup/code/riemann/validate.js), [primes_to_zeros.py](writeup/code/riemann/primes_to_zeros.py), [zeros_to_primes.js](writeup/code/riemann/zeros_to_primes.js)).

**Status:** computed. Nothing here is proved; every claim is a measurement with the test that produced it. The decoders assume nothing about RH; the zero finder and the backwards fit assume the zeros are simple and on the critical line, which holds at these heights.

## 1. The readout

The explicit formula ψ(x) = x − Σ_ρ x^ρ/ρ − ln 2π − ½ ln(1 − x⁻²) has a step of height Λ(n) at every integer n. Cut the sum over zeros at height T = Ωx and differentiate: what is left is the sum of all steps, each smeared by a sinc of width π/Ω,

```
r(x) = (π/Ω) (1 − (2/√x) Σ_{γ ≤ Ωx} w(γ/Ωx) cos(γ ln x))  ≈  Σ_m Λ(m) √(x/m) K(Ω x ln(m/x))
```

with K = sinc for a sharp cut. The decoder walks up from small n. At each candidate n everything below n is already decided, so its echo Σ_{m decided} Λ(m) √(x/m) K(…) is subtracted, and only the undecided numbers remain in r(x).

The weight w fades out the top seventh of the zeros with a cosine. Its kernel, used in the code for the echo, is

```js
// prime.html, decRT: kernel of the faded cut, z = Ω x ln(m/x), a = 7π
function decRT(x,m,om){const z=om*x*Math.log(m/x),a=7*Math.PI,az=Math.abs(z),r=Math.sqrt(x/m);if(az<1e-9)return r*13/14;if(Math.abs(az-a)<1e-6)return r/28;return r*(Math.sin(az)+Math.sin(6*az/7))*a*a/(2*az*(a*a-az*az));}

// prime.html, tapSum: the faded zero sum (the page uses block moments for speed, same value)
function tapSum(Z,lim,T,L,tap){
 let s=0;if(!tap){for(let i=0;i<lim&&Z[i]<=T;i++)s+=Math.cos(Z[i]*L);return s;}
 const T0=6*T/7,w=7*Math.PI/T;let i=0;for(;i<lim&&Z[i]<=T0;i++)s+=Math.cos(Z[i]*L);for(;i<lim&&Z[i]<=T;i++)s+=0.5*(1+Math.cos(w*Z[i]))*Math.cos(Z[i]*L);return s;
}

// prime.html, decSample: r(x) minus the echo of everything decided (undecided numbers at their mean)
function decSample(st,x){
 let v=st.D.get(x);if(v!==undefined)return v;
 const om=decOm(st,x),T=om*x,L=Math.log(x),G=st.G;let s=0;
 s=st.tap?zeroSumTaper(G,st.K,T,L):zeroSum(G,st.K,T,L);
 v=(Math.PI/om)*(1-2*s/Math.sqrt(x));
 for(let m=Math.max(2,x-st.P);m<=x+st.P;m++){const f=decFill(st,m);if(f!==0)v-=f*decK(st,x,m,om);}
 st.D.set(x,v);return v;
}
```

A sharp cut leaks the whole unknown future into every sample through the 1/d tails of the sinc. The fade removes that leak: with every prime known, the model error falls from 0.14–0.66 to an rms of 0.006.

## 2. The decoder

For each candidate n (odd numbers for 134; numbers coprime to 210 for 139; coprime to the first k primes for 140):

1. Take the next 20 undecided candidates U as unknowns, each Λ = 0 or ln m.
2. Use the samples r(x) for x from n − 96 to the end of the window.
3. Model every candidate beyond the window as a coin (the next 160 numbers for 134, 300 for 139). Its variance is p(1 − p) ln²m, with p the candidates' mean density (1/ρ)/ln m. Add a noise of 0.01, and whiten the samples with the Cholesky factor of this covariance.
4. QR the 20 echo columns, and find the 0/1 pattern with the least whitened residual by depth-first search with bounds (a sphere decoder, 5,000–20,000 nodes).
5. Keep only the answer for n, and slide one candidate on.

A prime power is not a separate unknown: its value is ln of its root, fixed only once the root is decided. The same solver runs the class formulas 135–137, where each residue class is decoded on its own channel.

```js
// prime.html, winSolve: exact window decoder (r: samples minus the echo, U: unknowns, far: modelled future)
function winSolve(xs,oms,r,U,far,mean,kern){
 const n=xs.length,k=U.length,C=[],col=new Float64Array(n);
 for(let i=0;i<n;i++){const row=new Float64Array(n);row[i]=CLS.sig;C.push(row);}
 for(const m of far){const L=Math.log(m),p=Math.min(0.95,mean/L),s2=p*(1-p)*L*L;for(let i=0;i<n;i++)col[i]=kern(xs[i],m,oms[i]);for(let i=0;i<n;i++){const ci=s2*col[i];if(ci!==0){const row=C[i];for(let q=0;q<=i;q++)row[q]+=ci*col[q];}}}
 for(let i=0;i<n;i++){const ri=C[i];for(let q=0;q<=i;q++){let acc=ri[q];const rq=C[q];for(let t=0;t<q;t++)acc-=ri[t]*rq[t];ri[q]=q===i?Math.sqrt(Math.max(acc,1e-300)):acc/rq[q];}}
 const white=v=>{const y=new Float64Array(n);for(let i=0;i<n;i++){let acc=v[i];const ri=C[i];for(let t=0;t<i;t++)acc-=ri[t]*y[t];y[i]=acc/ri[i];}return y;};
 const t0=Float64Array.from(r),A=[];for(let j=0;j<k;j++){const m=U[k-1-j],L=Math.log(m),v=new Float64Array(n);for(let i=0;i<n;i++){const q=kern(xs[i],m,oms[i]);v[i]=L*q;t0[i]+=mean*q;}A.push(v);}
 const t=white(t0),Q=[],R=[];for(let j=0;j<k;j++)R.push(new Float64Array(k));
 for(let j=0;j<k;j++){const v=white(A[j]);for(let i=0;i<j;i++){let d=0;const qi=Q[i];for(let z=0;z<n;z++)d+=qi[z]*v[z];R[i][j]=d;for(let z=0;z<n;z++)v[z]-=d*qi[z];}
   let nn=0;for(let z=0;z<n;z++)nn+=v[z]*v[z];nn=Math.sqrt(nn)||1e-150;R[j][j]=nn;for(let z=0;z<n;z++)v[z]/=nn;Q.push(v);}
 const y=new Float64Array(k);for(let j=0;j<k;j++){let d=0;const qj=Q[j];for(let z=0;z<n;z++)d+=qj[z]*t[z];y[j]=d;}
 const cur=new Uint8Array(k),best=new Uint8Array(k);let bestCost=Infinity,nodes=0;
 const rec=(lev,cost)=>{if(cost>=bestCost||nodes>300000)return;if(lev<0){bestCost=cost;best.set(cur);return;}nodes++;
   let s0=y[lev];const rl=R[lev];for(let j=lev+1;j<k;j++)if(cur[j])s0-=rl[j];
   const e0=s0*s0,e1=(s0-rl[lev])*(s0-rl[lev]);
   if(e0<=e1){cur[lev]=0;rec(lev-1,cost+e0);cur[lev]=1;rec(lev-1,cost+e1);}else{cur[lev]=1;rec(lev-1,cost+e1);cur[lev]=0;rec(lev-1,cost+e0);}cur[lev]=0;};
 rec(k-1,0);return best[k-1]===1;
}
```

## 3. Heights and zeros per prime

The zero height each method needs to read every prime up to n. Every row is 0 wrong on its validation set (section 9).

| Nr | Method | Zeros up to |
|---|---|---|
| 125 | step of ψ between n − ½ and n + ½ | about 2πn |
| 126 | Λ̂(n) = 2 − (4/√n) Σ cos(γ ln n) at odd n | 0.5·πn |
| 131 | 126 minus the echo of its own earlier outputs | 0.4·πn |
| 127 | 126 on 6k ± 1 | 0.333·πn |
| 133 | two samples below n: a(n) = (9π/8) r(n − 2) + (15π/8) r(n − 6) | 0.25·πn |
| 128 | all on/off patterns of the next 12 odd numbers | 0.2·πn |
| 130 | 128 on sieve survivors | 0.6·π·ρ·n |
| **134** | exact window decoder on odd numbers, no sieve | **0.075·πn**, 0.06·πn above 200,000 |
| **139** | 134 on numbers coprime to 210 | **0.045·πn** |
| **140** | 139 with the first k primes as divisibility facts | **0.23·π·ρ·n** (0.028 for k = 25) |
| 135 | 134 per class mod 4, ζ and L(s, χ₋₄) | 0.045·πn |
| 136 | per class mod 24, 8 real characters | 0.008·πn |
| 137 | per class mod 840, 192 characters | 0.0004·πn (0.00032 at 80,000) |
| 138 | picks 137, 136 or 139 | |

Zeros needed for every prime below 10⁶ (π(10⁶) = 78,498):

| Method | Zeros | Per prime | Fewer than plain ψ |
|---|---|---|---|
| plain explicit formula (zeros to about 3n) | 5,799,141 | 74 | 1× |
| 134, odd numbers, no sieve | 279,269 | 3.6 | 21× |
| 139, coprime to 210 | 202,979 | 2.6 | 29× |
| 140, k = 25 (primes to 97) | 118,225 | 1.5 | 49× |

140 leaves the easy part to a sieve; the zeros decide only what no sieve can see, whether a survivor is prime or the product of two large primes. Without a sieve the best is 134.

## 4. Laws: why it works and where it stops

- **Blur and candidate spacing.** Zeros up to T = Ωn blur the primes over π/Ω numbers. The decoder reads at a blur of about 5 candidate spacings, so the height needed is about 0.2·πn ÷ spacing. The measured spacings and heights: 134 (spacing 2) 0.075; 139 (4.4) 0.045; 135 (4) 0.045; 136 (24) 0.008; 137 (840) 0.0004.
- **Divisibility facts.** Leaving the multiples of the first k primes out lowers the height in proportion to the survivors' share ρ: 0.23·π·ρ·n, checked for k = 2 to 25.
- **Where it fails.** With the true past given, every wrong decision sits where the next candidate is 2 away and exactly one of the two is prime (twin positions). At spacing 4 or more there were 0 errors down to 0.10·πn. The height is set by the tightest cluster of candidates, not by the density of primes.
- **Where the information ends.** Oracle test, with the true decisions known to the right:
  - the true pattern beats every single flip down to 0.02·πn, and every twin swap with 80 numbers known to the right;
  - started from the truth, a second pass with windows of 16 stays right at 0.03 and 0.025 and drifts at 0.02.

  So the zeros carry the primes down to about 0.02–0.025·πn. The one-pass decoder stops near 0.035 because one error corrupts the echo for every later decision, and errors come in avalanches.
- **Counting.** At height Ωn the zeros give Ω/π independent real numbers per integer. At 0.045·πn that is 0.045 numbers per integer against 0.23 candidates per integer, so each number settles about 5 yes/no questions. That works only because the unknowns are 0 or ln m and the far future has a known spread.
- **Small n needs more:** 0.3·πn below 600 for 134 and 139. 136 needs 0.008·πn (not 0.0067) because of one spot near 4651, where primes are dense.
- **Precision.** Rounding the zeros to ε adds noise with rms (π/Ω)(2/√x)·ln x·(ε/√12)·√(N/2), where N is the number of zeros used (checked within 20 %). Against the model noise 0.1 that allows ε ≈ 10⁻³:
  - rounding to 10⁻⁴, 3·10⁻⁴ or 10⁻³ changes nothing at 0.045;
  - 3·10⁻³ needs 0.055;
  - 10⁻² fails at 0.055.

  So 3 decimals per zero are enough, and more decimals cannot lower the height.

## 5. The zeros

The page computes its own zeros: 2,000,000 of them, up to height 1,131,944, in about 32 s in the browser engine. The count equals N(T). The method:

- Riemann–Siegel Z(t), with each term advanced by rotation and resynchronised every 150 steps;
- a sign scan in steps of 0.06, with a close-pair check: under RH there is exactly one extremum between consecutive zeros, so three same-sign samples whose middle one is smallest in size prove a missed pair;
- each zero located by a Taylor step from the scan state, with Illinois false position as a fallback.

Until the last change in section 7 the remainder had only the C0 term, and zeros at heights 100–1500 were off by up to 4·10⁻⁴. With C1 and C2, the worst of the zeros checked against `mpmath.zetazero` (indices 31 to 2,000,000) is 6.5·10⁻⁷. The coefficients are the Taylor series in z = 2p − 1, fitted with mpmath to 10⁻¹⁹; the C1 coefficients agree with the published ones.

```js
// prime.html, rsRem: Riemann–Siegel remainder with C0, C1, C2 (a = sqrt(t/2π), m = floor(a))
function rsRem(a,m){const p=a-m,z=2*p-1,z2=z*z;let c1=0,c2=0;for(let i=13;i>=0;i--){c1=c1*z2+RSC1[i];c2=c2*z2+RSC2[i];}return (m%2?1:-1)*Math.pow(a,-0.5)*(Math.cos(2*Math.PI*(p*p-p-1/16))/Math.cos(2*Math.PI*p)+z*c1/a+c2/(a*a));}
```

The first 30 zeros are a table with 9 decimals. The sums over up to 10⁶ zeros at many nearby x use block moments: blocks of width √(0.3·top), with 8 Taylor moments around a reference point.

## 6. Backwards: zeros from the primes

Sum the prime powers up to X with the same faded weight (w = 1 on [0, 6/7], then a cosine fade to 0 at 1):

```
N_X(t) = θ(t)/π + 1 − (1/π) Σ_{n ≤ X} Λ(n) w(ln n / ln X) sin(t ln n) / (√n ln n)
```

By Weil's explicit formula this is C + Σ_γ Kc(ln X·(t − γ)), where Kc(s) = ½ + (1/π)∫₀¹ w(v) sin(sv)/v dv is a smoothed unit step. On top comes a pole term, about 1 + O(√X/t³). So N_X is a known sum of smoothed steps at the zeros, and a least-squares fit of the step positions (a deconvolution) finds them:

- the crossings of k − ½ alone give the zeros to only about 10⁻²;
- the fit gives 10⁻⁵ to 5·10⁻⁴.

```python
# primes_to_zeros.py: the fit of the step positions (Kc and its derivative k1 tabulated)
def fit(ts, Ns, p0, U):
    x0 = np.append(p0, np.median(Ns - np.sum(kc(U * (ts[:, None] - p0[None, :])), axis=1)))
    def f(p):
        return p[-1] + np.sum(kc(U * (ts[:, None] - p[None, :-1])), axis=1) - Ns
    def J(p):
        return np.hstack([-U * k1(U * (ts[:, None] - p[None, :-1])), np.ones((len(ts), 1))])
    r = least_squares(f, x0, jac=J, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=300)
    return np.sort(r.x[:-1])
```

- **Threshold.** Every zero comes out to under 10⁻³ when X ≥ 1.6·√(T/2π). This holds for X = 30, 100 and 300, at heights up to 565,000. At 1.3 it partly works, at 1.0 it fails. This is a count: N_X has ln X/π independent values per unit of t, against ln(T/2π)/2π zeros per unit of t. So X must be about √(T/2π), the same length as in the Riemann–Siegel formula.
- **Chain.** The 25 primes below 100 (35 prime powers) give all 10,142 zeros below height 10,000: rms error 6.7·10⁻⁵, max 1.7·10⁻³. The four zeros off by more than 10⁻³ are two close pairs, #4087/4088 and #6709/6710. Fed into 139, these zeros give all 6,935 primes below 70,000 with 0 wrong, the same as the true zeros.
- **Mirror.** The weak spots mirror each other: close pairs of zeros in this direction, twin positions in the other.
- **Weights.** The Fejér or mollifier weight 1 − ln n/ln X, which Selberg and Zhang use in proofs, gives only 10⁻² here. A positive kernel helps proofs, not resolution.

## 7. Change log: every change, what works, how it was tested

All changes are to [prime.html](writeup/code/riemann/prime.html). "Circles" are stretches of 100–1000 numbers starting at a random S, with the true primes below S known, decoded with just enough zeros for the stated height. All checks are against trial division.

| Date | Change | Works | How it was tested |
|---|---|---|---|
| 6 Nov | Explicit-formula page: ψ rebuilt from zeros, strip of real against predicted primes, method selector (zeros, bootstrap sieve, Fermat, Miller–Rabin, Wilson, sieve), ranking formula, recursive bootstrap sieve | yes (display) | visual; no error counts claimed |
| 10 Nov | Prime Connection Laboratory; ψ ring with zeros slider; 125 step of ψ | yes, about 2πn | first misread integer at n ≈ 0.19–0.22·T for 100–2000 zeros |
| 10 Nov | 126 pointwise at odd n; 127 on 6k ± 1 | yes, 0.5 and 0.333·πn | runs over all odd n / 6k ± 1 in range |
| 10 Nov | 128 joint on/off fit of the next 12 odd numbers | yes, 0.2·πn | 0 errors on all odd n ≤ 15,600 from 10,000 zeros; first error near n = 2γ_K |
| 10 Nov | 129 least witness of 1 + 8k·p! | yes, but it is Wilson's theorem in disguise | theorem; no zeros |
| 10 Nov | decoder against 24 decided positions; 130 on sieve survivors | yes, 0.6·π·ρ·n | runs on survivors of the first k primes |
| 10 Nov | 131 recursion (126 minus its own echo) | yes, 0.4·πn | runs over odd n |
| 10 Nov | 132 Golomb's ζ limit made exact with s = p + 16 | yes, but no zeros (needs the primes so far) | exact big-integer arithmetic |
| 11 Nov | 133 two-point formula; 134 linear left-sample fit | yes, 0.25 and 0.2·πn | 133: 0 errors on odd n 500–12,400; 134: 0 on all odd n ≤ 15,600 from 10,000 zeros |
| 11 Nov | 135 with L(s, χ₋₄); 136 with the 8 real characters mod 24; packed tables to height 2198 | yes, 0.1 and 0.02·πn | 136: every prime to 34,800 in one chain |
| 11 Nov | mod-4 prime race ring; Liouville λ(n) ring | yes, λ at 0.25·πn | 0 wrong on 675 odd n (300–1050, 1500–2100); fails at 0.2 |
| 11 Nov | 134 self-correction at close calls | yes, 0.167·πn | every prime to 18,700 from 10,000 zeros, to 49,000 from 30,000 |
| 11 Nov | 137 with the 192 characters mod 840; 138 automatic; self-correction for 136 and 137 | yes | 137 every prime 4000–31,700 from zeros below height 80; 136 0.0167·πn to 41,500 |
| 11 Nov | zero finder 5× faster; block-moment sums (about 80× faster); 1,000,000 zeros with the close-pair check | yes | count equals N(T); a plain scan at 0.02 missed 5 pairs in the first 500,000 |
| 11 Nov | 134 at 0.15·πn above 50,000; block sums sized to the height in use | yes | 58 stretches from 50,000 to 1,269,000, about 4,400 primes, 0 wrong |
| 11 Nov | 139 on numbers coprime to 210 | yes, 0.13·πn | 150 random stretches, 7,141 primes, 0 wrong |
| 11 Nov | 2,000,000 zeros, each by a Taylor step | yes | agrees with false position within 2.6·10⁻⁷; count equals N(T); 32 s |
| 11 Nov | 139 with the top tenth of the zeros faded, 8 unknowns, two-level look-ahead | yes, 0.1·πn | 150 stretches clean at 0.095 |
| 11 Nov | full model of the unknown future (every candidate a coin), top seventh faded | yes: 139 0.065, 134 0.1 (0.3 below 600); above 200,000 0.055 and 0.085 | 400 random circles to 232,000 (17,725 primes) and 40 to 5.5 million, 0 wrong |
| 11 Nov | class formulas 135–137 with the full model | yes: 0.06, 0.01; 137 to 80,000 | 137: 200 circles, 0 wrong |
| 11 Nov | exact window decoder (sphere decoder over 20 candidates) for 134 and 139 | yes: 0.075 / 0.06 and 0.045 | 100 random circles to 300,000 clean at 0.04; page uses 0.045 for margin |
| 11 Nov | 140 divisibility facts | yes: 0.23·π·ρ·n | k = 2 … 25: 0 wrong from 1 to 3000 and on 30 circles each to 300,000 |
| 11 Nov | Riemann–Siegel C1 and C2; exact window decoder in 135–137; note on primes → zeros | yes: zeros to 6.5·10⁻⁷; 135 0.045, 136 0.008, 137 0.0004 | see section 9: validate.js, zeros_check.py; 135 1–8,100 and 136 4,000–85,000 in one run each |

## 8. What was tried and did not work

| Idea | Result | Test |
|---|---|---|
| more decimals per zero | no gain: 10⁻³ is enough | rounding sweep 10⁻⁴ … 10⁻², 24 circles; noise formula checked |
| coarse zeros (fewer bits) | costs height: 0.01 → 0.08·πn, 0.03 → > 0.11 | same bench |
| Hann, Welch, Fejér, wider fades | cost more resolution than they save | 36 fixed stretches |
| 40–60 candidate windows, K-best search | no gain over 20 searched exactly | same |
| density prior 1/ln n in the search | no change | 60 stretches at 0.035 and 0.03 |
| second pass with all neighbours fixed | no gain (0.035: 19 → 23 wrong; 0.03: 228 → 225) | 60 stretches |
| samples right of the window | no gain | 60 stretches |
| re-deciding the last 4 or 8 decisions | worse | 60 stretches |
| backtracking where the residual energy rises | worse (19 failing stretches instead of 1 at 0.035) | 60 stretches; raw energy also prefers wrong patterns in correct stretches |
| synthetic error-free data | the same errors as real zeros | at every band: the limit is the decoder's information, not the zeros' precision |
| exact 0/1 program over 30 unknowns (milp) | never error-free below the linear formula | scipy milp |
| prime-rich polynomials, Kronecker witness, Golomb | all sieves in disguise | survey of 132 maps against a deep sieve baseline |
| Fejér weight in the backwards fit | 10⁻² instead of 10⁻⁴ | windows at X = 30, 100, 300 |

**Works but not on the page yet:** the M-algorithm keeps the 4 best decision histories. It compares them by the summed whitened window cost, and branches when the other answer costs less than 50 more. At 0.035·πn it takes 19 wrong to 0; at 0.03 it takes 228 to 61 (60 stretches); 8 histories are no better. In 9 of 10 failing stretches the truth fits the data 4–10 times better than the decoder's answer. So the gap between 0.045 and about 0.025 is a search problem, not missing information. On the page it would need decisions delayed until the histories agree, and 3–4 times the time. Expected gain: 139 from 0.045 to about 0.04.

## 9. Validation

Run from the repository root (node 22, Python with numpy, scipy and mpmath):

```
node writeup/code/riemann/zeros.js 1000000 /tmp/z1m.f64
node writeup/code/riemann/zeros.js 10142 /tmp/z10k.f64
python3 writeup/code/riemann/zeros_check.py /tmp/z1m.f64 20
node writeup/code/riemann/validate.js
python3 writeup/code/riemann/primes_to_zeros.py threshold /tmp/z1m.f64
python3 writeup/code/riemann/primes_to_zeros.py all 100 10150 10000 /tmp/z10k.f64 /tmp/z_from100.f64
node writeup/code/riemann/zeros_to_primes.js /tmp/z_from100.f64 70000
```

Output on 11 November 2026:

```
$ node writeup/code/riemann/zeros.js 1000000 /tmp/z1m.f64
1000000 zeros, the last at 600269.677012 [13.1 s]
$ node writeup/code/riemann/zeros.js 10142 /tmp/z10k.f64
10142 zeros, the last at 9998.850397 [0.1 s]
$ python3 writeup/code/riemann/zeros_check.py /tmp/z1m.f64 20
1000000 zeros, 43 checked against mpmath.zetazero: max error 6.51e-07 at zero #31 [19 s]
$ node writeup/code/riemann/validate.js
134: 1-4000: 0 wrong of 551 primes ; 100 random circles 1-199000: 0 wrong of 3404  [148 s]
135: 1-4000: 0 wrong of 551 primes ; 100-8100 in one run: 0 wrong of 993 ; 60 random circles 100-23000: 0 wrong of 2605  [252 s]
136: 1-4000: 0 wrong of 551 primes ; 4000-85000 in one run: 0 wrong of 7727 ; 60 random circles 100-85000: 0 wrong of 2261  [92 s]
137: 4000-20000 in one run: 0 wrong of 1712 ; 60 random circles 4000-79500: 0 wrong of 2321  [21 s]
138: 1-4000: 0 wrong of 551 primes ; 60 random circles 1-199000: 0 wrong of 1979  [25 s]
139: 1-4000: 0 wrong of 551 primes ; 100 random circles 1-199000: 0 wrong of 3548  [53 s]
140k2: 1-3000: 0 wrong of 431 primes ; 30 random circles 600-196442: 0 wrong of 1059  [25 s]
140k4: 1-3000: 0 wrong of 431 primes ; 30 random circles 600-286890: 0 wrong of 1038  [18 s]
140k11: 1-3000: 0 wrong of 431 primes ; 30 random circles 600-300000: 0 wrong of 1031  [22 s]
140k25: 1-3000: 0 wrong of 431 primes ; 30 random circles 600-300000: 0 wrong of 1001  [28 s]
$ python3 writeup/code/riemann/primes_to_zeros.py threshold /tmp/z1m.f64
X = 30: X = 2.0 sqrt(T/2pi) (T = 1414): max error 4.5e-04 (crossings 2.5e-01); X = 1.6 sqrt(T/2pi) (T = 2209): max error 3.5e-04 (crossings 2.2e-01); X = 1.3 sqrt(T/2pi) (T = 3346): max error 2.9e-04 (crossings 1.5e-01); X = 1.0 sqrt(T/2pi) (T = 5655): max error 1.7e-01 (crossings 2.2e-01)
X = 100: X = 2.0 sqrt(T/2pi) (T = 15708): max error 3.4e-04 (crossings 1.7e-01); X = 1.6 sqrt(T/2pi) (T = 24544): max error 5.5e-04 (crossings 1.3e-01); X = 1.3 sqrt(T/2pi) (T = 37179): max error 1.1e-02 (crossings 1.8e-01); X = 1.0 sqrt(T/2pi) (T = 62832): max error 1.1e-01 (crossings 2.0e-01)
X = 300: X = 2.0 sqrt(T/2pi) (T = 141372): max error 2.8e-04 (crossings 1.2e-01); X = 1.6 sqrt(T/2pi) (T = 220893): max error 3.7e-04 (crossings 9.7e-02); X = 1.3 sqrt(T/2pi) (T = 334608): max error 2.7e-03 (crossings 1.1e-01); X = 1.0 sqrt(T/2pi) (T = 565487): max error 1.3e-01 (crossings 1.1e-01)
$ python3 writeup/code/riemann/primes_to_zeros.py all 100 10150 10000 /tmp/z10k.f64 /tmp/z_from100.f64
prime powers up to 100: 35; zeros recovered up to 10000: 10142 (true count 10142) [425 s]
max error 1.74e-03, rms 6.68e-05; above 1e-3: [4087, 4088, 6709, 6710]
$ node writeup/code/riemann/zeros_to_primes.js /tmp/z_from100.f64 70000
10142 zeros (the last at 9998.85): read 6935 primes up to 69997, 0 wrong; primes up to 70000: 6935 [80 s]
```

## 10. Hypotheses and what would settle them

| Hypothesis | Evidence so far | Confirms it | Refutes it | What is needed |
|---|---|---|---|---|
| **H1** Blur law: the height needed is about 0.2·πn ÷ candidate spacing (a blur of about 5 spacings) | spacings 2, 4, 4.4, 24 and 840 (section 4) | a new spacing on the predicted line, e.g. the 32 characters mod 120 (spacing 120 → about 0.0017·πn) | a spacing that needs clearly more or less | zeros of the 32 L-functions mod 120 (built off the page before), and the class decoder of 136 with 32 sets |
| **H2** Divisibility law 0.23·π·ρ·n | k = 2 … 25 up to 300,000 | the same factor for k > 25 and for n of 10⁶–10⁷ | a drift of the factor with k or n | more zeros: n = 10⁷ at k = 25 needs height about 8.7·10⁵ (the page has 1.13·10⁶) |
| **H3** The zeros carry the primes down to about 0.02–0.025·πn; the gap from 0.045 is a search problem | oracle tests; truth fits better in 9 of 10 failures; M-algorithm 19 → 0 errors at 0.035 | a decoder that is clean at 0.035 or lower on 100 random circles to 300,000 | many failures where the wrong pattern has the lower whitened cost (then the information is missing), or any decoder clean below 0.02 (then the limit is lower) | the M-algorithm on the page decoder with delayed decisions; the oracle test on more stretches |
| **H4** Backwards: the prime powers up to c·√(T/2π), c ≈ 1.6, fix every zero up to height T; the threshold is a count | X = 30, 100, 300, heights to 5.7·10⁵, one window of 60 zeros per point | many windows per (X, c), larger X, much larger T; the threshold moves with the fade as the count predicts (fading from 0.9 instead of 6/7 should lower c) | success well below c = 1 (better than the count) or failure above 1.6 at large X | true zeros at large heights, e.g. A. Odlyzko's published tables near the 10¹²-th zero (height about 2.7·10¹¹): the count predicts that the primes up to about 3.3·10⁵ fix those zeros |
| **H5** Close pairs of zeros and twin positions are the same obstacle seen from both sides | the only zeros above 10⁻³ in the chain are two close pairs | the backwards error grows as the gap shrinks, in the same way as the forward error at spacing 2 | close pairs recovered as well as the rest | the gap of every zero in the chain against its error (data already there) |

**The Riemann Hypothesis itself** is not tested by any of this. The zero finder assumes the zeros are on the line, and checks it in the usual way: it counts sign changes of Z(t) and compares with N(T). Both agree up to height 1.13·10⁶; RH is known by others to hold much higher. Computation can only refute RH, by a height where the sign changes fall short of N(T), and that would have to be far above anything a page can reach. A confirmation needs a proof, not a computation.

## 11. Open

- The M-algorithm on the page (section 8).
- A prior for clusters from the Hardy–Littlewood singular series. The admissible tuples of the bounded-gaps work (Zhang, Annals 179, 2014) are exactly the clusters where the decoders fail; this is untested.
- A second round of the chain: the primes below 70,000 would fix the zeros up to about 2π(70000/1.6)² ≈ 1.2·10¹⁰, beyond what can be computed here.

## References

- B. Riemann, Über die Anzahl der Primzahlen unter einer gegebenen Grösse (1859); H. von Mangoldt (1895) for ψ(x).
- A. Weil's explicit formula, used here for the kernel of N_X.
- The Riemann–Siegel formula with the correction terms C0, C1, C2; the coefficients here were fitted with mpmath and agree with the published ones for C1.
- Y. Zhang, Discrete mean estimates and the Landau–Siegel zero, arXiv:2211.02515 (v1, 4 Nov 2022). Prop. 2.2: under a Siegel zero the zeros in his box have gaps π/log P, the resolution scale π/ln X of section 6. Lemma 3.3: the large sieve over characters, the reason more characters act as near-independent channels (135–137). Eq. (2.23)–(2.25): the mollifier weights tested in section 6.
- Y. Zhang, Bounded gaps between primes, Ann. of Math. 179 (2014), 1121–1174; only the abstract was read.
- Details and dates of each step are in [IDEAS.txt](IDEAS.txt) and [FINDINGS.txt](FINDINGS.txt) (entries of 10–11 November 2026).
