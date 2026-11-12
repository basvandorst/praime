# Primes from the zeros of zeta, and zeros from the primes

> **In short:** Riemann's explicit formula reads Λ(n) from the zeros of ζ, but the plain step in ψ needs zeros up to about 2πn. An exact window decoder (echo of decided numbers subtracted, undecided numbers as coins, top seventh of the zeros faded) needs 0.075·πn on odd numbers and 0.045·πn coprime to 210: 3.6 and 2.6 zeros per prime below 10⁶ instead of 74. Any band of zeros works too; its width and the pattern of its phases count, not its height: the 30,000 zeros Odlyzko published near the 10¹²-, 10²¹- and 10²²-th zero give every prime below 30,831 (64,853 with the wheel mod 210). Backwards, the prime powers up to 1.6–2·√(T/2π) fix the zeros near height T to 10⁻³, tested up to T = 2.7·10¹¹; the 25 primes below 100 give all zeros below 10,000, and those all primes below 70,000. The primes confirm Odlyzko's zeros to 10⁻⁸ and our own to 3·10⁻⁸.

**Page:** [prime.html](writeup/code/riemann/prime.html), an interactive page (transformations 125–140, section "Riemann's explicit formula")

**Code:** [code/riemann/](writeup/code/riemann/): the page itself ([prime.html](writeup/code/riemann/prime.html), the version of 11 November 2026) and the scripts that check every number below ([zeros.js](writeup/code/riemann/zeros.js), [zeros_check.py](writeup/code/riemann/zeros_check.py), [validate.js](writeup/code/riemann/validate.js), [primes_to_zeros.py](writeup/code/riemann/primes_to_zeros.py), [zeros_to_primes.js](writeup/code/riemann/zeros_to_primes.js); for sections 12 and 13 [band_decoder.py](writeup/code/riemann/band_decoder.py), [band_check.py](writeup/code/riemann/band_check.py), [zeros_vs_primes.py](writeup/code/riemann/zeros_vs_primes.py), [zeros_high.py](writeup/code/riemann/zeros_high.py), [band_phase.py](writeup/code/riemann/band_phase.py), [band_validate.sh](writeup/code/riemann/band_validate.sh), [band_validate_m.sh](writeup/code/riemann/band_validate_m.sh)). Odlyzko's tables of zeros near the 10¹²-, 10²¹- and 10²²-th zero are in [data/riemann/](writeup/data/riemann/).

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

- **Threshold.** Every zero comes out to under 10⁻³ when X ≥ 1.6·√(T/2π). This holds for X = 30, 100 and 300, at heights up to 565,000. At 1.3 it partly works, at 1.0 it fails. Near the 10¹²-th zero (height 2.68·10¹¹, X up to 412,787) the threshold is c ≈ 2 (section 13). This is a count: N_X has ln X/π independent values per unit of t, against ln(T/2π)/2π zeros per unit of t. So X must be about √(T/2π), the same length as in the Riemann–Siegel formula.
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
| band decoder: exact window of 16 candidates instead of 12 | worse: 6,031 → 5,483 | band near the 10²²-th zero |
| band decoder: samples at a third of the envelope scale with the noise per sample unchanged | worse: 6,043 → 5,483; with the noise times the step no loss | same band, scratch version of the decoder |
| band decoder: far coins counted at the end of the run | the reach moved with the end (lowpass 1,342: 8,267 at N = 16,000, 8,699 at 24,000); now 12 local envelope scales | lowpass of our own zeros |
| band decoder: 16 histories instead of 8, commit delay 48 instead of 24 | no gain (7,615), slightly worse (7,557) | band near the 10²²-th zero, scratch version of the decoder |

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

For sections 12 and 13 (12 November 2026; about 45 and 40 minutes on this machine, side by side). The decoders print the centre of each band, the first wrong decisions and the reach; band_check.py then checks every list against sympy.isprime:

```
sh writeup/code/riemann/band_validate.sh
sh writeup/code/riemann/band_validate_m.sh
```

Output of band_validate.sh:

```
2000000 zeros, the last at 1131944.471825 [50.0 s]
== bands of Odlyzko's zeros, odd candidates
1 band(s) at centre 2.67653e+11, total width 2567, 10000 zeros, wheel (2,): 1415 primes decided; first wrong [11791, 11793, 11799]; reach 11791 = 4.59 x width [33 s]
1 band(s) at centre 1.44177e+20, total width 1409, 10000 zeros, wheel (2,): 777 primes decided; first wrong [5901, 5903, 5905]; reach 5901 = 4.19 x width [21 s]
1 band(s) at centre 1.37092e+21, total width 1342, 10000 zeros, wheel (2,): 788 primes decided; first wrong [6031, 6033, 6035]; reach 6031 = 4.50 x width [21 s]
3 band(s) at centre 2.67653e+11, 1.44177e+20, 1.37092e+21, total width 5318, 30000 zeros, wheel (2,): 2557 primes decided; first wrong [22875, 22879, 22881]; reach 22875 = 4.30 x width [363 s]
== wheels
1 band(s) at centre 1.37092e+21, total width 1342, 10000 zeros, wheel (2, 3): 1232 primes decided; first wrong [10013, 10015, 10021]; reach 10013 = 7.46 x width [31 s]
1 band(s) at centre 1.37092e+21, total width 1342, 10000 zeros, wheel (2, 3, 5, 7): 1444 primes decided; first wrong [12107, 12109, 12113]; reach 12107 = 9.02 x width [23 s]
3 band(s) at centre 2.67653e+11, 1.44177e+20, 1.37092e+21, total width 5318, 30000 zeros, wheel (2, 3, 5, 7): 5204 primes decided; first wrong [50783, 50791, 50831]; reach 50783 = 9.55 x width [448 s]
== same width, low and high (own zeros)
1 band(s) at centre 0, total width 1341, 933 zeros, wheel (2,): 1084 primes decided; first wrong [8699, 8701, 8703]; reach 8699 = 6.49 x width [15 s]
1 band(s) at centre 500671, total width 1341, 2410 zeros, wheel (2,): 798 primes decided; first wrong [6129, 6131, 6133]; reach 6129 = 4.57 x width [17 s]
1 band(s) at centre 1.00067e+06, total width 1341, 2558 zeros, wheel (2,): 848 primes decided; first wrong [6569, 6571, 6573]; reach 6569 = 4.90 x width [15 s]
1 band(s) at centre 0, total width 2999, 2469 zeros, wheel (2,): 2246 primes decided; first wrong [19853, 19855, 19857]; reach 19853 = 6.62 x width [34 s]
1 band(s) at centre 501500, total width 2999, 5389 zeros, wheel (2,): 1700 primes decided; first wrong [14531, 14533, 14537]; reach 14531 = 4.85 x width [67 s]
1 band(s) at centre 1.0015e+06, total width 2999, 5719 zeros, wheel (2,): 1774 primes decided; first wrong [15193, 15195, 15197]; reach 15193 = 5.07 x width [63 s]
1 band(s) at centre 0, total width 9999, 10142 zeros, wheel (2,): 5724 primes decided; first wrong [56467, 56469, 56471]; reach 56467 = 5.65 x width [191 s]
1 band(s) at centre 505000, total width 9999, 17975 zeros, wheel (2,): 4218 primes decided; first wrong [40153, 40155, 40159]; reach 40153 = 4.02 x width [250 s]
1 band(s) at centre 1.005e+06, total width 9999, 19070 zeros, wheel (2,): 3985 primes decided; first wrong [37615, 37617, 37619]; reach 37615 = 3.76 x width [184 s]
2 band(s) at centre 0, 1.37092e+21, total width 2683, 10933 zeros, wheel (2,): 1746 primes decided; first wrong [14921, 14923, 14929]; reach 14921 = 5.56 x width [76 s]
== the phase rule at the reach of each band
centre 2.67653e+11, width 2567: reach 11791 = 4.59 W; chirp over the window 5.5e+05 rad: random phases
centre 1.44177e+20, width 1409: reach 5901 = 4.19 W; chirp over the window 1.2e+15 rad: random phases
centre 1.37092e+21, width 1342: reach 6031 = 4.49 W; chirp over the window 1.1e+16 rad: random phases
centre 500671, width 1341: reach 6129 = 4.57 W; step 52.00 pi; nearest bad point (j = 52) 6130, reach -1 from it; bad points 118 apart; chirp over the window 3.84 rad
centre 1.00067e+06, width 1341: reach 6569 = 4.90 W; step 96.98 pi; nearest bad point (j = 97) 6567, reach +2 from it; bad points 68 apart; chirp over the window 6.68 rad
centre 501500, width 2999: reach 14531 = 4.85 W; step 21.97 pi; nearest bad point (j = 22) 14512, reach +19 from it; bad points 661 apart; chirp over the window 0.68 rad
centre 1.0015e+06, width 2999: reach 15193 = 5.07 W; step 41.97 pi; nearest bad point (j = 42) 15180, reach +13 from it; bad points 362 apart; chirp over the window 1.25 rad
centre 505000, width 9999: reach 40153 = 4.02 W; step 8.01 pi; nearest bad point (j = 8) 40187, reach -34 from it; bad points 5015 apart; chirp over the window 0.09 rad
centre 1.005e+06, width 9999: reach 37615 = 3.76 W; step 17.01 pi; nearest bad point (j = 17) 37635, reach -20 from it; bad points 2211 apart; chirp over the window 0.20 rad
6 bands with a small chirp; distance to the nearest bad point in half spacings: 0.009, 0.045, 0.057, 0.070, 0.013, 0.019. If the failures did not care about bad points these would be uniform on [0, 1]; Fisher's combined test: chi2 = 43.2 with 12 degrees of freedom, p = 2.1e-05
== cold start (nothing known below S - 6 scales), and the oracle start (true values from S - 6 scales to S)
1 band(s) at centre 1.37092e+21, total width 1342, 10000 zeros, wheel (2,), cold start at 1300: 79 primes decided; first wrong []; clean to 1900 = 1.42 x width [5 s]
1 band(s) at centre 1.37092e+21, total width 1342, 10000 zeros, wheel (2,), cold start at 2000: 75 primes decided; first wrong []; clean to 2600 = 1.94 x width [5 s]
1 band(s) at centre 1.37092e+21, total width 1342, 10000 zeros, wheel (2,), cold start at 2700: 70 primes decided; first wrong []; clean to 3300 = 2.46 x width [7 s]
1 band(s) at centre 1.37092e+21, total width 1342, 10000 zeros, wheel (2,), cold start at 3300: 3 primes decided; first wrong [3309]; reach 3309 = 2.47 x width [5 s]
3 band(s) at centre 2.67653e+11, 1.44177e+20, 1.37092e+21, total width 5318, 30000 zeros, wheel (2,), cold start at 5000: 114 primes decided; first wrong []; clean to 6000 = 1.13 x width [31 s]
3 band(s) at centre 2.67653e+11, 1.44177e+20, 1.37092e+21, total width 5318, 30000 zeros, wheel (2,), cold start at 8000: 110 primes decided; first wrong []; clean to 9000 = 1.69 x width [40 s]
3 band(s) at centre 2.67653e+11, 1.44177e+20, 1.37092e+21, total width 5318, 30000 zeros, wheel (2,), cold start at 11000: 2 primes decided; first wrong [11023]; reach 11023 = 2.07 x width [32 s]
1 band(s) at centre 1.37092e+21, total width 1342, 10000 zeros, wheel (2,), cold start at 3300: 71 primes decided; first wrong [3853]; reach 3853 = 2.87 x width [6 s]
1 band(s) at centre 1.37092e+21, total width 1342, 10000 zeros, wheel (2,), cold start at 4000: 72 primes decided; first wrong []; clean to 4600 = 3.43 x width [7 s]
1 band(s) at centre 1.37092e+21, total width 1342, 10000 zeros, wheel (2,), cold start at 5000: 25 primes decided; first wrong [5213]; reach 5213 = 3.89 x width [6 s]
== independent checks of the lists (sympy.isprime)
p3.txt: 1412 primes in [3, 11791), all prime: True; primes in range other than the wheel primes: 1412; complete: True
p4.txt: 775 primes in [3, 5901), all prime: True; primes in range other than the wheel primes: 775; complete: True
p5.txt: 785 primes in [3, 6031), all prime: True; primes in range other than the wheel primes: 785; complete: True
p345.txt: 2553 primes in [3, 22875), all prime: True; primes in range other than the wheel primes: 2553; complete: True
p5w6.txt: 1229 primes in [3, 10013), all prime: True; primes in range other than the wheel primes: 1229; complete: True
p5w210.txt: 1444 primes in [3, 12107), all prime: True; primes in range other than the wheel primes: 1444; complete: True
p345w210.txt: 5199 primes in [3, 50783), all prime: True; primes in range other than the wheel primes: 5199; complete: True
c5_1300.txt: 79 primes in [1300, 1901), all prime: True; primes in range other than the wheel primes: 79; complete: True
c5_2000.txt: 75 primes in [2000, 2601), all prime: True; primes in range other than the wheel primes: 75; complete: True
c5_2700.txt: 70 primes in [2700, 3301), all prime: True; primes in range other than the wheel primes: 70; complete: True
c345_5000.txt: 114 primes in [5000, 6001), all prime: True; primes in range other than the wheel primes: 114; complete: True
c345_8000.txt: 110 primes in [8000, 9001), all prime: True; primes in range other than the wheel primes: 110; complete: True
== the formula, zero side against prime side
x = 101: zero side -88.286872-3.678997j, prime side -88.286846-3.679005j, difference 2.7e-05
x = 997: zero side -46.823037+3.980941j, prime side -46.823040+3.980933j, difference 8.2e-06
x = 1000: zero side 5.450268+17.740316j, prime side 5.450271+17.740316j, difference 2.3e-06
x = 2003: zero side -73.466546+6.403458j, prime side -73.466543+6.403463j, difference 6.5e-06
x = 4999: zero side -7.792995+1.675283j, prime side -7.792998+1.675287j, difference 5.5e-06
largest difference 2.7e-05
== the primes check the zeros: Odlyzko's tables and our own zeros
odlyzko_zeros3.txt: 10000 zeros, height 2.6765e+11; x = 500..5000 step 100
  residual |S - P|: rms 1.1e-06, max 2.5e-06 (|P| rms 46.7)
  as zero errors: rms per zero <= 1.7e-09; common shift <= 1.9e-07
  zeros moved by rms 1e-7: residual rms 8.2e-05
  zeros moved by rms 1e-6: residual rms 7.9e-04
  zeros moved by shift 1e-7: residual rms 3.6e-05
odlyzko_zeros4.txt: 10000 zeros, height 1.4418e+20; x = 500..5000 step 100
  residual |S - P|: rms 6.2e-06, max 1.3e-05 (|P| rms 35.3)
  as zero errors: rms per zero <= 9.0e-09; common shift <= 9.4e-08
  zeros moved by rms 1e-7: residual rms 6.9e-05
  zeros moved by rms 1e-6: residual rms 6.6e-04
  zeros moved by shift 1e-7: residual rms 2.7e-05
odlyzko_zeros5.txt: 10000 zeros, height 1.3709e+21; x = 500..5000 step 100
  residual |S - P|: rms 6.7e-06, max 1.3e-05 (|P| rms 40.0)
  as zero errors: rms per zero <= 9.7e-09; common shift <= 2.3e-07
  zeros moved by rms 1e-7: residual rms 6.7e-05
  zeros moved by rms 1e-6: residual rms 6.4e-04
  zeros moved by shift 1e-7: residual rms 3.1e-05
band:2000:1000:: 952 zeros, height 2.5000e+03; x = 500..5000 step 100
  residual |S - P|: rms 6.5e-06, max 1.2e-05 (|P| rms 30.9)
  as zero errors: rms per zero <= 3.0e-08; common shift <= 1.9e-07
  zeros moved by rms 1e-7: residual rms 2.2e-05
  zeros moved by rms 1e-6: residual rms 2.3e-04
  zeros moved by shift 1e-7: residual rms 2.4e-05
band:20000:1000:: 1288 zeros, height 2.0500e+04; x = 500..5000 step 100
  residual |S - P|: rms 4.6e-06, max 1.1e-05 (|P| rms 33.2)
  as zero errors: rms per zero <= 1.8e-08; common shift <= 1.5e-07
  zeros moved by rms 1e-7: residual rms 2.6e-05
  zeros moved by rms 1e-6: residual rms 2.7e-04
  zeros moved by shift 1e-7: residual rms 2.6e-05
band:200000:1000:: 1651 zeros, height 2.0050e+05; x = 500..5000 step 100
  residual |S - P|: rms 7.3e-06, max 1.6e-05 (|P| rms 35.2)
  as zero errors: rms per zero <= 2.5e-08; common shift <= 4.0e-07
  zeros moved by rms 1e-7: residual rms 2.9e-05
  zeros moved by rms 1e-6: residual rms 3.2e-04
  zeros moved by shift 1e-7: residual rms 2.8e-05
band:1100000:1000:: 1922 zeros, height 1.1005e+06; x = 500..5000 step 100
  residual |S - P|: rms 8.0e-06, max 1.5e-05 (|P| rms 30.3)
  as zero errors: rms per zero <= 2.6e-08; common shift <= 3.7e-07
  zeros moved by rms 1e-7: residual rms 3.3e-05
  zeros moved by rms 1e-6: residual rms 3.0e-04
  zeros moved by shift 1e-7: residual rms 2.5e-05
== H4 at height 2.68e11
c = 2.0: X = 412787, 35004 prime powers [30 s]
   zeros 500-559: max 5.2e-04, rms 1.5e-04 (crossings 5.1e-02); worst zero's nearest gap 0.25 mean spacings
   zeros 2500-2559: max 3.7e-04, rms 1.0e-04 (crossings 4.6e-02); worst zero's nearest gap 0.48 mean spacings
   zeros 4500-4559: max 5.9e-04, rms 1.1e-04 (crossings 5.4e-02); worst zero's nearest gap 0.30 mean spacings
   zeros 6500-6559: max 5.4e-04, rms 1.4e-04 (crossings 5.3e-02); worst zero's nearest gap 0.14 mean spacings
   zeros 8500-8559: max 4.4e-04, rms 1.5e-04 (crossings 5.3e-02); worst zero's nearest gap 0.55 mean spacings
nearest gap 0-0.3 spacings: 16 zeros, median error 1.7e-04, max 5.4e-04
nearest gap 0.3-0.6 spacings: 94 zeros, median error 7.1e-05, max 5.9e-04
nearest gap 0.6-1.0 spacings: 142 zeros, median error 4.4e-05, max 4.3e-04
nearest gap 1.0-9 spacings: 48 zeros, median error 4.5e-05, max 1.5e-04
c = 1.6: X = 330230, 28584 prime powers [26 s]
   zeros 500-559: max 2.4e-04, rms 6.9e-05 (crossings 4.9e-02); worst zero's nearest gap 0.25 mean spacings
   zeros 2500-2559: max 7.6e-04, rms 2.3e-04 (crossings 4.8e-02); worst zero's nearest gap 0.25 mean spacings
   zeros 4500-4559: max 3.9e-03, rms 1.3e-03 (crossings 5.3e-02); worst zero's nearest gap 0.40 mean spacings
   zeros 6500-6559: max 1.5e-03, rms 4.4e-04 (crossings 5.3e-02); worst zero's nearest gap 0.14 mean spacings
   zeros 8500-8559: max 3.1e-04, rms 1.1e-04 (crossings 5.7e-02); worst zero's nearest gap 0.36 mean spacings
nearest gap 0-0.3 spacings: 16 zeros, median error 3.9e-04, max 1.5e-03
nearest gap 0.3-0.6 spacings: 94 zeros, median error 1.6e-04, max 3.9e-03
nearest gap 0.6-1.0 spacings: 142 zeros, median error 8.8e-05, max 1.5e-03
nearest gap 1.0-9 spacings: 48 zeros, median error 1.7e-04, max 1.2e-03
c = 1.3: X = 268312, 23661 prime powers [62 s]
   zeros 500-559: max 5.2e-03, rms 1.8e-03 (crossings 4.7e-02); worst zero's nearest gap 0.25 mean spacings
   zeros 2500-2559: max 1.5e-02, rms 4.5e-03 (crossings 5.0e-02); worst zero's nearest gap 0.25 mean spacings
   zeros 4500-4559: max 1.9e-02, rms 5.6e-03 (crossings 5.2e-02); worst zero's nearest gap 0.30 mean spacings
   zeros 6500-6559: max 3.1e-02, rms 1.2e-02 (crossings 5.4e-02); worst zero's nearest gap 0.14 mean spacings
   zeros 8500-8559: max 1.1e-02, rms 3.5e-03 (crossings 6.0e-02); worst zero's nearest gap 0.55 mean spacings
nearest gap 0-0.3 spacings: 16 zeros, median error 4.9e-03, max 3.1e-02
nearest gap 0.3-0.6 spacings: 94 zeros, median error 3.4e-03, max 2.8e-02
nearest gap 0.6-1.0 spacings: 142 zeros, median error 1.9e-03, max 1.9e-02
nearest gap 1.0-9 spacings: 48 zeros, median error 3.0e-03, max 1.1e-02
c = 1.0: X = 206394, 18639 prime powers [49 s]
   zeros 500-559: max 6.6e-03, rms 2.2e-03 (crossings 4.7e-02); worst zero's nearest gap 0.45 mean spacings
   zeros 2500-2559: max 5.7e-02, rms 1.9e-02 (crossings 5.2e-02); worst zero's nearest gap 0.43 mean spacings
   zeros 4500-4559: max 1.4e-02, rms 4.0e-03 (crossings 5.2e-02); worst zero's nearest gap 0.30 mean spacings
   zeros 6500-6559: max 3.4e-02, rms 1.3e-02 (crossings 5.7e-02); worst zero's nearest gap 0.49 mean spacings
   zeros 8500-8559: max 7.3e-02, rms 2.1e-02 (crossings 6.3e-02); worst zero's nearest gap 0.55 mean spacings
nearest gap 0-0.3 spacings: 16 zeros, median error 1.2e-02, max 3.3e-02
nearest gap 0.3-0.6 spacings: 94 zeros, median error 6.7e-03, max 7.3e-02
nearest gap 0.6-1.0 spacings: 142 zeros, median error 2.9e-03, max 3.0e-02
nearest gap 1.0-9 spacings: 48 zeros, median error 2.1e-03, max 1.6e-02
```

Output of band_validate_m.sh:

```
2000000 zeros, the last at 1131944.471825 [35.2 s]
== with 8 decision histories
1 band(s) at centre 2.67653e+11, total width 2567, 10000 zeros, wheel (2,), M = 8: 1688 primes decided; first wrong [14399, 14401, 14409]; reach 14399 = 5.61 x width [31 s]
1 band(s) at centre 1.44177e+20, total width 1409, 10000 zeros, wheel (2,), M = 8: 969 primes decided; first wrong [7633, 7635, 7639]; reach 7633 = 5.42 x width [19 s]
1 band(s) at centre 1.37092e+21, total width 1342, 10000 zeros, wheel (2,), M = 8: 968 primes decided; first wrong [7615, 7617, 7621]; reach 7615 = 5.68 x width [18 s]
3 band(s) at centre 2.67653e+11, 1.44177e+20, 1.37092e+21, total width 5318, 30000 zeros, wheel (2,), M = 8: 3325 primes decided; first wrong [30831, 30833, 30841]; reach 30831 = 5.80 x width [193 s]
3 band(s) at centre 2.67653e+11, 1.44177e+20, 1.37092e+21, total width 5318, 30000 zeros, wheel (2, 3, 5, 7), M = 8: 6478 primes decided; first wrong [64853, 64859, 64867]; reach 64853 = 12.19 x width [242 s]
1 band(s) at centre 0, total width 1341, 933 zeros, wheel (2,), M = 8: 1084 primes decided; first wrong [8699, 8701, 8703]; reach 8699 = 6.49 x width [12 s]
1 band(s) at centre 500671, total width 1341, 2410 zeros, wheel (2,), M = 8: 954 primes decided; first wrong [7527, 7529, 7531]; reach 7527 = 5.61 x width [13 s]
1 band(s) at centre 0, total width 9999, 10142 zeros, wheel (2,), M = 8: 5724 primes decided; first wrong [56467, 56469, 56471]; reach 56467 = 5.65 x width [109 s]
1 band(s) at centre 505000, total width 9999, 17975 zeros, wheel (2,), M = 8: 4743 primes decided; first wrong [45775, 45777, 45779]; reach 45775 = 4.58 x width [197 s]
== cold start with 8 decision histories
1 band(s) at centre 1.37092e+21, total width 1342, 10000 zeros, wheel (2,), cold start at 3300, M = 8: 76 primes decided; first wrong []; clean to 3900 = 2.91 x width [5 s]
1 band(s) at centre 1.37092e+21, total width 1342, 10000 zeros, wheel (2,), cold start at 4000, M = 8: 72 primes decided; first wrong []; clean to 4600 = 3.43 x width [6 s]
1 band(s) at centre 1.37092e+21, total width 1342, 10000 zeros, wheel (2,), cold start at 4600, M = 8: 0 primes decided; first wrong [4603]; reach 4603 = 3.43 x width [4 s]
1 band(s) at centre 1.37092e+21, total width 1342, 10000 zeros, wheel (2,), cold start at 5200, M = 8: 1 primes decided; first wrong [5201]; reach 5201 = 3.88 x width [4 s]
3 band(s) at centre 2.67653e+11, 1.44177e+20, 1.37092e+21, total width 5318, 30000 zeros, wheel (2,), cold start at 9500, M = 8: 1 primes decided; first wrong [9501]; reach 9501 = 1.79 x width [20 s]
3 band(s) at centre 2.67653e+11, 1.44177e+20, 1.37092e+21, total width 5318, 30000 zeros, wheel (2,), cold start at 10000, M = 8: 106 primes decided; first wrong []; clean to 11000 = 2.07 x width [29 s]
3 band(s) at centre 2.67653e+11, 1.44177e+20, 1.37092e+21, total width 5318, 30000 zeros, wheel (2,), cold start at 11000, M = 8: 0 primes decided; first wrong [11003]; reach 11003 = 2.07 x width [23 s]
== two bands of width 5,000 whose bad points coincide (c2 = 2 c1) or not
2 band(s) at centre 502500, 1.0025e+06, total width 9998, 18516 zeros, wheel (2,): 4192 primes decided; first wrong [39869, 39871, 39873]; reach 39869 = 3.99 x width [139 s]
2 band(s) at centre 502500, 702500, total width 9998, 18233 zeros, wheel (2,): 4462 primes decided; first wrong [42675, 42679, 42681]; reach 42675 = 4.27 x width [162 s]
centre 505000, width 9999: reach 45775 = 4.58 W; step 7.02 pi; nearest bad point (j = 7) 45928, reach -153 from it; bad points 6518 apart; chirp over the window 0.07 rad
2 bands, width 9998: reach 39869 = 3.99 W; centre 502500: step 8.02 pi, 119 from the nearest bad point (j = 8); centre 1.0025e+06: step 16.01 pi, 19 from the nearest bad point (j = 16)
2 bands, width 9998: reach 42675 = 4.27 W; centre 502500: step 7.50 pi, 2687 from the nearest bad point (j = 8); centre 702500: step 10.48 pi, 2018 from the nearest bad point (j = 11)
1 bands with a small chirp; distance to the nearest bad point in half spacings: 0.047. If the failures did not care about bad points these would be uniform on [0, 1]; Fisher's combined test: chi2 = 6.1 with 2 degrees of freedom, p = 4.7e-02
== independent checks of the lists (sympy.isprime)
m3.txt: 1685 primes in [3, 14399), all prime: True; primes in range other than the wheel primes: 1685; complete: True
m4.txt: 967 primes in [3, 7633), all prime: True; primes in range other than the wheel primes: 967; complete: True
m5.txt: 966 primes in [3, 7615), all prime: True; primes in range other than the wheel primes: 966; complete: True
m345.txt: 3322 primes in [3, 30831), all prime: True; primes in range other than the wheel primes: 3322; complete: True
m345w210.txt: 6476 primes in [3, 64853), all prime: True; primes in range other than the wheel primes: 6476; complete: True
mc5_3300.txt: 76 primes in [3300, 3901), all prime: True; primes in range other than the wheel primes: 76; complete: True
mc5_4000.txt: 72 primes in [4000, 4601), all prime: True; primes in range other than the wheel primes: 72; complete: True
mc345_10000.txt: 106 primes in [10000, 11001), all prime: True; primes in range other than the wheel primes: 106; complete: True
```

## 10. Hypotheses and what would settle them

| Hypothesis | Evidence so far | Confirms it | Refutes it | What is needed |
|---|---|---|---|---|
| **H1** Blur law: the height needed is about 0.2·πn ÷ candidate spacing (a blur of about 5 spacings) | spacings 2, 4, 4.4, 24 and 840 (section 4) | a new spacing on the predicted line, e.g. the 32 characters mod 120 (spacing 120 → about 0.0017·πn) | a spacing that needs clearly more or less | zeros of the 32 L-functions mod 120 (built off the page before), and the class decoder of 136 with 32 sets |
| **H2** Divisibility law 0.23·π·ρ·n | k = 2 … 25 up to 300,000 | the same factor for k > 25 and for n of 10⁶–10⁷ | a drift of the factor with k or n | more zeros: n = 10⁷ at k = 25 needs height about 8.7·10⁵ (the page has 1.13·10⁶) |
| **H3** The zeros carry the primes down to about 0.02–0.025·πn; the gap from 0.045 is a search problem | oracle tests; truth fits better in 9 of 10 failures; M-algorithm 19 → 0 errors at 0.035 | a decoder that is clean at 0.035 or lower on 100 random circles to 300,000 | many failures where the wrong pattern has the lower whitened cost (then the information is missing), or any decoder clean below 0.02 (then the limit is lower) | the M-algorithm on the page decoder with delayed decisions; the oracle test on more stretches |
| **H4** Backwards: the prime powers up to c·√(T/2π), c ≈ 1.6–2, fix every zero near height T; the threshold is a count | X = 30, 100, 300, heights to 5.7·10⁵; near the 10¹²-th zero (2.68·10¹¹, Odlyzko's table) c = 2.0 fixes 5 of 5 windows of 60 zeros to < 6·10⁻⁴, c = 1.6 3 of 5, c ≤ 1.3 none (section 13) | many windows per (X, c); the threshold moves with the fade as the count predicts (fading from 0.9 instead of 6/7 should lower c) | success well below c = 1 (better than the count), or c growing without bound with T | the zeros near 10²¹ need the primes to 7.7·10⁹: a segmented prime sum in C; or a test of the oversampling law 2 ln X / ln(T/2π) against fixed c |
| **H5** Close pairs of zeros and twin positions are the same obstacle seen from both sides | the only zeros above 10⁻³ in the chain are two close pairs; near the 10¹²-th zero the median error is 1.7·10⁻⁴ for a neighbour within 0.3 spacings against 4.4·10⁻⁵ beyond 0.6 (section 13) | the forward errors at spacing 2 and the backward errors follow one curve in the normalised gap | close pairs recovered as well as the rest | the forward error against the gap to the next candidate, on the same scale |
| **H6** Bandwidth and phase law: a band of width W reads the primes up to about 4–4.6·W (5.4–5.8·W with 8 histories) when its phases are random (c ≫ n²); just above the bottom it fails at the bad points 2c/(πj); bands add | Odlyzko's 3 bands at 2.7·10¹¹–1.4·10²¹; 6 own bands at 5·10⁵–10⁶ all failing at a bad point (Fisher p = 2·10⁻⁵); sums of bands within 4 % of the sum of reaches (section 12) | more pairs of bands with non-aligned bad points (one pair so far: 42,675 against 39,869 aligned); longer random bands (10⁵ zeros) on the same factor | a failure far from every bad point in a band with a small chirp, or a reach that depends on the height in the random regime | a long band far up, e.g. from the LMFDB; a test of combined bands whose bad points differ |
| **H7** Random access stops at 2–3·W mainly because the cold start avalanches | cold starts clean at S ≤ 2.0·W (one band) and 1.5·W (three bands), total failure beyond (section 12). Oracle: with the true values of the 6 scales before S given (nothing below), the band near 10²² is clean from 4,000 to 4,600 (3.4·W) and fails at 3,853 (from 3,300) and 5,213 (from 5,000); from 3 upwards it is clean to 6,031. With 8 histories the cold start is clean from 3,300 and 4,000 (3.0·W) | a joint start (the first 30–40 candidates decided together) that is clean near 3.4·W | the oracle start no better than the cold one | the joint start; more oracle stretches |

**The Riemann Hypothesis itself** is not tested by any of this. The zero finder assumes the zeros are on the line, and checks it in the usual way: it counts sign changes of Z(t) and compares with N(T). Both agree up to height 1.13·10⁶; RH is known by others to hold much higher. Computation can only refute RH, by a height where the sign changes fall short of N(T), and that would have to be far above anything a page can reach. A confirmation needs a proof, not a computation.

## 11. Open

- The M-algorithm on the page (section 8).
- A prior for clusters from the Hardy–Littlewood singular series. The admissible tuples of the bounded-gaps work (Zhang, Annals 179, 2014) are exactly the clusters where the decoders fail; this is untested.
- A second round of the chain: the primes below 70,000 would fix the zeros up to about 2π(70000/1.6)² ≈ 1.2·10¹⁰, beyond what can be computed here.
- Bands on the page: a slider for the height of the band, with Odlyzko's three tables (section 12).
- Bands of L-function zeros: with the class formulas 135–137 a band of each L(s, χ) mod q should give the same factor q.
- A better cold start (H7). The M-algorithm (8 histories) is in the band decoder; on the page it would lift 134 and 139 the same way.
- Bands that avoid each other's bad points: one pair tested (section 12); a band split in two with a chosen ratio of centres.

## 12. Every band of zeros carries the primes

Landau's formula (1911) says that for fixed x > 1, Σ_{0<γ≤T} x^ρ = −(T/2π) Λ(x) + O(log T). So every band of zeros has a spike at every prime power, not only the band at the bottom. Odlyzko (1987) saw the peaks at ln pᵏ in the zeros near the 10¹²-th zero. Here a band is read exactly, prime by prime.

**The formula.** Weight the zeros of a band [c − W/2, c + W/2] with a Tukey window w: flat, with cosine fades of length W/7 at both ends. Let Ŵ be its Fourier transform (real, even). Then

```
r(x) = −(2π √x / Ŵ(0)) Σ_γ w(γ) e^{i(γ − c) ln x}  =  Σ_m Λ(m) √(x/m) e^{−ic ln m} Ŵ(ln(x/m)) / Ŵ(0)
```

Derivation: N(t) = θ(t)/π + 1 + S(t), and S′(t) = −(1/π) Σ_m Λ(m) m^{−1/2} cos(t ln m) as a distribution (the explicit formula). Multiply by w(t) e^{it ln x} and integrate. The half of the cosine with e^{+it ln m} oscillates away; the θ′ term is a smooth function under a window without low frequencies and vanishes. What is left is −(1/2π) Σ_m Λ(m) m^{−1/2} ∫ w(t) e^{it ln(x/m)} dt.

So, after demodulation by e^{−ic ln x}, a band is the lowpass view of section 1 in which every prime power m carries the known phase e^{−ic ln m}. The lowpass is the case c = 0 with the mirror zeros: a band of width 2T. Check (computed): with every prime power below 20x + 1000, the zero side of Odlyzko's band near the 10²²-th zero, summed in mpmath from the decimal strings of the table, equals the prime side, with the kernel integrated numerically, to at most 2.7·10⁻⁵ at x = 101, 997, 1000, 2003, 4999, against sums of size 7 to 90 ([band_check.py](writeup/code/riemann/band_check.py)). Cutting the prime side at |ln(x/m)| < 0.3 leaves 2·10⁻³: the far prime powers matter.

**The decoder** ([band_decoder.py](writeup/code/riemann/band_decoder.py)) is the exact window decoder of section 2 on complex samples:
- the next 12 candidates searched exhaustively;
- the candidates after them, over 12 envelope scales 2πn/W, modelled as coins;
- samples every third of the envelope scale, from 4.4 scales left of n to 3 scales right of the window;
- noise 0.01 per sample times the sampling step. Without that factor the coarser sampling costs reach (6,043 → 5,483 in a test); with it, the reach is the same as with a sample at every integer, at a quarter of the time. An exact window of 16 candidates instead of 12 was worse (6,031 → 5,483);
- with `--M 8`, the M-algorithm of section 8: the 8 best decision histories by summed window cost, a decision committed when it is 24 candidates old. The covariance, its Cholesky factor and the QR do not depend on the history, only the residual does, so 8 histories cost about twice the time of one.

Every prime below the first wrong decision is read, and each list is checked against `sympy.isprime` (section 9). Odd candidates unless a wheel is named; "own zeros" are those of zeros.js:

| Zeros | Centre c | Zeros used | Width W | First wrong | ÷ W | With 8 histories | ÷ W |
|---|---|---|---|---|---|---|---|
| Odlyzko, numbers 10¹² + 1 … 10¹² + 10⁴ | 2.68·10¹¹ | 10,000 | 2,567 | 11,791 | 4.59 | 14,399 | 5.61 |
| Odlyzko, numbers 10²¹ + 1 … 10²¹ + 10⁴ | 1.44·10²⁰ | 10,000 | 1,409 | 5,901 | 4.19 | 7,633 | 5.42 |
| Odlyzko, numbers 10²² + 1 … 10²² + 10⁴ | 1.37·10²¹ | 10,000 | 1,342 | 6,031 | 4.50 | 7,615 | 5.68 |
| the three together | | 30,000 | 5,318 | 22,875 | 4.30 | 30,831 | 5.80 |
| the three together, coprime to 210 | | 30,000 | 5,318 | 50,783 | 9.55 | 64,853 | 12.19 |
| band near 10²², coprime to 6 | 1.37·10²¹ | 10,000 | 1,342 | 10,013 | 7.46 | | |
| band near 10²², coprime to 210 | 1.37·10²¹ | 10,000 | 1,342 | 12,107 | 9.02 | | |
| lowpass, own zeros up to 1,342 | 0 | 933 | 1,341 | 8,699 | 6.49 | 8,699 | 6.49 |
| lowpass, own zeros up to 3,000 | 0 | 2,469 | 2,999 | 19,853 | 6.62 | | |
| lowpass, own zeros up to 10,000 | 0 | 10,142 | 9,999 | 56,467 | 5.65 | 56,467 | 5.65 |
| own zeros in (500,000, 501,342] | 500,671 | 2,410 | 1,341 | 6,129 | 4.57 | 7,527 | 5.61 |
| own zeros in (1,000,000, 1,001,342] | 1,000,671 | 2,558 | 1,341 | 6,569 | 4.90 | | |
| own zeros in (500,000, 503,000] | 501,500 | 5,389 | 2,999 | 14,531 | 4.85 | | |
| own zeros in (1,000,000, 1,003,000] | 1,001,500 | 5,719 | 2,999 | 15,193 | 5.07 | | |
| own zeros in (500,000, 510,000] | 505,000 | 17,975 | 9,999 | 40,153 | 4.02 | 45,775 | 4.58 |
| own zeros in (1,000,000, 1,010,000] | 1,005,000 | 19,070 | 9,999 | 37,615 | 3.76 | | |
| lowpass up to 1,342 and the band near 10²² | | 10,933 | 2,683 | 14,921 | 5.56 | | |

**Laws (computed):**
- **The height itself does not count.** The three Odlyzko bands, at heights 2.7·10¹¹ to 1.4·10²¹, read odd candidates to 4.2–4.6·W (5.4–5.8·W with 8 histories), although near the 10²²-th zero a band of width 1,342 holds 10,000 zeros against 933 at the bottom.
- **The phase rule.** Around n the phase of m is c ln(m/n) ≈ c(m − n)/n − c(m − n)²/(2n²). The step between odd neighbours is 2c/n (mod 2π), and over the window of 12 candidates the chirp is c·24²/(2n²). Where the chirp is small and the step is a multiple of π, the phases of a window are real up to one common factor: the real and the imaginary part carry the same signal, and the band gives only half its information there. The bad points are n_j = 2c/(πj). Every band of our own zeros (c = 5·10⁵ – 10⁶, chirp 0.1–6.7 rad) fails at a bad point: at −1, +2, +19, +13, −34 and −20 from it, with bad points 68 to 5,015 apart ([band_phase.py](writeup/code/riemann/band_phase.py)). If the failures did not care about bad points, the distances in half spacings would be uniform on [0, 1]; Fisher's combined test gives p = 2·10⁻⁵. The rule was stated after the first runs. One run was predicted before it was made, with an earlier version of the decoder (far coins counted at the end of the run): the band (1,000,000, 1,010,000] failed 5 before the bad point j = 19 (33,674), after passing j = 20 and 21. The final decoder passes j = 19 and fails at j = 17. Which bad point stops the decoder is not predicted. Odlyzko's bands have chirps of 10⁵–10¹⁶ rad: random phases, no bad points. With 8 histories the band at c = 500,671 passes its bad point (7,527 instead of 6,129), and the band at c = 505,000 passes j = 8 and fails 153 before j = 7 (45,775 instead of 40,153): the histories carry the decoder past some bad points, not past the rule.
- **A band against the lowpass of the same width.** Width 1,341: lowpass 8,699, bands 6,031–6,569 (0.69–0.76); width 3,000: 19,853 against 14,531–15,193 (0.73–0.77); width 10⁴: 56,467 against 37,615–40,153 (0.67–0.71, both at a bad point). With 8 histories the random bands reach 5.4–5.8·W against 6.49 for the lowpass (0.84–0.89), and the band of width 10⁴ 45,775 against 56,467 (0.81); the lowpass gains nothing from more histories (widths 1,341 and 10⁴). So most of the gap is search, not information: both give W/π real numbers per unit of ln x (section 4).
- **Bands add.** The three Odlyzko bands together read to 22,875, against 23,723 for the sum of their separate reaches; the lowpass and the band near 10²² together to 14,921, against 14,730.
- **Bands whose bad points differ avoid them.** Two bands of width 5,000 with c₂ ≈ 2c₁, (500,000, 505,000] and (1,000,000, 1,005,000], share every bad point of the first (j and 2j) and fail at one: 39,869 = 3.99·W, 119 and 19 from the bad points j = 8 and j = 16. With c₂/c₁ ≈ 1.4, (500,000, 505,000] and (700,000, 705,000], they read to 42,675 = 4.27·W and fail away from every bad point (steps 7.50π and 10.48π, 2,687 and 2,018 from the nearest), as the random bands far up do. One band of width 10⁴ at the same heights: 37,615–40,153.
- **Divisibility facts work as on the page:** coprime to 6 gives ×1.66, coprime to 210 ×2.01 (134 → 139 gave ×1.67).
- **The reach per width falls with the width,** for the lowpass too: 6.49, 6.62, 5.65 at widths 1,341, 3,000, 10⁴, with one history or eight.
- **Random access.** With `--cold S` nothing is known below S − 6 envelope scales; the numbers there are coins too.
  - One band (near 10²²): clean from S = 1,300, 2,000 and 2,700 (each to S + 600), fails from 3,300 = 2.5·W. With 8 histories clean from 3,300 and 4,000 (3.0·W), fails from 4,600 and 5,200.
  - The three bands: clean from S = 5,000 and 8,000 (each to S + 1,000), fail from 11,000 = 2.1·W. With 8 histories clean from 10,000, fail from 9,500 and 11,000.
  - So the 110 primes in [8,000, 9,000] come from the zeros near the 10¹²-, 10²¹- and 10²²-th zero without any prime below 7,775.
  - A failure is total: the first decisions go wrong and the echo of the errors spoils the rest. Oracle: with the true values of the 6 scales before S given (nothing below), the band near 10²² is clean from 4,000 to 4,600 (3.4·W) and fails at 3,853 (from 3,300) and 5,213 (from 5,000).

**What it means.** The explicit formula stores the primes in every part of the spectrum. A band of zeros at any height gives the primes up to about 4–4.6·W, 5.4–5.8·W with a better search, and 9–12·W with the wheel mod 210; only the width and the phase pattern count. Just above the bottom (c below about n²/3, a chirp below 100 rad over the window) the phases form a slow carrier, and the band loses half its information wherever the carrier steps by a multiple of π; far up the phases are random. For zeros per prime the bottom stays cheapest, because the zeros get denser with height (ln(T/2π)/2π per unit). The chain of section 6 does not close at height: the band near T is fixed by the prime powers up to about 2√(T/2π) (section 13) and gives back the primes up to about 4.4·W, so it returns more than it needs only for W > 0.45√(T/2π), near the 10¹²-th zero a band of width 94,000 (366,000 zeros).

## 13. The primes check the zeros, at heights the page cannot reach

**Accuracy of published zeros.** The formula of section 12 holds to the precision of the zeros. An error ε_j in zero j changes the zero side by i w_j e^{iγ_j ln x} ε_j ln x, so random errors of rms ε give a residual of about ε ln x (Σ w_j²)^{1/2}, and a common shift d gives about d ln x times the size of the sum. Over x = 500, 600, …, 5,000 (46 values; prime powers below 20x + 1000; [zeros_vs_primes.py](writeup/code/riemann/zeros_vs_primes.py)):

| Table | Stated accuracy | Residual rms (max) | rms error per zero ≤ | Common shift ≤ |
|---|---|---|---|---|
| numbers 10¹² + 1 … (10 decimals) | 10⁻⁸, guaranteed | 1.1·10⁻⁶ (2.5·10⁻⁶) | 1.7·10⁻⁹ | 1.9·10⁻⁷ |
| numbers 10²¹ + 1 … (8 decimals) | "probably 10⁻⁶" | 6.2·10⁻⁶ (1.3·10⁻⁵) | 9.0·10⁻⁹ | 9.4·10⁻⁸ |
| numbers 10²² + 1 … (8 decimals) | "probably 10⁻⁶" | 6.7·10⁻⁶ (1.3·10⁻⁵) | 9.7·10⁻⁹ | 2.3·10⁻⁷ |

Our own zeros (zeros.js), in bands of width 1,000 starting at heights 2,000, 20,000, 200,000 and 1,100,000: residual rms 4.6–8.0·10⁻⁶, so an rms error per zero of at most 1.8–3.0·10⁻⁸, and a common shift of at most 1.5–4.0·10⁻⁷. Until now only 43 of the 2,000,000 zeros were compared with mpmath (section 5).

Calibration: moving every zero by random errors of rms 10⁻⁷ raises the residual of Odlyzko's bands to 6.7–8.2·10⁻⁵, rms 10⁻⁶ to 6.4–7.9·10⁻⁴, and a common shift of 10⁻⁷ to 2.7–3.6·10⁻⁵, linear in the error as the formula says. So the prime powers below 101,000 confirm the zeros near the 10²¹-th and 10²²-th zero to about 10⁻⁸ rms, a hundred times better than the tables claim; rounding to 8 decimals alone accounts for 2.9·10⁻⁹. The same residual shows that no zero is missing from the three tables and none is spurious: one missing zero would add w(γ) to the zero side, about 1 except within 0.4 of the ends of a band. What it does not show: a single zero off by up to about 1.5·10⁻⁶ (one error ε adds at most ε ln x ≤ 8.5ε).

**Backwards at height 2.68·10¹¹ (H4).** The fit of section 6 with the prime powers up to X = c√(T/2π) (here √(T/2π) = 206,394), in five windows of 60 zeros spread over the 10⁴ zeros near the 10¹²-th zero, scored against Odlyzko's table ([zeros_high.py](writeup/code/riemann/zeros_high.py)). The phases B ln n mod 2π (B = 267,653,395,647) are taken in mpmath:

| c | X | Prime powers | Max error per window | rms per window |
|---|---|---|---|---|
| 2.0 | 412,787 | 35,004 | 3.7–5.9·10⁻⁴ | 1.0–1.5·10⁻⁴ |
| 1.6 | 330,230 | 28,584 | 2.4·10⁻⁴ – 3.9·10⁻³ (3 of 5 below 10⁻³) | 0.7–13·10⁻⁴ |
| 1.3 | 268,312 | 23,661 | 5.2·10⁻³ – 3.1·10⁻² | 1.8–12·10⁻³ |
| 1.0 | 206,394 | 18,639 | 6.6·10⁻³ – 7.3·10⁻² | 2.2–21·10⁻³ |

The count X ≈ c√(T/2π) holds 5·10⁵ times higher than in section 6, with c ≈ 2 instead of 1.6 for every window. In the count of section 6 the oversampling 2 ln X / ln(T/2π) = 1 + 2 ln c / ln(T/2π) falls with height at fixed c, which may explain the drift.

**Close pairs (H5).** At c = 2.0 the median error of the 300 zeros grows as the gap to the nearest neighbour shrinks: 1.7·10⁻⁴ below 0.3 mean spacings (16 zeros), 7.1·10⁻⁵ at 0.3–0.6 (94), 4.4·10⁻⁵ at 0.6–1.0 (142), 4.5·10⁻⁵ above (48). Close pairs are the weak spots of this direction, as twin positions are of the other.

## References

- B. Riemann, Über die Anzahl der Primzahlen unter einer gegebenen Grösse (1859); H. von Mangoldt (1895) for ψ(x).
- E. Landau, Über die Nullstellen der Zetafunktion, Math. Ann. 71 (1911), 548–564: Σ_{0<γ≤T} x^ρ = −(T/2π) Λ(x) + O(log T); S. M. Gonek, An explicit formula of Landau and its applications to the theory of the zeta-function, Contemp. Math. 143 (1993), for the uniform version (section 12).
- A. M. Odlyzko, On the distribution of spacings between zeros of the zeta function, Math. Comp. 48 (1987), 273–308: the peaks at ln pᵏ in the zeros near the 10¹²-th zero. His tables of zeros (https://www-users.cse.umn.edu/~odlyzko/zeta_tables/, files zeros3, zeros4, zeros5) are copied in [writeup/data/riemann/](writeup/data/riemann/) with their source line.
- A. Weil's explicit formula, used here for the kernel of N_X.
- The Riemann–Siegel formula with the correction terms C0, C1, C2; the coefficients here were fitted with mpmath and agree with the published ones for C1.
- Y. Zhang, Discrete mean estimates and the Landau–Siegel zero, arXiv:2211.02515 (v1, 4 Nov 2022). Prop. 2.2: under a Siegel zero the zeros in his box have gaps π/log P, the resolution scale π/ln X of section 6. Lemma 3.3: the large sieve over characters, the reason more characters act as near-independent channels (135–137). Eq. (2.23)–(2.25): the mollifier weights tested in section 6.
- Y. Zhang, Bounded gaps between primes, Ann. of Math. 179 (2014), 1121–1174; only the abstract was read.
- Details and dates of each step are in [IDEAS.txt](IDEAS.txt) and [FINDINGS.txt](FINDINGS.txt) (entries of 10–12 November 2026).
