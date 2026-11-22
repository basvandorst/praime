/* Quartics F(y) = (a y^4 + b y^3 + c y^2 + d y + e)/k, k in {1,2,4}, 0 <= b < 4a (shift-normalised), c in [CL, CH],
   |d| <= DR, |e| <= ER. Conditions: N(y) = numerator = k mod 2k for every y (F integer and odd); N has no root mod every
   odd prime p <= PS. e enumerated by residue classes mod 2k*15015; 17..PS by table. Longest run of consecutive y in
   [-Y, Y] with |F(y)| prime is printed when >= MINRUN. usage: qq a k CL CH DR ER PS Y MINRUN */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef unsigned __int128 u128;
static uint8_t *bits; static uint64_t LIM;
static uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m){ return (u128)a*b % m; }
static uint64_t powmod(uint64_t a, uint64_t e, uint64_t m){ uint64_t r=1; a%=m; while(e){ if(e&1) r=mulmod(r,a,m); a=mulmod(a,a,m); e>>=1;} return r; }
static int mr(uint64_t n){ if(n<2) return 0; static const uint64_t P[]={2,3,5,7,11,13,17,19,23,29,31,37};
  for(int i=0;i<12;i++){ if(n%P[i]==0) return n==P[i]; }
  uint64_t d=n-1; int s=0; while(!(d&1)){d>>=1;s++;}
  for(int i=0;i<12;i++){ uint64_t x=powmod(P[i],d,n); if(x==1||x==n-1) continue; int ok=0; for(int r=1;r<s;r++){ x=mulmod(x,x,n); if(x==n-1){ok=1;break;} } if(!ok) return 0; }
  return 1; }
static inline int isp(int64_t v){ uint64_t n = v<0 ? (uint64_t)(-v) : (uint64_t)v; if(n<LIM){ if(n==2) return 1; if(!(n&1)) return 0; return (bits[n>>4]>>((n>>1)&7))&1; } return mr(n); }
static inline int64_t md(int64_t x, int64_t p){ x%=p; return x<0?x+p:x; }
int main(int argc, char **argv){
  int64_t a=atoll(argv[1]), k=atoll(argv[2]), CL=atoll(argv[3]), CH=atoll(argv[4]), DR=atoll(argv[5]), ER=atoll(argv[6]);
  int PS=atoi(argv[7]), Y=atoi(argv[8]), MINRUN=atoi(argv[9]); int64_t V=argc>10?atoll(argv[10]):(int64_t)4e18; int H=argc>11?atoi(argv[11]):18;
  LIM=1ULL<<29; bits=malloc(LIM>>4); memset(bits,0xff,LIM>>4); bits[0]&=~2;
  for(uint64_t i=3;i*i<LIM;i+=2) if((bits[i>>4]>>((i>>1)&7))&1) for(uint64_t j=i*i;j<LIM;j+=2*i) bits[j>>4]&=~(1<<((j>>1)&7));
  int pr[32],np=0; for(int p=3;p<=PS;p+=2){ int q=1; for(int d=3;d*d<=p;d+=2) if(p%d==0) q=0; if(q) pr[np++]=p; }
  int nsm=0; while(nsm<np && pr[nsm]<=13) nsm++;
  int64_t K2=2*k, M=K2*15015; int64_t mods[8]; mods[0]=K2; for(int i=0;i<nsm;i++) mods[i+1]=pr[i];
  int64_t crt[8]; for(int i=0;i<=nsm;i++){ int64_t p=mods[i], Mi=M/p, inv=1; while((Mi%p)*inv%p!=1) inv++; crt[i]=Mi*inv%M; }
  static uint8_t allow[32][64]; int64_t *cls=malloc(sizeof(int64_t)*M), *buf=malloc(sizeof(int64_t)*M);
  int W=2*Y+1; uint8_t *ok=malloc(W); long long cand=0, tested=0; long long hist[128]={0};
  for(int64_t b=0;b<4*a;b++) for(int64_t c=CL;c<=CH;c++) for(int64_t d=-DR;d<=DR;d++){
    if(8*a*a*d-4*a*b*c+b*b*b==0) continue;   /* symmetric: every value twice */
    /* size: |N(y)| <= V*k on the core window |y| <= H gives an interval for e */
    int64_t gmax=-(int64_t)4e18, gmin=(int64_t)4e18;
    for(int64_t y=-H;y<=H;y++){ int64_t g=(((a*y+b)*y+c)*y+d)*y; if(g>gmax) gmax=g; if(g<gmin) gmin=g; }
    int64_t elo=-V*k-gmin, ehi=V*k-gmax; if(elo<-ER) elo=-ER; if(ehi>ER) ehi=ER; if(elo>ehi) continue;
    /* mod 2k: Q(y) mod 2k must be constant over y */
    int64_t q0=md(((a*0+b)*0+c)*0*0+d*0,K2); int dead=0;
    for(int64_t y=1;y<K2;y++){ int64_t q=md((((a*y+b)*y+c)*y+d)*y,K2); if(q!=q0){dead=1;break;} }
    if(dead) continue;
    int64_t e2=md(k-q0,K2);
    for(int i=0;i<np;i++){ int64_t p=pr[i]; memset(allow[i],1,p); int64_t A=md(a,p),B=md(b,p),C=md(c,p),D=md(d,p);
      for(int64_t y=0;y<p;y++){ int64_t v=((((A*y+B)%p*y+C)%p*y+D)%p*y)%p; allow[i][(p-v)%p]=0; }
      int any=0; for(int r=0;r<p;r++) any|=allow[i][r]; if(!any){dead=1;break;} }
    if(dead) continue;
    int nc=1; cls[0]=crt[0]*e2%M;
    for(int i=0;i<nsm;i++){ int p=pr[i]; int nb=0;
      for(int j=0;j<nc;j++) for(int r=0;r<p;r++) if(allow[i][r]) buf[nb++]=(cls[j]+crt[i+1]*r)%M;
      memcpy(cls,buf,sizeof(int64_t)*nb); nc=nb; }
    for(int j=0;j<nc;j++) for(int64_t e=cls[j]+((elo-cls[j])/M-1)*M; e<=ehi; e+=M){
      if(e<elo) continue; cand++;
      int pass=1; for(int i=nsm;i<np;i++) if(!allow[i][md(e,pr[i])]){pass=0;break;}
      if(!pass) continue; tested++;
      int best=0,bs=0,cur=0;
      for(int y=-Y;y<=Y;y++){ int64_t N=(((a*y+b)*y+c)*y+d)*y+e; int t=isp(N/k); cur=t?cur+1:0; if(cur>best){best=cur;bs=y-cur+1;} }
      hist[best<127?best:127]++;
      if(best>=MINRUN) printf("%lld %lld %lld %lld %lld %lld %d %d\n",(long long)a,(long long)k,(long long)b,(long long)c,(long long)d,(long long)e,bs,best);
    }
  }
  fprintf(stderr,"a=%lld k=%lld cand %lld tested %lld | runs: ",(long long)a,(long long)k,cand,tested);
  for(int L=30;L<128;L++) if(hist[L]) fprintf(stderr,"%d:%lld ",L,hist[L]); fprintf(stderr,"\n"); return 0; }
