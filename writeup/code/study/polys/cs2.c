/* Cubics F(y) = a y^3 + b y^2 + c y + d, 0 <= b < 3a, c in [CL*a, CH*a], |d| <= DR*a. Filter: F has no root mod every
   prime p <= PS. d enumerated by residue classes mod 30030 (product over p <= 13 of the allowed sets), p = 17..PS by table.
   Longest run of consecutive y in [-Y, Y] with |F(y)| prime (no distinctness here; checked afterwards).
   usage: cs2 AMIN AMAX CL CH DR PS Y MINRUN */
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
  int64_t AMIN=atoll(argv[1]), AMAX=atoll(argv[2]), CL=atoll(argv[3]), CH=atoll(argv[4]), DR=atoll(argv[5]); int PS=atoi(argv[6]), Y=atoi(argv[7]), MINRUN=atoi(argv[8]);
  LIM=1ULL<<29; bits=malloc(LIM>>4); memset(bits,0xff,LIM>>4); bits[0]&=~2;
  for(uint64_t i=3;i*i<LIM;i+=2) if((bits[i>>4]>>((i>>1)&7))&1) for(uint64_t j=i*i;j<LIM;j+=2*i) bits[j>>4]&=~(1<<((j>>1)&7));
  int pr[32],np=0; for(int p=2;p<=PS;p++){ int q=1; for(int d=2;d*d<=p;d++) if(p%d==0) q=0; if(q) pr[np++]=p; }
  int nsm=0; while(nsm<np && pr[nsm]<=13) nsm++;
  const int64_t M=30030; int64_t crt[6]; for(int i=0;i<nsm;i++){ int64_t p=pr[i], Mi=M/p; int64_t inv=1; while((Mi*inv)%p!=1) inv++; crt[i]=Mi*inv%M; }
  int W=2*Y+1; uint8_t *ok=malloc(W); int64_t *cls=malloc(sizeof(int64_t)*M);
  static uint8_t allow[32][64]; long long tested=0, cand=0;
  for(int64_t a=AMIN;a<=AMAX;a++){ long long runs[100]={0};
   for(int64_t b=0;b<3*a;b++) for(int64_t c=CL*a;c<=CH*a;c++){
    int dead=0;
    for(int i=0;i<np;i++){ int64_t p=pr[i]; memset(allow[i],1,p); int64_t A=md(a,p),B=md(b,p),C=md(c,p);
      for(int64_t y=0;y<p;y++){ int64_t v=((A*y%p*y%p*y)+(B*y%p*y)+C*y)%p; allow[i][(p-v)%p]=0; }
      int any=0; for(int r=0;r<p;r++) any|=allow[i][r]; if(!any){dead=1;break;} }
    if(dead) continue;
    /* residues mod 30030 by CRT over the allowed sets of 2..13 */
    int nc=1; cls[0]=0;
    for(int i=0;i<nsm;i++){ int p=pr[i]; int n0=nc; int k=0; int64_t *tmp=cls+M/2; /* scratch */
      static int64_t buf[30030]; int nb=0;
      for(int j=0;j<n0;j++) for(int r=0;r<p;r++) if(allow[i][r]) buf[nb++]=(cls[j]+crt[i]*r)%M;
      memcpy(cls,buf,sizeof(int64_t)*nb); nc=nb; (void)k; (void)tmp; }
    int64_t D=DR*a; int h=(MINRUN+1)/2;
    for(int j=0;j<nc;j++) for(int64_t d=cls[j]-((D+cls[j])/M+1)*M; d<=D; d+=M){
      if(d<-D) continue; cand++;
      int pass=1; for(int i=nsm;i<np;i++){ if(!allow[i][md(d,pr[i])]){pass=0;break;} }
      if(!pass) continue; tested++;
      int any=0; for(int y=-(Y/h)*h;y+h<=Y;y+=h){ if(isp(((a*y+b)*y+c)*y+d) && isp(((a*(y+h)+b)*(y+h)+c)*(y+h)+d)){any=1;break;} }
      if(!any) continue;
      for(int y=-Y;y<=Y;y++) ok[y+Y]=isp(((a*y+b)*y+c)*y+d);
      int best=0,bs=0,cur=0; for(int q=0;q<W;q++){ cur=ok[q]?cur+1:0; if(cur>best){best=cur;bs=q-cur+1-Y;} }
      if(best<100) runs[best]++;
      if(best>=MINRUN) printf("%lld %lld %lld %lld %d %d\n",(long long)a,(long long)b,(long long)c,(long long)d,bs,best);
    }
   }
   fprintf(stderr,"a=%lld runs>=30..: ",(long long)a); for(int L=30;L<60;L++) if(runs[L]) fprintf(stderr,"%d:%lld ",L,runs[L]); fprintf(stderr,"\n"); fflush(stdout);
  }
  fprintf(stderr,"candidates %lld tested %lld\n",cand,tested); return 0; }
