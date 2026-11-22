/* Cubics F(y) = a y^3 + b y^2 + c y + d, 0 <= b < 3a (shift-normalised), |c| <= CM, |d| <= DM. Filter: F has no root
   mod p for p in 2..PS (necessary for long runs). Then the longest run of consecutive y in [-Y, Y] with |F| prime.
   usage: cs AMIN AMAX CM DM PS Y MINRUN */
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
int main(int argc, char **argv){
  int64_t AMIN=atoll(argv[1]), AMAX=atoll(argv[2]), CM=atoll(argv[3]), DM=atoll(argv[4]); int PS=atoi(argv[5]), Y=atoi(argv[6]), MINRUN=atoi(argv[7]);
  LIM=1ULL<<29; bits=malloc(LIM>>4); memset(bits,0xff,LIM>>4); bits[0]&=~2;
  for(uint64_t i=3;i*i<LIM;i+=2) if((bits[i>>4]>>((i>>1)&7))&1) for(uint64_t j=i*i;j<LIM;j+=2*i) bits[j>>4]&=~(1<<((j>>1)&7));
  int pr[64],np=0; for(int p=2;p<=PS;p++){ int q=1; for(int d=2;d*d<=p;d++) if(p%d==0) q=0; if(q) pr[np++]=p; }
  const int64_t M=210; int W=2*Y+1; uint8_t *ok=malloc(W); long long tested=0;
  static uint8_t allow[64][128]; int32_t cls[210];
  for(int64_t a=AMIN;a<=AMAX;a++) for(int64_t b=0;b<3*a;b++) for(int64_t c=-CM;c<=CM;c++){
    int dead=0;
    for(int i=0;i<np;i++){ int p=pr[i]; for(int r=0;r<p;r++) allow[i][r]=1;
      for(int y=0;y<p;y++){ int64_t v=((a%p)*y%p*y%p*y + (b%p)*y%p*y + ((c%p)+p)%p*y)%p; allow[i][(p-v)%p]=0; }
      int any=0; for(int r=0;r<p;r++) any|=allow[i][r]; if(!any){dead=1;break;} }
    if(dead) continue;
    int nc=0; for(int r=0;r<M;r++){ int pass=1; for(int i=0;i<np&&pr[i]<=7;i++) if(!allow[i][r%pr[i]]){pass=0;break;} if(pass) cls[nc++]=r; }
    int h=(MINRUN+1)/2;
    for(int64_t k=(-DM)/M-1;k<=DM/M+1;k++) for(int j=0;j<nc;j++){
      int64_t d=k*M+cls[j]; if(d<-DM||d>DM) continue;
      int pass=1; for(int i=0;i<np;i++){ if(pr[i]<=7) continue; int p=pr[i]; if(!allow[i][((d%p)+p)%p]){pass=0;break;} }
      if(!pass) continue; tested++;
      int any=0; for(int y=-(Y/h)*h;y+h<=Y;y+=h){ if(isp(((a*y+b)*y+c)*y+d) && isp(((a*(y+h)+b)*(y+h)+c)*(y+h)+d)){any=1;break;} }
      if(!any) continue;
      for(int y=-Y;y<=Y;y++) ok[y+Y]=isp(((a*y+b)*y+c)*y+d);
      int best=0,bs=0,cur=0; for(int q=0;q<W;q++){ cur=ok[q]?cur+1:0; if(cur>best){best=cur;bs=q-cur+1-Y;} }
      if(best>=MINRUN) printf("%lld %lld %lld %lld %d %d\n",(long long)a,(long long)b,(long long)c,(long long)d,bs,best);
    }
  }
  fprintf(stderr,"tested %lld\n",tested); return 0; }
