// Puzzle 1139, monic case: v(N) = |N^2 - D| (or N^2 + D with a 3rd argument "neg"), halved when N is even.
// Usage: ./1139_monic LO HI [neg]   -- scans every admissible D in [LO, HI], D = 2 mod 4.
// gcc -O2 -o 1139_monic 1139_monic.c
// a=1 family: v(N)=|N^2-D|, halved when N even (D = 2 mod 4). A run of >=110 values needs
// (D/p) = -1 for every odd p <= 53. Enumerate such D by CRT and test anchors N=0,55,110,...
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
static uint64_t mulmod(uint64_t a,uint64_t b,uint64_t m){return (unsigned __int128)a*b%m;}
static uint64_t pw(uint64_t a,uint64_t e,uint64_t m){uint64_t r=1;a%=m;while(e){if(e&1)r=mulmod(r,a,m);a=mulmod(a,a,m);e>>=1;}return r;}
static int isprime(uint64_t n){ if(n<2) return 0; static const int sp[]={2,3,5,7,11,13,17,19,23,29,31,37};
  for(int i=0;i<12;i++){ if(n==(uint64_t)sp[i]) return 1; if(n%sp[i]==0) return 0; }
  uint64_t d=n-1;int s=0;while(!(d&1)){d>>=1;s++;}
  for(int i=0;i<12;i++){uint64_t x=pw(sp[i],d,n); if(x==1||x==n-1) continue; int ok=0;
    for(int r=1;r<s;r++){x=mulmod(x,x,n); if(x==n-1){ok=1;break;}} if(!ok) return 0;} return 1; }
static uint64_t D; static int NEG=0;
static uint64_t val(long N){ __int128 v=NEG?(__int128)N*N+(__int128)D:(__int128)N*N-(__int128)D; if(v<0)v=-v; uint64_t u=v; if(!(N&1)&&u!=2) u>>=1; return u; }
int main(int argc,char**argv){
  uint64_t lo=strtoull(argv[1],0,10), hi=strtoull(argv[2],0,10); NEG=argc>3;
  const int P[]={3,5,7,11,13,17,19,23,29,31,37,41,43,47,53}; const int NP=15, NB=7; // base primes 3..19
  static char qr[15][64];
  for(int i=0;i<NP;i++){ for(int r=0;r<P[i];r++) qr[i][r]=0; for(int x=0;x<P[i];x++) qr[i][x*x%P[i]]=1; }
  uint64_t M=4; for(int i=0;i<NB;i++) M*=P[i];
  // residues mod M: D = 2 mod 4 and QNR mod base primes
  uint32_t *cls=malloc(sizeof(uint32_t)*20000); int nc=0;
  for(uint64_t r=2;r<M;r+=4){ int ok=1; for(int i=0;i<NB;i++) if(qr[i][(NEG?P[i]-r%P[i]:r%P[i])%P[i]]){ok=0;break;} if(ok) cls[nc++]=r; }
  fprintf(stderr,"M=%lu classes=%d\n",M,nc);
  uint64_t tested=0; int best=0;
  for(uint64_t base=lo/M*M; base<=hi; base+=M){
    for(int k=0;k<nc;k++){ uint64_t d=base+cls[k]; if(d<lo||d>hi) continue;
      int ok=1; for(int i=NB;i<NP;i++) if(qr[i][(NEG?P[i]-d%P[i]:d%P[i])%P[i]]){ok=0;break;} if(!ok) continue;
      tested++; D=d;
      int bestrun=0,bs=0; long lastend=-1;
      for(long an=0;an<=300;an+=55){ if(an<=lastend||!isprime(val(an))) continue;
        long l=an,h=an; while(l>0&&isprime(val(l-1))) l--; while(h<300&&isprime(val(h+1))) h++;
        lastend=h; if(h-l+1>bestrun){bestrun=h-l+1;bs=l;} }
      int total = (bs==0) ? 2*bestrun-1 : bestrun;
      if(total>best || total>=90){ if(total>best) best=total; printf("D=%lu runN=%d..%d total=%d\n",d,bs,bs+bestrun-1,total); fflush(stdout);} }
  }
  fprintf(stderr,"tested=%lu best=%d\n",tested,best); return 0; }
