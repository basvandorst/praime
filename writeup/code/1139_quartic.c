// Puzzle 1139, even quartics: v(N) = |N^4 + B N^2 + G|, halved when N is even.
// Usage: ./1139_quartic Bmin Bmax Gmin Gmax K0
// gcc -O2 -o 1139_quartic 1139_quartic.c
// even quartic v(N)=|N^4+B N^2+G|, halved when N even (B even, G = 2 mod 4). Symmetric run N=0..K.
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
static long B,G;
static uint64_t val(long N){ __int128 t=(__int128)N*N; __int128 v=t*t+B*t+G; if(v<0)v=-v; uint64_t u=v; if(!(N&1)&&u!=2)u>>=1; return u; }
int main(int argc,char**argv){
  long Bmin=atol(argv[1]),Bmax=atol(argv[2]),Gmin=atol(argv[3]),Gmax=atol(argv[4]); int K0=atoi(argv[5]);
  const int P[]={3,5,7,11,13,17,19,23}; const int NP=8;
  static char bad[8][23][23]; // bad[i][B mod p][G mod p]: some x has x^4+Bx^2+G = 0 mod p
  for(int i=0;i<NP;i++){int p=P[i]; for(int b=0;b<p;b++) for(int g=0;g<p;g++){ bad[i][b][g]=0; for(int x=0;x<p;x++) if(((long)x*x%p*x%p*x+(long)b*x*x+g)%p==0){bad[i][b][g]=1;break;} }}
  int best=0;
  for(long b=Bmin;b<=Bmax;b+=2){ B=b; int bm[8]; for(int i=0;i<NP;i++) bm[i]=((b%P[i])+P[i])%P[i];
    long g0=Gmin-((Gmin%4+4)%4)+2;
    for(long g=g0;g<=Gmax;g+=4){
      int ok=1; for(int i=0;i<NP;i++){ int gm=((g%P[i])+P[i])%P[i]; if(bad[i][bm[i]][gm]){ok=0;break;} }
      if(!ok) continue; G=g;
      int K=-1; while(K<400 && isprime(val(K+1))) K++;
      if(K>=K0 || 2*K+1>best){ if(2*K+1>best) best=2*K+1; printf("B=%ld G=%ld K=%d total=%d\n",b,g,K,2*K+1); fflush(stdout);} } }
  return 0; }
