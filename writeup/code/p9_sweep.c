/* Puzzle 9: every perfect power among S(k) = p_1 + ... + p_k and T(k) = S(k) - 2 for all primes p_k <= LIMIT.
   Segmented sieve (4 MB segments); squares and cubes tested with exact integer roots, prime exponents 5..53
   with one pointer per exponent (enough for S < 2^64: 2^59 > S). Prints one line per hit.
   gcc -O2 -o p9_sweep writeup/code/p9_sweep.c -lm && nice -n 19 ./p9_sweep 4000000000   (~10 min, 1 core) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
typedef unsigned long long u64;
typedef unsigned __int128 u128;
static u64 iroot(u64 x,int m){u64 r=(u64)powl((long double)x,1.0L/m);while(r>0){u128 p=1;for(int i=0;i<m;i++)p*=r;if(p<=x)break;r--;}for(;;){u128 p=1;for(int i=0;i<m;i++)p*=(r+1);if(p>x)break;r++;}return r;}
static int ispow(u64 x,int m){u64 r=iroot(x,m);u128 p=1;for(int i=0;i<m;i++)p*=r;return p==x;}
static const int HM[]={5,7,11,13,17,19,23,29,31,37,41,43,47,53}; enum{NM=14};
static u64 pS[NM],pT[NM];
static u128 pw(u64 r,int m){u128 p=1;for(int i=0;i<m;i++)p*=r;return p;}
static void check(const char*tag,u64 k,u64 p,u64 s){
  static const int M[]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53};
  for(int i=0;i<16;i++){int m=M[i]; if(s>1 && ispow(s,m)){printf("%s k=%llu p=%llu value=%llu = %llu^%d\n",tag,k,p,s,iroot(s,m),m);fflush(stdout);}}
}
int main(int argc,char**argv){
  u64 LIM=strtoull(argv[1],0,10);
  u64 sq=(u64)sqrtl((long double)LIM)+1;
  char*small=calloc(sq+1,1); for(u64 i=2;i*i<=sq;i++) if(!small[i]) for(u64 j=i*i;j<=sq;j+=i) small[j]=1;
  u64 *sp=malloc(sizeof(u64)*sq); u64 ns=0; for(u64 i=3;i<=sq;i++) if(!small[i]) sp[ns++]=i;
  const u64 SEG=1<<22; char*seg=malloc(SEG);
  u64 k=1,S=2; for(int i=0;i<NM;i++){pS[i]=1;pT[i]=1;} u64 sqn=1; u64 cbn=1;
  for(u64 lo=3;lo<=LIM;lo+=SEG){
    u64 hi=lo+SEG-1; if(hi>LIM) hi=LIM; memset(seg,0,SEG);
    for(u64 i=0;i<ns;i++){u64 q=sp[i]; if(q*q>hi)break; u64 st=(lo+q-1)/q*q; if(st<q*q) st=q*q; for(u64 j=st;j<=hi;j+=q) seg[j-lo]=1;}
    for(u64 x=lo;x<=hi;x++){ if(!(x&1)||seg[x-lo]) continue;
      k++; S+=x; u64 T=S-2;
      while((u128)sqn*sqn<S) sqn++;
      if((u128)sqn*sqn==S) check("S",k,x,S);
      while((u128)cbn*cbn*cbn<T) cbn++;
      u64 r=iroot(T,2); if(r*r==T) check("T",k,x,T); else if((u128)cbn*cbn*cbn==T) check("T",k,x,T);
      { u64 c=iroot(S,3); if((u128)c*c*c==S) check("S",k,x,S); }
      for(int i=0;i<NM;i++){ while(pw(pS[i],HM[i])<S) pS[i]++; if(pw(pS[i],HM[i])==S) check("S",k,x,S); while(pw(pT[i],HM[i])<T) pT[i]++; if(pw(pT[i],HM[i])==T) check("T",k,x,T); }
    }
  }
  fprintf(stderr,"done k=%llu S=%llu\n",k,S); return 0;
}
