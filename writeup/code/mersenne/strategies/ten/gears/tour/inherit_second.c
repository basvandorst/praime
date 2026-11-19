/* A second program for the rounds of machine 23, written separately from inherit.c (by a reviewing agent, from the
   description of the rule only): primality and factors come from one segmented sieve over all integers; no Miller-Rabin.
   usage: inherit_second X1 X2   rounds are stored for primes <= X1 (needs X1 >= (X2+1)/2), computed for primes <= X2;
   prints counts per round, class means, and every prime of 5 or more rounds.
   L can exceed 2^64 for N above 6.07e9, so L > N is tested before L * L (the square would not fit in 128 bits).   gcc -O2 -o inherit_second inherit_second.c -lm */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
typedef unsigned long long u64; typedef unsigned __int128 u128;
#define W 32768
#define MAXF 10
static unsigned char *R;                 /* nibble per odd number */
static u64 X1;
static int getr(u64 p){ if(p==2) return 0; if(p>X1){fprintf(stderr,"lookup beyond X1: %llu\n",p); exit(2);} u64 i=p>>1; int v=(R[i>>1]>>(4*(i&1)))&15; if(v==0){fprintf(stderr,"round of %llu not set\n",p); exit(3);} return v; }
static void setr(u64 p,int v){ if(p>X1) return; u64 i=p>>1; R[i>>1]|=(unsigned char)(v<<(4*(i&1))); }
static u64 g64(u64 a,u64 b){ while(b){u64 t=a%b;a=b;b=t;} return a; }
static u64 rem[W+2], fpw[W+2][MAXF]; static int fr[W+2][MAXF], nf[W+2]; static u64 fprime[W+2][MAXF];
int main(int argc,char**argv){
  X1=strtoull(argv[1],0,10); u64 X2=strtoull(argv[2],0,10);
  if(2*X1 < X2+1){fprintf(stderr,"X1 too small\n");return 1;}
  R=calloc(X1/4+8,1);
  u64 lim=(u64)sqrt((double)(X2+2))+2; unsigned char *c=calloc(lim+1,1); u64 *pr=malloc(8*(lim+1)); int np=0;
  for(u64 i=2;i<=lim;i++) if(!c[i]){ pr[np++]=i; for(u64 j=i*i;j<=lim;j+=i) c[j]=1; }
  u64 cnt[16]={0}, first[16]={0}; cnt[0]=1; first[0]=2; double sum[64]={0}; u64 num[64]={0};
  /* window holds the integers A .. A+W+1; primes N = A+1 .. A+W are decided in it */
  for(u64 A=2; A<=X2; A+=W){
    for(int i=0;i<W+2;i++){ rem[i]=A+i; nf[i]=0; }
    for(int k=0;k<np;k++){ u64 p=pr[k]; if(p*p>A+W+1) break; u64 s=((A+p-1)/p)*p;
      for(u64 m=s;m<=A+W+1;m+=p){ int i=(int)(m-A); u64 pe=1; while(rem[i]%p==0){rem[i]/=p;pe*=p;} if(nf[i]>=MAXF){fprintf(stderr,"more than %d prime factors: raise MAXF\n",MAXF);exit(7);} fprime[i][nf[i]]=p; fpw[i][nf[i]]=pe; nf[i]++; } }
    for(int i=1;i<=W;i++){ u64 N=A+i; if(N>X2) break; if(!(N&1)) continue;
      if(!((nf[i]==0 && rem[i]==N) || (nf[i]==1 && fprime[i][0]==N))) continue;   /* N prime: no smaller sieve prime divides it */
      /* gather the prime powers of N-1 (side 0) and N+1 (side 1) with the round of their prime */
      u64 pw[2*MAXF+2]; int rd[2*MAXF+2], sd[2*MAXF+2]; int n=0, maxr=0;
      for(int side=0; side<2; side++){ int j=side? i+1 : i-1;
        for(int t=0;t<nf[j];t++){ pw[n]=fpw[j][t]; rd[n]=getr(fprime[j][t]); sd[n]=side; if(rd[n]>maxr)maxr=rd[n]; n++; }
        if(rem[j]>1){ pw[n]=rem[j]; rd[n]=getr(rem[j]); sd[n]=side; if(rd[n]>maxr)maxr=rd[n]; n++; } }
      int res=-1;
      for(int g=0; g<=maxr; g++){ u64 F1=1,F2=1; for(int t=0;t<n;t++) if(rd[t]<=g){ if(sd[t]) F2*=pw[t]; else F1*=pw[t]; }
        u128 L=(u128)F1/g64(F1,F2)*F2; if(L>(u128)N||L*L>(u128)N){ res=g+1; break; } }
      if(res<0||res>15){fprintf(stderr,"no round for %llu\n",N);return 4;}
      setr(N,res); cnt[res]++; if(!first[res]) first[res]=N;
      if(res>=5) printf("deep %d %llu\n",res,N);
      int lg=63-__builtin_clzll(N); sum[lg]+=res; num[lg]++; }
  }
  printf("X2 %llu\n",X2); u64 tot=0;
  for(int d=0;d<16;d++) if(cnt[d]){ printf("round %d: %llu primes, first %llu\n",d,cnt[d],first[d]); tot+=cnt[d]; }
  printf("total primes %llu\n",tot);
  for(int l=1;l<64;l++) if(num[l]) printf("class 2^%d: %llu primes, mean %.4f\n",l,num[l],sum[l]/num[l]);
  return 0; }
