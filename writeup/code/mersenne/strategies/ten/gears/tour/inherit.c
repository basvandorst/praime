/* Rounds needed to prove every prime up to X from the prime 2 alone, when a prime N may use the primes already
   proved that divide N - 1 or N + 1: F1 | N - 1 and F2 | N + 1 are the parts built from proved primes (full powers),
   and N is proved in the round after lcm(F1, F2)^2 > N holds.  depth2 X [side]   side 0: both, 1: N - 1 only */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
typedef uint64_t u64; typedef unsigned __int128 u128;
static uint8_t *D; static inline int getd(u64 p){if(p==2)return 0;u64 i=p>>1;return (D[i>>1]>>((i&1)*4))&15;}
static inline void setd(u64 p,int d){u64 i=p>>1;D[i>>1]|=(uint8_t)(d<<((i&1)*4));}
#define S 65536
#define MF 12
static u64 rem[S+2]; static u64 fp[S+2][MF], fe[S+2][MF]; static int fc[S+2];
static u64 gcd64(u64 a,u64 b){while(b){u64 t=a%b;a=b;b=t;}return a;}
int main(int argc,char**argv){u64 X=strtoull(argv[1],0,10);int only=argc>2?atoi(argv[2]):0;
  int thr=argc>4?atoi(argv[3]):99;FILE*dump=argc>4?fopen(argv[4],"w"):0;
  D=calloc(X/4+16,1);u64 H=X/2+2;u64 r=(u64)sqrtl((long double)H)+2;
  uint8_t*sm=calloc(r+1,1);u64*pr=malloc((r+1)*8);int np=0;for(u64 i=2;i<=r;i++){if(!sm[i]){pr[np++]=i;for(u64 j=i*i;j<=r;j+=i)sm[j]=1;}}
  u64 cnt[16]={0},first[16]={0},sumd[40]={0},numd[40]={0};first[0]=2;cnt[0]=1;
  /* m runs over [A, A+S]; N = 2m+1 for m in [A, A+S) uses m and m+1 */
  for(u64 A=1;A<H;A+=S){u64 n=S+1;for(u64 i=0;i<n;i++){rem[i]=A+i;fc[i]=0;}
    for(int k=0;k<np;k++){u64 p=pr[k];if(p*p>A+n)break;u64 s=(A+p-1)/p*p;for(u64 m=s;m<A+n;m+=p){u64 i=m-A,pe=1;while(rem[i]%p==0){rem[i]/=p;pe*=p;}if(fc[i]>=MF-1){fprintf(stderr,"too many prime factors: raise MF\n");return 2;}fp[i][fc[i]]=p;fe[i][fc[i]++]=pe;}}
    for(u64 i=0;i<n;i++)if(rem[i]>1){fp[i][fc[i]]=rem[i];fe[i][fc[i]++]=rem[i];}
    for(u64 i=0;i+1<n;i++){u64 m=A+i,N=2*m+1;if(N>X)break;
      /* N prime? use a small test: N is prime iff ... we need primality: check via the nibble of factors is not possible; do trial by Miller-Rabin */
      /* deterministic Miller-Rabin for N < 3.3e9..2^32 with bases 2,3,5,7 (valid below 3,215,031,751) plus 11 */
      {u64 d=N-1;int s2=0;while(!(d&1)){d>>=1;s2++;}int isp=1;static const u64 bs[]={2,3,5,7,11,13};if(N<14){isp=(N==3||N==5||N==7||N==11||N==13);}else{if(N%3==0||N%5==0||N%7==0||N%11==0||N%13==0)isp=0;
        for(int b=0;b<6&&isp;b++){u64 x=1,bb=bs[b],e=d;while(e){if(e&1)x=(u128)x*bb%N;bb=(u128)bb*bb%N;e>>=1;}if(x==1||x==N-1)continue;int ok=0;for(int j=1;j<s2;j++){x=(u128)x*x%N;if(x==N-1){ok=1;break;}}if(!ok)isp=0;}}
        if(!isp)continue;}
      /* items: (depth, prime power, side) from N-1 = 2m and N+1 = 2(m+1) */
      int dd[2*MF+4];u64 pe[2*MF+4];int sd[2*MF+4];int c=0;
      for(int side=0;side<(only?1:2);side++){u64 j=i+side;u64 two=2;int has2=0;for(int t=0;t<fc[j];t++){if(fp[j][t]==2){two=2*fe[j][t];has2=1;}else{dd[c]=getd(fp[j][t]);pe[c]=fe[j][t];sd[c++]=side;}}(void)has2;dd[c]=0;pe[c]=two;sd[c++]=side;}
      for(int a=1;a<c;a++){int d0=dd[a];u64 p0=pe[a];int s0=sd[a];int b=a-1;while(b>=0&&dd[b]>d0){dd[b+1]=dd[b];pe[b+1]=pe[b];sd[b+1]=sd[b];b--;}dd[b+1]=d0;pe[b+1]=p0;sd[b+1]=s0;}
      u64 F1=1,F2=1;int g=-1;for(int a=0;a<c;a++){if(sd[a])F2*=pe[a];else F1*=pe[a];u128 L=(u128)F1*F2/gcd64(F1,F2);if(L*L>N){g=dd[a];break;}}
      if(g<0){fprintf(stderr,"no depth for %llu\n",(unsigned long long)N);return 1;}
      int dep=g+1;if(dep>15){fprintf(stderr,"depth overflow\n");return 1;}setd(N,dep);cnt[dep]++;if(!first[dep]){first[dep]=N;printf("first of depth %d: %llu\n",dep,(unsigned long long)N);fflush(stdout);}
      if(dep>=5&&cnt[dep]<=12)printf("depth %d prime %llu\n",dep,(unsigned long long)N);
      if(dump&&dep>=thr)fprintf(dump,"%d %llu\n",dep,(unsigned long long)N);
      int lg=63-__builtin_clzll(N);sumd[lg]+=dep;numd[lg]++;}
  }
  printf("X %llu side %d\n",(unsigned long long)X,only);for(int d=0;d<16;d++)if(cnt[d])printf("depth %d: %llu primes, first %llu\n",d,(unsigned long long)cnt[d],(unsigned long long)first[d]);
  for(int l=1;l<40;l++)if(numd[l])printf("primes in [2^%d, 2^%d): %llu, mean depth %.4f\n",l,l+1,(unsigned long long)numd[l],(double)sumd[l]/numd[l]);
  if(dump)fclose(dump);printf("done\n");return 0;}
