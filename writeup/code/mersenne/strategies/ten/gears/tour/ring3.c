#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef unsigned __int128 u128; typedef uint64_t u64;
#define MR 128
static u64 N,NI,ONE,R2;
static inline u64 mont_inv(u64 n){u64 x=n;for(int i=0;i<6;i++)x*=2-n*x;return x;}
static inline u64 redc(u128 t){u64 m=(u64)t*NI;u128 mn=(u128)m*N;u64 r=(u64)((t-mn)>>64);return (t<mn)?r+N:r;}
static inline u64 ad(u64 a,u64 b){u64 s=a+b;return (s>=N||s<a)?s-N:s;}
static inline u64 sb(u64 a,u64 b){return a>=b?a-b:a+(N-b);}
static void setn(u64 n){N=n;NI=mont_inv(n);ONE=(u64)(((u128)1<<64)%n);R2=(u64)(((u128)ONE*ONE)%n);}
static inline u64 fin(u128 lo,u64 hi){u128 s=lo+(u128)hi*R2;if(s<lo)s+=R2;u64 sh=(u64)(s>>64);if(sh>=N)sh%=N;return redc(((u128)sh<<64)|(u64)s);}
typedef struct{int np[2],ns[2];unsigned char p[2][MR][2],s[2][MR];}Out;
typedef struct{int len,nout,kind,r;Out o[MR];}Ring;
static Ring*mk(int kind,int r){Ring*R=calloc(1,sizeof(Ring));R->kind=kind;R->r=r;int len=r;R->len=len;R->nout=kind==1?(r+1)/2:len;
  for(int k=0;k<R->nout;k++){Out*o=&R->o[k];for(int i=0;i<len;i++){int j,sg=0;
      if(kind==2){j=k-i;if(j<0){j+=len;sg=1;}}else{j=((k-i)%len+len)%len;}
      if(i<j){o->p[sg][o->np[sg]][0]=i;o->p[sg][o->np[sg]][1]=j;o->np[sg]++;}else if(i==j){o->s[sg][o->ns[sg]++]=i;}}}
  return R;}
static inline u64 acc(Out*o,int sg,const u64*a){u128 lo=0;u64 hi=0;for(int t=0;t<o->np[sg];t++){u128 p=(u128)a[o->p[sg][t][0]]*a[o->p[sg][t][1]];lo+=p;hi+=lo<p;}
  hi=(hi<<1)|(u64)(lo>>127);lo<<=1;for(int t=0;t<o->ns[sg];t++){u64 v=a[o->s[sg][t]];u128 p=(u128)v*v;lo+=p;hi+=lo<p;}return fin(lo,hi);}
static void sq(Ring*R,const u64*a,u64*c){int len=R->len;for(int k=0;k<R->nout;k++){Out*o=&R->o[k];u64 v=acc(o,0,a);if(R->kind==2)v=sb(v,acc(o,1,a));c[k]=v;}
  if(R->kind==1)for(int k=R->nout;k<len;k++)c[k]=c[len-k];}
static int DUMP=0;
static int test(Ring*R){int len=R->len,kind=R->kind;u64 cur[MR],t[MR];memset(cur,0,sizeof(u64)*len);
  if(kind==1){cur[0]=ad(ONE,ONE);cur[1]=ONE;cur[len-1]=ONE;}else{cur[0]=ONE;cur[1]=ONE;}
  u64 e=N;int top=63-__builtin_clzll(e);
  for(int i=top-1;i>=0;i--){sq(R,cur,t);
    if((e>>i)&1){
      if(kind==0)for(int k=0;k<len;k++)cur[k]=ad(t[k],t[k==0?len-1:k-1]);
      else if(kind==1)for(int k=0;k<len;k++)cur[k]=ad(ad(t[k],t[k]),ad(t[k==0?len-1:k-1],t[k==len-1?0:k+1]));
      else{cur[0]=sb(t[0],t[len-1]);for(int k=1;k<len;k++)cur[k]=ad(t[k],t[k-1]);}}
    else memcpy(cur,t,sizeof(u64)*len);}
  if(DUMP){printf("%llu %d %d",(unsigned long long)N,R->r,kind);for(int k=0;k<len;k++)printf(" %llu",(unsigned long long)redc(cur[k]));printf("\n");}
  u64 want[MR];memset(want,0,sizeof(u64)*len);
  if(kind==0){want[0]=ONE;int j=(int)(N%len);want[j]=ad(want[j],ONE);}
  else if(kind==1){want[0]=ad(ONE,ONE);int j=(int)(N%len);want[j]=ad(want[j],ONE);want[(len-j)%len]=ad(want[(len-j)%len],ONE);}
  else{want[0]=ONE;int j=(int)(N%(2*len));if(j<len)want[j]=ad(want[j],ONE);else want[j-len]=sb(want[j-len],ONE);}
  for(int k=0;k<len;k++)if(cur[k]!=want[k])return 0;return 1;}
static u64 gcd64(u64 a,u64 b){while(b){u64 t=a%b;a=b;b=t;}return a;}
int main(int argc,char**argv){int ai=1,K=1,I=0,mode=0;
  while(ai<argc&&argv[ai][0]=='-'){if(!strcmp(argv[ai],"-s")){K=atoi(argv[ai+1]);I=atoi(argv[ai+2]);ai+=3;}else if(!strcmp(argv[ai],"-fit")){mode=1;ai++;}else if(!strcmp(argv[ai],"-plain")){mode=2;ai++;}else if(!strcmp(argv[ai],"-dump")){DUMP=1;ai++;}else return 1;}
  int nr=argc-ai,rs[64];Ring*F[64],*P[64];u64 T[64]={0},S[64]={0},L[64]={0};
  for(int i=0;i<nr;i++){int r=atoi(argv[ai+i]);rs[i]=r;if(r<3||r>=MR)return 1;F[i]=mk(0,r);if(r&1)P[i]=mk(1,r);else if((r&(r-1))==0)P[i]=mk(2,r/2);else return 1;}
  char buf[64];u64 line=0,tested=0,last=0;
  while(fgets(buf,64,stdin)){u64 n=strtoull(buf,0,10);line++;if((int)(line%K)!=I%K)continue;if(n%2==0||n<3)continue;last=n;tested++;int set=0;
    for(int i=0;i<nr;i++){int r=rs[i];if(gcd64(n,r)!=1)continue;
      if(mode==0){u64 q=n%r;if((q*q)%r==1)continue;int cov=0;for(int d=3;d<r;d++)if(r%d==0){u64 q2=n%d;if((q2*q2)%d!=1){cov=1;break;}}if(cov)continue;}
      if(!set){setn(n);set=1;}T[i]++;
      if(mode!=2){if(!test(P[i]))continue;S[i]++;}
      if(test(F[i])){L[i]++;printf("ring %d liar %llu\n",r,(unsigned long long)n);fflush(stdout);}}
    if((tested&0x3FFFFF)==0){fprintf(stderr,"  .. %llu numbers, at %llu\n",(unsigned long long)tested,(unsigned long long)n);}}
  fprintf(stderr,"numbers read %llu, last %llu\n",(unsigned long long)tested,(unsigned long long)last);
  for(int i=0;i<nr;i++)fprintf(stderr,"  ring %d: tested %llu, passed the first step %llu, came home %llu\n",rs[i],(unsigned long long)T[i],(unsigned long long)S[i],(unsigned long long)L[i]);
  return 0;}
