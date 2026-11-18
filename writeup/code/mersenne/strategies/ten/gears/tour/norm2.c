/* The two-cell wheel from a start of norm 2: (a + sqrt c)^n = a - sqrt c (mod n),
   c = a*a - 2 or a*a + 2, the first with Jacobi symbol (c/n) = -1 (a = 1, 2, 3, ...; a*a - 2 first, from a = 2).
   norm2 list [-s K]        reads numbers on stdin; with -s K also prints every K-th line as "n a c x y home"
   norm2 range LO HI        every odd n in [LO, HI): counts the numbers that come home
   norm2 rangeold LO HI     the same for machine 11's wheel (1 + sqrt c, c the smallest odd prime with (c/n) = -1) */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
typedef unsigned __int128 u128; typedef uint64_t u64;
static u64 N,NI,ONE,R2;
static inline u64 mont_inv(u64 n){u64 x=n;for(int i=0;i<6;i++)x*=2-n*x;return x;}
static inline u64 redc(u128 t){u64 m=(u64)t*NI;u128 mn=(u128)m*N;u64 r=(u64)((t-mn)>>64);return (t<mn)?r+N:r;}
static inline u64 mm(u64 a,u64 b){return redc((u128)a*b);}
static inline u64 ad(u64 a,u64 b){u64 s=a+b;return (s>=N||s<a)?s-N:s;}
static void setn(u64 n){N=n;NI=mont_inv(n);ONE=(u64)(((u128)1<<64)%n);R2=(u64)(((u128)ONE*ONE)%n);}
static int jac(u64 a,u64 n){a%=n;int r=1;while(a){while(!(a&1)){a>>=1;u64 k=n&7;if(k==3||k==5)r=-r;}u64 t=a;a=n;n=t;if((a&3)==3&&(n&3)==3)r=-r;a%=n;}return n==1?r:0;}
static u64 gcd64(u64 a,u64 b){while(b){u64 t=a%b;a=b;b=t;}return a;}
static int issq(u64 n){u64 r=(u64)sqrtl((long double)n);while(r*r>n)r--;while((r+1)*(r+1)<=n&&(r+1)*(r+1)>r*r)r++;return r*r==n;}
static const int OP[]={3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,0};
/* returns 1 home, 0 not home, -1 perfect square, -2 shares a factor with a small c; fills *pa,*pc,*px,*py */
static int wheel(u64 n,int old,u64*pa,u64*pc,u64*px,u64*py){
  if(issq(n))return -1;
  u64 a=0,c=0;
  if(old){for(int i=0;OP[i];i++){int j=jac(OP[i],n);if(j==0){if(n==(u64)OP[i])continue;return -2;}if(j==-1){c=OP[i];break;}}if(!c)return -1;a=1;}
  else for(u64 t=1;;t++){u64 cs[2];int nc=0;if(t>=2)cs[nc++]=t*t-2;cs[nc++]=t*t+2;int f=0;
    for(int i=0;i<nc;i++){int j=jac(cs[i],n);if(j==-1){a=t;c=cs[i];f=1;break;}if(j==0&&cs[i]%n){*pa=t;*pc=cs[i];return -2;}}
    if(f)break;if(t>100000)return -1;}
  setn(n);u64 aM=mm(a%n,R2),cM=mm(c%n,R2),x=aM,y=ONE;
  int top=63-__builtin_clzll(n);
  for(int i=top-1;i>=0;i--){u64 xx=mm(x,x),yy=mm(y,y),xy=mm(x,y);x=ad(xx,mm(cM,yy));y=ad(xy,xy);
    if((n>>i)&1){u64 nx=ad(mm(aM,x),mm(cM,y)),ny=ad(x,mm(aM,y));x=nx;y=ny;}}
  *pa=a;*pc=c;*px=redc(x);*py=redc(y);
  return x==aM&&y==N-ONE;}
static void pw(u64 n,u64 a,u64 c,u64*px,u64*py){u64 aM=mm(a%n,R2),cM=mm(c%n,R2),x=aM,y=ONE;int top=63-__builtin_clzll(n);
  for(int i=top-1;i>=0;i--){u64 xx=mm(x,x),yy=mm(y,y),xy=mm(x,y);x=ad(xx,mm(cM,yy));y=ad(xy,xy);if((n>>i)&1){u64 nx=ad(mm(aM,x),mm(cM,y)),ny=ad(x,mm(aM,y));x=nx;y=ny;}}*px=x;*py=y;}
static int sprp(u64 n,u64 b){setn(n);u64 d=n-1;int s=0;while(!(d&1)){d>>=1;s++;}u64 x=ONE,B=mm(b%n,R2);if(B==0)return 1;for(int i=63-__builtin_clzll(d);i>=0;i--){x=mm(x,x);if((d>>i)&1)x=mm(x,B);}
  if(x==ONE||x==N-ONE)return 1;for(int i=1;i<s;i++){x=mm(x,x);if(x==N-ONE)return 1;}return 0;}
static const u64 B7[]={2,325,9375,28178,450775,9780504,1795265022ULL};
int main(int argc,char**argv){
  if(argc>=5&&!strcmp(argv[1],"census")){u64 A=strtoull(argv[2],0,10),K=strtoull(argv[3],0,10),I=strtoull(argv[4],0,10),SK=argc>=6?strtoull(argv[5],0,10):0;
    u64 sa[64],sc[64];int ns=0;for(u64 t=1;t<=A;t++){if(t>=2){sa[ns]=t;sc[ns++]=t*t-2;}sa[ns]=t;sc[ns++]=t*t+2;}
    static u64 mt[64],mh[64],ft[64],fp[64],fm[64],sh[64];char buf[64];u64 line=0,mine=0,last=0;
    while(fgets(buf,64,stdin)){u64 n=strtoull(buf,0,10);line++;last=n;if(line%K!=I||n<3||n%2==0)continue;mine++;setn(n);int smp=SK&&line%SK==I;
      for(int i=0;i<ns;i++){int j=jac(sc[i],n);if(j==0){sh[i]++;continue;}u64 x,y,aM=mm(sa[i]%n,R2);pw(n,sa[i],sc[i],&x,&y);int minus=(x==aM&&y==N-ONE),plus=(x==aM&&y==ONE);
        if(j==-1){mt[i]++;if(minus){mh[i]++;printf("misfit home %llu a %llu c %llu\n",(unsigned long long)n,(unsigned long long)sa[i],(unsigned long long)sc[i]);}}else{ft[i]++;fp[i]+=plus;fm[i]+=minus;}
        if(smp)printf("sample %llu %llu %llu %d %llu %llu\n",(unsigned long long)n,(unsigned long long)sa[i],(unsigned long long)sc[i],j,(unsigned long long)redc(x),(unsigned long long)redc(y));}}
    fprintf(stderr,"census of starts up to a = %llu, lines with number mod %llu = %llu: lines read %llu, last %llu, mine %llu\n",(unsigned long long)A,(unsigned long long)K,(unsigned long long)I,(unsigned long long)line,(unsigned long long)last,(unsigned long long)mine);
    for(int i=0;i<ns;i++)fprintf(stderr,"start %llu + sqrt %llu: does not fit %llu, of these home %llu | fits %llu, of these (a + sqrt c)^n = a + sqrt c %llu, = a - sqrt c %llu | shares a factor %llu\n",(unsigned long long)sa[i],(unsigned long long)sc[i],(unsigned long long)mt[i],(unsigned long long)mh[i],(unsigned long long)ft[i],(unsigned long long)fp[i],(unsigned long long)fm[i],(unsigned long long)sh[i]);
    fprintf(stderr,"done\n");return 0;}
  if(argc>=2&&!strcmp(argv[1],"bench")){static u64 v[2000000];int m=0;char buf[64];while(m<2000000&&fgets(buf,64,stdin))v[m++]=strtoull(buf,0,10);
    for(int mode=0;mode<4;mode++){struct timespec t0,t1;clock_gettime(CLOCK_MONOTONIC,&t0);u64 yes=0;
      for(int i=0;i<m;i++){u64 n=v[i],a,c,x,y;if(mode==0)yes+=sprp(n,2);else if(mode==1){int ok=1;for(int j=0;j<7&&ok;j++)ok=sprp(n,B7[j]);yes+=ok;}else yes+=wheel(n,mode==3,&a,&c,&x,&y)==1;}
      clock_gettime(CLOCK_MONOTONIC,&t1);double ns=((t1.tv_sec-t0.tv_sec)*1e9+(t1.tv_nsec-t0.tv_nsec))/m;
      fprintf(stderr,"%s: %d numbers, said prime %llu, %.0f ns per number\n",mode==0?"doubling wheel (strong test, base 2)":mode==1?"seven strong tests (bases 2, 325, 9375, 28178, 450775, 9780504, 1795265022)":mode==2?"norm-2 wheel":"machine 11's wheel",m,(unsigned long long)yes,ns);}
    return 0;}
  if(argc>=2&&!strcmp(argv[1],"list")){u64 K=0;if(argc>=4&&!strcmp(argv[2],"-s"))K=strtoull(argv[3],0,10);
    char buf[64];u64 lines=0,tested=0,home=0,sq=0,fac=0,last=0,amax=0,amaxn=0,hist[64]={0};
    while(fgets(buf,64,stdin)){u64 n=strtoull(buf,0,10);lines++;if(n<3||n%2==0)continue;last=n;u64 a=0,c=0,x=0,y=0;int r=wheel(n,0,&a,&c,&x,&y);
      if(r==-1){sq++;printf("square %llu\n",(unsigned long long)n);}else if(r==-2)fac++;else{tested++;if(a>amax){amax=a;amaxn=n;}hist[a<63?a:63]++;if(r==1){home++;printf("home %llu a %llu c %llu\n",(unsigned long long)n,(unsigned long long)a,(unsigned long long)c);}}
      if(K&&lines%K==0)printf("sample %llu %llu %llu %llu %llu %d\n",(unsigned long long)n,(unsigned long long)a,(unsigned long long)c,(unsigned long long)x,(unsigned long long)y,r);
      if((lines&0xFFFFFF)==0){fprintf(stderr,"  .. %llu lines, at %llu\n",(unsigned long long)lines,(unsigned long long)n);fflush(stderr);}}
    fprintf(stderr,"lines %llu, last %llu, on a wheel %llu, came home %llu, perfect squares %llu, sharing a factor with a small c %llu, largest a %llu (n = %llu)\n",(unsigned long long)lines,(unsigned long long)last,(unsigned long long)tested,(unsigned long long)home,(unsigned long long)sq,(unsigned long long)fac,(unsigned long long)amax,(unsigned long long)amaxn);
    fprintf(stderr,"a used:");for(int i=1;i<64;i++)if(hist[i])fprintf(stderr," %d:%llu",i,(unsigned long long)hist[i]);fprintf(stderr,"\ndone\n");return 0;}
  if(argc>=4&&(!strcmp(argv[1],"range")||!strcmp(argv[1],"rangeold"))){int old=!strcmp(argv[1],"rangeold");u64 lo=strtoull(argv[2],0,10),hi=strtoull(argv[3],0,10),home=0,sq=0,fac=0,amax=0,amaxn=0;if(lo<3)lo=3;if(lo%2==0)lo++;
    for(u64 n=lo;n<hi;n+=2){u64 a=0,c=0,x=0,y=0;int r=wheel(n,old,&a,&c,&x,&y);if(r==1){home++;if(a>amax){amax=a;amaxn=n;}}else if(r==-1)sq++;else if(r==-2)fac++;}
    fprintf(stderr,"%s odd n in [%llu, %llu): came home %llu, perfect squares %llu, sharing a factor with a small c %llu, largest a among those home %llu (n = %llu)\ndone\n",old?"old wheel,":"norm-2 wheel,",(unsigned long long)lo,(unsigned long long)hi,(unsigned long long)home,(unsigned long long)sq,(unsigned long long)fac,(unsigned long long)amax,(unsigned long long)amaxn);return 0;}
  return 1;}
