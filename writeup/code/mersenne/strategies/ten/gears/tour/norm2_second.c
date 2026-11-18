/* Second program for the norm-2 wheel, written separately from norm2.c: no Montgomery form (plain 128-bit remainders),
   the power computed from the lowest bit up, the Jacobi symbol by the binary algorithm on signed steps, the square test by Newton.
   norm2_second K I   reads the list on stdin and tests the lines with (line number mod K) == I. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned __int128 u128; typedef uint64_t u64;
static u64 mulm(u64 a,u64 b,u64 n){return (u64)(((u128)a*b)%n);}
static u64 addm(u64 a,u64 b,u64 n){return (u64)(((u128)a+b)%n);}
static u64 isqrt64(u64 n){if(n<2)return n;u64 x=n,y=(x+1)/2;while(y<x){x=y;y=(x+n/x)/2;}return x;}
static int jacobi2(u64 a,u64 n){int t=1;a%=n;while(a!=0){int e=0;while((a&1)==0){a/=2;e++;}if(e&1){u64 r=n%8;if(r==3||r==5)t=-t;}if(a%4==3&&n%4==3)t=-t;u64 tmp=a;a=n%tmp;n=tmp;}return n==1?t:0;}
int main(int argc,char**argv){u64 K=strtoull(argv[1],0,10),I=strtoull(argv[2],0,10);char buf[64];
  u64 line=0,mine=0,tested=0,home=0,sq=0,fac=0,amax=0,last=0,sumc=0;
  while(fgets(buf,64,stdin)){u64 n=strtoull(buf,0,10);line++;last=n;if(line%K!=I)continue;mine++;
    u64 r=isqrt64(n);if(r*r==n){sq++;continue;}
    u64 a=0,c=0;int stop=0;
    for(u64 t=1;!stop&&t<200000;t++){for(int w=0;w<2&&!stop;w++){if(w==0&&t<2)continue;u64 cc=w==0?t*t-2:t*t+2;int j=jacobi2(cc,n);if(j==-1){a=t;c=cc;stop=1;}else if(j==0&&cc%n!=0){stop=2;}}}
    if(stop==2){fac++;continue;}if(stop!=1){fprintf(stderr,"no step for %llu\n",(unsigned long long)n);continue;}
    tested++;if(a>amax)amax=a;sumc+=c;
    u64 cm=c%n,bx=a%n,by=1,rx=1,ry=0,e=n;
    while(e){if(e&1){u64 nx=addm(mulm(rx,bx,n),mulm(cm,mulm(ry,by,n),n),n),ny=addm(mulm(rx,by,n),mulm(ry,bx,n),n);rx=nx;ry=ny;}
      u64 sx=addm(mulm(bx,bx,n),mulm(cm,mulm(by,by,n),n),n),sy=mulm(2%n,mulm(bx,by,n),n);bx=sx;by=sy;e>>=1;}
    if(rx==a%n&&ry==n-1){home++;printf("home %llu\n",(unsigned long long)n);fflush(stdout);}}
  fprintf(stderr,"second program, lines with number mod %llu = %llu: lines read %llu, last %llu, mine %llu, on a wheel %llu, came home %llu, perfect squares %llu, sharing a factor with a small c %llu, largest a %llu, sum of the steps c %llu\ndone\n",(unsigned long long)K,(unsigned long long)I,(unsigned long long)line,(unsigned long long)last,(unsigned long long)mine,(unsigned long long)tested,(unsigned long long)home,(unsigned long long)sq,(unsigned long long)fac,(unsigned long long)amax,(unsigned long long)sumc);return 0;}
