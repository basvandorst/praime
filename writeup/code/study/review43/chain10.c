#include <stdio.h>
#include <stdlib.h>
/* simulate a 4-state first-order chain; write states as bytes to stdout */
int main(int argc,char**argv){ long n=atol(argv[1]); unsigned long long x=88172645463325252ULL^atol(argv[2]);
 double P[4][4]; for(int i=0;i<4;i++)for(int j=0;j<4;j++)P[i][j]=atof(argv[3+4*i+j]);
 int s=0; unsigned char*buf=malloc(n);
 for(long t=0;t<n;t++){ x^=x<<13;x^=x>>7;x^=x<<17; double u=(x>>11)*(1.0/9007199254740992.0);
  double c=0; int j=0; for(j=0;j<3;j++){c+=P[s][j]; if(u<c)break;} s=j; buf[t]=(unsigned char)s;}
 fwrite(buf,1,n,stdout); return 0;}
