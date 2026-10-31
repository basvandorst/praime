#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define L 200000000
static unsigned char comp[L+1];
int main(int c,char**v){ int MM=atoi(v[1]); comp[0]=comp[1]=1;
 for(long i=2;i*i<=L;i++) if(!comp[i]) for(long j=i*i;j<=L;j+=i) comp[j]=1;
 static int tw[20000000]; int nt=0; for(int p=3;p+2<=L;p+=2) if(!comp[p]&&!comp[p+2]) { if(nt<20000000) tw[nt++]=p; }
 long worst=0; int wm=0, fail=0;
 for(int m=1;m<=MM;m++){ long h=6L*m; int found=0; for(int i=0;i<nt;i++){ long P=tw[i]; if(P+h+2>L) break; if(!comp[P+h]&&!comp[P+h+2]){ if(P>worst){worst=P;wm=m;} found=1; break;} } if(!found){ fail++; printf("no P for m=%d\n",m);} }
 printf("m<=%d fails=%d worst m=%d P=%ld\n",MM,fail,wm,worst);
 for(int m=1;m<=5;m++){ long h=6L*m; for(int i=0;i<nt;i++){ long P=tw[i]; if(!comp[P+h]&&!comp[P+h+2]){ printf("m=%d P=%ld n=%ld\n",m,P,P*(P+h+2)); break;} } }
}
