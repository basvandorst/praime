// Puzzle 882: A Z..Z B (Z single digit, i copies) prime for i = 0..K; all N_i < 2^64.
// usage: azb Z b Amax  -> prints runs of length >= 10 (i = 0..9 or more); stats of run lengths to stderr
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef unsigned __int128 u128; typedef uint64_t u64;
static u64 mulm(u64 a,u64 b,u64 m){return (u128)a*b%m;}
static u64 powm(u64 a,u64 e,u64 m){u64 r=1;a%=m;while(e){if(e&1)r=mulm(r,a,m);a=mulm(a,a,m);e>>=1;}return r;}
static int mr(u64 n){ if(n<2)return 0; static const u64 sp[]={2,3,5,7,11,13,17,19,23,29,31,37};
  for(int i=0;i<12;i++){if(n%sp[i]==0)return n==sp[i];}
  u64 d=n-1;int s=0;while(!(d&1)){d>>=1;s++;}
  static const u64 B[]={2,325,9375,28178,450775,9780504,1795265022};
  for(int i=0;i<7;i++){u64 a=B[i]%n; if(!a)continue; u64 x=powm(a,d,n); if(x==1||x==n-1)continue; int ok=0;
    for(int r=1;r<s;r++){x=mulm(x,x,n); if(x==n-1){ok=1;break;}} if(!ok)return 0;} return 1;}
#define K 11            /* we want i = 0..10 */
#define NP 600
int main(int c,char**v){ int Z=atoi(v[1]), b=atoi(v[2]); u64 Amax=strtoull(v[3],0,10);
  u64 p10[25]; p10[0]=1; for(int i=1;i<25;i++) p10[i]=p10[i-1]*10;
  int pr[NP]; int np=0; for(int q=3;np<NP;q+=2){int ok=1; for(int d=3;d*d<=q;d+=2) if(q%d==0){ok=0;break;} if(ok&&q!=5) pr[np++]=q;}
  const u64 BL=1<<18; unsigned char*sv=malloc(BL); long hist[K+2]={0}; u64 tested=0;
  u64 Blo=(b==1)?1:p10[b-1], Bhi=p10[b];
  for(u64 B=Blo;B<Bhi;B++){ if(B%2==0||B%5==0) continue; if((u64)(B/p10[b-1])==(u64)Z) continue; // B must not start with Z
    // Amax limited so that N_10 < 2^64 : A*10^(10+b) + ... < 1.8e19
    u64 Alim=Amax; u64 lim=(u64)(18000000000000000000ULL/p10[10+b]); if(lim-1<Alim) Alim=lim-1;
    // forbidden A residues mod q: A == -(Z*R_i*10^b + B) * 10^-(i+b)
    static int forb[NP][K]; static int nf[NP];
    for(int j=0;j<np;j++){ int q=pr[j]; nf[j]=0; u64 Ri=0;
      for(int i=0;i<K;i++){ u64 t=(Z*Ri%q*(p10[b]%q)+B)%q; u64 inv=powm(powm(10,i+b,q),q-2,q); int r=(int)((q-t)%q*inv%q);
        int dup=0; for(int k=0;k<nf[j];k++) if(forb[j][k]==r) dup=1; if(!dup) forb[j][nf[j]++]=r; Ri=(Ri*10+1)%q; } }
    for(u64 lo=1;lo<=Alim;lo+=BL){ u64 n=BL; if(lo+n-1>Alim) n=Alim-lo+1; memset(sv,1,n);
      for(int j=0;j<np;j++){ int q=pr[j]; for(int k=0;k<nf[j];k++){ u64 r=forb[j][k]; u64 st=(r+q-lo%q)%q; for(u64 x=st;x<n;x+=q) sv[x]=0; } }
      for(u64 x=0;x<n;x++) if(sv[x]){ u64 A=lo+x; if(A%10==(u64)Z) continue; tested++;
        u64 N=A*p10[b]+B, Rz=0; int run=0;   // N_i = A*10^(i+b) + Z*R_i*10^b + B
        for(int i=0;i<K;i++){ u64 Ni=A*p10[i+b]+Z*Rz*p10[b]+B; if(!mr(Ni)) break; run++; Rz=Rz*10+1; }
        hist[run]++; if(run>=10){ printf("%d A=%llu Z=%d B=%llu b=%d\n",run,(unsigned long long)A,Z,(unsigned long long)B,b); fflush(stdout);} (void)N; } } }
  fprintf(stderr,"Z=%d b=%d tested=%llu hist:",Z,b,(unsigned long long)tested); for(int i=0;i<=K;i++) fprintf(stderr," %ld",hist[i]); fprintf(stderr,"\n"); }
