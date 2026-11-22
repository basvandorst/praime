import json, math
E = {int(k): v for k, v in json.load(open('explicit.json')).items()}
C = json.load(open('block_counts.json'))
rep = json.load(open('/home/bas/praime/writeup/data/study/r35.json'))
T = 2948.8675772093948
edges = sorted(E)
def sig_repo(x): return (2*math.sqrt(x)/math.log(x))*math.sqrt(math.log(T/(2*math.pi))/(2*math.pi*T))
def sig_rp(x):   # random phases: <cos^2>=1/2, sum_{g>T} 1/g^2 * density = (ln(T/2pi)+1)/(2 pi T)
    return (2*math.sqrt(x)/math.log(x))*math.sqrt(0.5*(math.log(T/(2*math.pi))+1)/(2*math.pi*T))
chi_r = chi_o = 0
print('blk  count  cnt-li  pred-li(mine) repo  diff  z_repo  z(rand-phase sigma, both ends)  n=2,3 terms')
for k in range(10):
    a, b = edges[k], edges[k+1]
    dli = E[b]['li']-E[a]['li']; pred = (E[b]['R']-E[a]['R']) - (E[b]['Z1']-E[a]['Z1'])
    corr = 0.5*(E[b]['Z2']-E[a]['Z2']) + (E[b]['Z3']-E[a]['Z3'])/3
    d = C[k]-pred; s_rep = math.sqrt(2)*sig_repo(b); s_rp = math.hypot(sig_rp(a), sig_rp(b))
    chi_r += (d/s_rep)**2; chi_o += (d/s_rp)**2
    print(k, C[k], round(C[k]-dli), round(pred-dli,1), round(rep['blocks'][k]['prediction_minus_li'],1), round(d,1), round(d/s_rep,2), round(d/s_rp,2), round(corr,2))
a, b = edges[0], edges[-1]
pred = (E[b]['R']-E[a]['R']) - (E[b]['Z1']-E[a]['Z1']); N = sum(C); dli = E[b]['li']-E[a]['li']
s_rp = math.hypot(sig_rp(a), sig_rp(b))
print('whole', N, round(N-dli), round(pred-dli,1), round(N-pred,1), round((N-pred)/(math.sqrt(2)*sig_repo(b)),2), round((N-pred)/s_rp,2))
print('chi2 of 10 blocks: repo sigma %.2f, random-phase sigma %.2f (10 dof)' % (chi_r, chi_o))
print('sigma repo(b=2e9)=%.1f, rand-phase=%.1f, ratio %.3f' % (math.sqrt(2)*sig_repo(2e9), math.sqrt(2)*sig_rp(2e9), sig_repo(2e9)/sig_rp(2e9)))
# R - li per block, and rms swings
import statistics
sw = [C[k]-(E[edges[k+1]]['li']-E[edges[k]]['li']) for k in range(10)]
print('rms(count-li) %.0f, sqrt(mean count) %.0f' % (math.sqrt(sum(x*x for x in sw)/10), math.sqrt(sum(C)/10)))
for k in range(10):
    a,b=edges[k],edges[k+1]; print('R-li per block', k, round((E[b]['R']-E[a]['R'])-(E[b]['li']-E[a]['li']),1))
# Montgomery-Soundararajan variance for H=1e8 at x=1.5e9
H, x = 1e8, 1.5e9; L = math.log(x); B = 1 - 0.5772156649 - math.log(2*math.pi)
print('MS var/mean H=1e8: %.3f -> sd %.0f' % ((math.log(x/H)+B)/L, math.sqrt((math.log(x/H)+B)/L*H/L)))
