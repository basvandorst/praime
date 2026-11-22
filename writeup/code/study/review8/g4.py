import json, math
from scipy import stats
D='/home/bas/praime/writeup/data/study/'
r40=json.load(open(D+'r40.json')); r41=json.load(open(D+'r41.json'))
W={'train':r40['train']['per_gap'],'A':r40['A']['per_gap'],'w89':r41['per_gap']}
for w,per in W.items():
    below=sum(per[g]['meas_u2']<per[g]['pred_u2'] for g in per)
    print(w,'gaps below pred:',below,'of',len(per), 'sign p', stats.binomtest(below,len(per),0.5,alternative='greater').pvalue)
    for g in ('2','4'):
        d=per[g]; n=d['rows']; x=n*d['meas_u2']/d['pred_u2']
        print('  g',g, round(d['meas_u2'],3), round(d['pred_u2'],3), 'p_low chi2_n', round(stats.chi2.cdf(x,n),4))
ps=[]
for w,per in W.items():
    for g in ('2','4'):
        d=per[g]; n=d['rows']; ps.append(stats.chi2.cdf(n*d['meas_u2']/d['pred_u2'],n))
c=-2*sum(map(math.log,ps)); print('Fisher g2+g4 three windows', round(c,2), 'p', stats.chi2.sf(c,12))
# class>=40 arithmetic
Yb=45086078/48; P=1-0.8733758938322979; eps=1-0.9973970854947489
print('Ybar',Yb,'P20',P,'eps2 P Ybar',eps**2*P*Yb,'E',P*Yb,'1/sqrt(E)',1/math.sqrt(P*Yb),'multinomial sd',math.sqrt((1-P)/(P*Yb)))
sd=0.002285924127880695
print('implied mean cell from sd and offset', (sd**2*47/48+eps**2)*P*Yb, 'variance ratio vs multinomial', (sd/math.sqrt((1-P)/(P*Yb)))**2, '+-', math.sqrt(2/47))
Ybt=5096875/48; P20t=465572/5096875; print('train eps2 P Ybar with same eps', eps**2*P20t*Ybt)
