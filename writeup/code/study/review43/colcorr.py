import numpy as np
for f in ['an_3e11.npz','an_6e11.npz','an_9e11.npz']:
    d=np.load(f)
    for k in [1,2]:
        P=d['Ph_3_%d'%k]; a=(P[0,1]+P[1,0])/2
        x=d['S1_3_1000000']; c=d['S2col_3_%d_1000000'%k]
        print(f,k,'corr(col, S1) = %.5f'%np.corrcoef(x,c)[0,1], ' mean col/S1 = %.4f  alpha/(1-alpha) = %.4f'%(c.mean()/x.mean(), a/(1-a)))
    # q=10: col vs quadratic form of class counts
