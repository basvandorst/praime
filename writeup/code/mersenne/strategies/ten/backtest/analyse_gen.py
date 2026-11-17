import json,sys,statistics as st,math
d=[r for r in json.load(open(sys.argv[1])) if r['in'] and r['out'] and r['in']['beyond'] and r['out']['beyond']]
print("maps scored:",len(d))
def nm(sp): return {"rows":lambda s:"rows w=%d"%s['w'],"spiral":lambda s:"spiral n^%.3f x %.4f"%(s['a'],s['c']),"phyl":lambda s:"turn %.6f"%s['g'],"quad":lambda s:"quad k=%d"%s['k'],"pow":lambda s:"n = k r^e + x, e=%s k=%d"%(s['e'],s['k']),"tri":lambda s:"triangular k=%d"%s['k'],"cube":lambda s:"cube s=%d shear %d,%d"%(s['s'],s['a'],s['b'])}[sp['t']](sp)
bi=[r['in']['beyond'] for r in d]; bo=[r['out']['beyond'] for r in d]
print("beyond, in-sample (n <= 5e5): mean %.2f, sd %.2f, max %.2f; out-of-sample (5e5 < n <= 1e6): mean %.2f, sd %.2f, max %.2f"%(st.mean(bi),st.pstdev(bi),max(bi),st.mean(bo),st.pstdev(bo),max(bo)))
def corr(a,b):
    ma,mb=st.mean(a),st.mean(b); return sum((x-ma)*(y-mb) for x,y in zip(a,b))/math.sqrt(sum((x-ma)**2 for x in a)*sum((y-mb)**2 for y in b))
print("correlation in/out: beyond %.2f, structure %.3f"%(corr(bi,bo),corr([r['in']['struct'] for r in d],[r['out']['struct'] for r in d])))
top=sorted(d,key=lambda r:-r['in']['beyond'])[:15]
print("the 15 best by 'beyond' in-sample, and the same maps out-of-sample:")
for r in top: print("   %-28s in %.2f (max z %+.1f)   out %.2f (max z %+.1f)"%(nm(r['sp']),r['in']['beyond'],r['in']['maxz'],r['out']['beyond'],r['out']['maxz']))
print("   mean of the 15: in %.2f, out %.2f; mean of all maps out-of-sample: %.2f"%(st.mean(r['in']['beyond'] for r in top),st.mean(r['out']['beyond'] for r in top),st.mean(bo)))
rat=[r['in']['struct']/r['in']['decoy'] for r in d if r['in']['decoy']>5]; print("structure primes / decoys (in-sample, structure > 5): median %.3f, 5%%..95%% %.3f..%.3f, n=%d"%(st.median(rat),sorted(rat)[len(rat)//20],sorted(rat)[-len(rat)//20],len(rat)))
print("the 12 with most structure (millibits, in / out):")
for r in sorted(d,key=lambda r:-r['in']['struct'])[:12]: print("   %-28s %.1f / %.1f   decoys %.1f / %.1f"%(nm(r['sp']),r['in']['struct'],r['out']['struct'],r['in']['decoy'],r['out']['decoy']))
for t in ("rows","spiral","phyl","quad","cube","pow","tri"):
    g=[r for r in d if r["sp"]["t"]==t]
    if not g: continue
    print("  family %-6s n=%3d  structure median %.1f max %.1f; beyond median %.2f"%(t,len(g),st.median(r['in']['struct'] for r in g),max(r['in']['struct'] for r in g),st.median(r['in']['beyond'] for r in g)))
