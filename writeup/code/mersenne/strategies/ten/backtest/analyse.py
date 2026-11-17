import json,sys
d=json.load(open(sys.argv[1])); CON={'struck','fermat','leftsacks','leftulam','leftklauber','fishbone','parity','gaptowers','bubbles','sunrays','fan','dna','stairs','returnmap','chaosgap','chaoslane','dejavu','goldfirst','famplus','famminus','expclass','per10','per2','mersieve','ord2','primroot','llres','wieferich','cunningham','classrows','survivors','atlas','scales','spectrum','lean','sqwalk','curve','pow2','riesel','primorial','mirror','shield','lucky','glue','omega','dickman','totient','liouville','mertens','pili','chebrace','race4','goldbach','quadform','eisenstein','twinwave'}
ok=[r for r in d if r.get('struct') is not None]
print("layouts scored:",len(ok),"of",len(d))
hon=[r for r in ok if r['layout'] not in CON and r['pairs']>200000]
print("\nplain position maps (place n by a formula, colour by prime): %d"%len(hon))
print("%-18s %8s %8s %7s %6s %6s %6s"%("layout","struct","decoy","ratio","beyond","maxz","expz"))
for r in sorted(hon,key=lambda r:-r['struct']):
    print("%-18s %8.1f %8.1f %7.3f %6.2f %6.1f %6s"%(r['tag'],r['struct'],r['structSieve'],r['struct']/max(r['structSieve'],1e-9),r['beyond'] or 0,r['maxz'],"-" if r['expz'] is None else "%.2f"%r['expz']))
import statistics as st
b=[r['beyond'] for r in hon if r['beyond']]; e=[r['expz'] for r in hon if r['expz'] is not None]
print("\nbeyond: median %.2f, max %.2f (%s); share above 1.5: %d of %d"%(st.median(b),max(b),max(hon,key=lambda r:r['beyond'] or 0)['tag'],sum(x>1.5 for x in b),len(b)))
print("exponent clustering z: mean %.2f sd %.2f min %.2f (%s) max %.2f (%s)"%(st.mean(e),st.pstdev(e),min(e),min((r for r in hon if r['expz'] is not None),key=lambda r:r['expz'])['tag'],max(e),max((r for r in hon if r['expz'] is not None),key=lambda r:r['expz'])['tag']))
rat=[r['struct']/r['structSieve'] for r in hon if r['structSieve']>5]; print("struct / decoy ratio where struct > 5: median %.3f, min %.3f, max %.3f"%(st.median(rat),min(rat),max(rat)))
oth=[r for r in ok if r['layout'] in CON]
print("\nlayouts whose positions or colours use more than n (not comparable):",", ".join(sorted(set(r['layout'] for r in oth))))
low=[r['tag'] for r in ok if r['layout'] not in CON and r['pairs']<=200000]; print("too few neighbouring pairs to score:",", ".join(low))
