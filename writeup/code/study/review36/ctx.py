import sys,re
pats=sys.argv[2:]; f=sys.argv[1]
txt=open(f).read().split('\n')
for p in pats:
    print('=====',p)
    n=0
    for i,l in enumerate(txt,1):
        for m in re.finditer(re.escape(p),l):
            print(f'{f}:{i}:', l[max(0,m.start()-220):m.end()+150].replace('\n',' ')); n+=1
            if n>=3: break
        if n>=3: break
