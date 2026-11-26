import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qcounts as Q
HERE = os.path.dirname(os.path.abspath(__file__))
for r2 in (int(a) for a in sys.argv[1].split(',')):
    Q.counts(os.path.join(HERE, 'data', 'disc_r%d' % r2), 5000, 20000, r2, float(sys.argv[2]), log=lambda s: print(s, flush=True))
