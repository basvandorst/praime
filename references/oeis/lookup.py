#!/usr/bin/env python3
"""Offline OEIS lookup against the local bulk files in this directory.

  python3 lookup.py 2,3,5,7,11,13      sequences containing these terms consecutively
  python3 lookup.py 2 3 5 7 11 13      same (space-separated)
  python3 lookup.py A000040            name of one A-number
  python3 lookup.py A000040 --full     the whole entry (terms, formulas, comments, programs) from the
                                       oeisdata clone; a missing entry is added to the sparse checkout
  python3 lookup.py -w Sophie Germain  sequences whose name contains all these words (any case)
  python3 lookup.py --update           download fresh stripped.gz and names.gz (OEIS rebuilds them daily)
  options: -n N  show at most N matches (default 30)

Needs stripped.gz and names.gz from https://oeis.org/ (not in git; refreshed daily):
  curl -O https://oeis.org/stripped.gz && curl -O https://oeis.org/names.gz
The term match is on the comma-separated string from stripped.gz, so
"3,5,7" matches ",3,5,7," and never "13,5,7". stripped.gz holds only the
first ~260 characters of terms that OEIS shows in %S/%T/%U.
"""
import gzip, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))

def need(fn):
    path = os.path.join(HERE, fn)
    if not os.path.exists(path):
        sys.exit(f"{path} is missing. It is not in git (large, refreshed daily).\n"
                 f"Download it first:  cd {HERE} && curl -O https://oeis.org/{fn}")
    return path

def names(wanted=None):
    out = {}
    with gzip.open(need('names.gz'), 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('A'):
                a, _, n = line.rstrip('\n').partition(' ')
                if wanted is None or a in wanted:
                    out[a] = n
    return out

def search(terms, limit):
    key = ',' + ','.join(str(int(t)) for t in terms) + ','
    hits = []
    with gzip.open(need('stripped.gz'), 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('A') and key in line:
                hits.append(line.split(' ', 1)[0])
    return hits

def update():
    clone = os.path.join(HERE, 'oeisdata')
    if os.path.isdir(clone):                           # the .seq files: from GitHub, always works
        subprocess.run(['git', '-C', clone, 'fetch', '-q', '--depth', '1', 'origin', 'main'], check=True)
        subprocess.run(['git', '-C', clone, 'reset', '-q', '--hard', 'origin/main'], check=True)
        print("oeisdata clone updated to", subprocess.run(['git', '-C', clone, 'log', '-1', '--format=%h %cd'],
              capture_output=True, text=True).stdout.strip())
    for fn in ('stripped.gz', 'names.gz'):             # the bulk files: oeis.org may answer with a bot challenge
        tmp = os.path.join(HERE, fn + '.part')
        r = subprocess.run(['curl', '-sSfL', '--max-time', '600', '-o', tmp, 'https://oeis.org/' + fn])
        try:
            if r.returncode:
                raise OSError(f"curl exit {r.returncode}")
            with gzip.open(tmp) as f:
                f.read()                               # integrity check before replacing the old file
        except Exception as e:
            if os.path.exists(tmp):
                os.remove(tmp)
            print(f"{fn}: download refused ({e}); oeis.org sometimes blocks scripts. Keeping the old file.\n"
                  f"  Download https://oeis.org/{fn} in a browser and save it in {HERE}")
            continue
        os.replace(tmp, os.path.join(HERE, fn))
        print(f"{fn}: {os.path.getsize(os.path.join(HERE, fn)):,} bytes")
    return 0

def full(a):
    clone = os.path.join(HERE, 'oeisdata')
    if not os.path.isdir(clone):
        sys.exit(f"{clone} is missing: run {os.path.join(HERE, 'clone_oeisdata.sh')} first")
    rel = f"seq/{a[:4]}/{a}.seq"
    path = os.path.join(clone, rel)
    if not os.path.exists(path):                       # fetch just this one file into the sparse clone
        sc = os.path.join(clone, '.git', 'info', 'sparse-checkout')
        with open(sc, 'a') as f:
            f.write('/' + rel + '\n')
        subprocess.run(['git', '-C', clone, 'checkout', '-q', 'main'], check=True)
        if not os.path.exists(path):
            print(f"{a}: no such entry in oeisdata"); return 1
    print(open(path, encoding='utf-8', errors='replace').read(), end='')
    return 0

def words(ws, limit):
    ws = [w.lower() for w in ws]
    hits = [(a, n) for a, n in names().items() if all(w in n.lower() for w in ws)]
    print(f"{len(hits)} names contain: {' '.join(ws)}")
    for a, n in hits[:limit]:
        print(f"{a}  {n}")
    if len(hits) > limit:
        print(f"... {len(hits) - limit} more (use -n)")
    return 0

def main(argv):
    limit = 30
    if '-n' in argv:
        i = argv.index('-n'); limit = int(argv[i + 1]); del argv[i:i + 2]
    if argv[:1] == ['--update']:
        return update()
    if argv[:1] == ['-w']:
        return words(argv[1:], limit)
    if len(argv) == 2 and argv[1] == '--full' and re.fullmatch(r'[Aa]\d{1,6}', argv[0]):
        return full('A%06d' % int(argv[0][1:]))
    args = [x for a in argv for x in re.split(r'[,\s]+', a) if x]
    if not args:
        print(__doc__); return 1
    if len(args) == 1 and re.fullmatch(r'[Aa]\d{1,6}', args[0]):
        a = 'A%06d' % int(args[0][1:])
        n = names({a})
        print(f"{a}  {n.get(a, '(not in names.gz)')}")
        return 0
    hits = search(args, limit)
    nm = names(set(hits[:limit]))
    print(f"{len(hits)} sequences contain {','.join(args)}")
    for a in hits[:limit]:
        print(f"{a}  {nm.get(a, '')}")
    if len(hits) > limit:
        print(f"... {len(hits) - limit} more (use -n)")
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
