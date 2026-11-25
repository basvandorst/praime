"""Load the repo's frozen libraries by path, unmodified: bklib.py (two-point function, Bogomolny-Keating 2013) and
tplib.py (three-point function, ratios conjecture). tplib.py reads grid.npz and gtab.npz from the directory of its
__file__, so its source is executed with __file__ set to this folder, where grid.npz is a byte copy of the repo's
amp/grid.npz (sha256 7af1eb14...) and gtab.npz is produced by the repo's gtab.py (executed the same way)."""
import os, sys, types, hashlib, importlib.util
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', '..', '..'))
LOW = REPO + '/writeup/code/study/lowzero'
HERE = os.path.dirname(os.path.abspath(__file__))
SHA = {'bklib.py': '9bfc1a3bf6c3b1bf006118660e7596eaab92be41b77d8eb33565ce1bdfa66161',
       'tplib.py': 'bd049701c58a938016d229421d63d82e63fc2ab3ffc2d824920a9ac076cc4f62',
       'gtab.py': 'ce2468ea8f451fa4da65c622a97355fb5685b4eb5d04f8e4f3bdd405abcea07d',
       'grid.npz': '7af1eb149cfca33513879667465b93f66b9d7b62a5925969ffb4106547685cd1'}


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def load_bklib():
    p = LOW + '/bk/bklib.py'; assert sha(p) == SHA['bklib.py']
    if 'bklib' in sys.modules: return sys.modules['bklib']
    spec = importlib.util.spec_from_file_location('bklib', p); m = importlib.util.module_from_spec(spec)
    sys.modules['bklib'] = m; spec.loader.exec_module(m); return m


def make_gtab():
    src = LOW + '/triple/lib/gtab.py'; assert sha(src) == SHA['gtab.py']
    if os.path.exists(HERE + '/gtab.npz'): return
    code = open(src).read(); g = {'__file__': HERE + '/gtab.py', '__name__': 'gtab'}
    exec(compile(code, src, 'exec'), g)


def load_tplib():
    load_bklib(); make_gtab()
    assert sha(HERE + '/grid.npz') == SHA['grid.npz']
    src = LOW + '/triple/lib/tplib.py'; assert sha(src) == SHA['tplib.py']
    if 'tplib' in sys.modules: return sys.modules['tplib']
    m = types.ModuleType('tplib'); m.__file__ = HERE + '/tplib.py'
    sys.modules['tplib'] = m
    exec(compile(open(src).read(), src, 'exec'), m.__dict__)
    return m
