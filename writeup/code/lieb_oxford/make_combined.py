# Build the single self-contained verification script lo_certify.py from the two parts.
import sys, json, textwrap

d = '/tmp/claude-1000/-home-bas-praime/9968315b-c905-4858-b310-2bb45fd31cf2/scratchpad/lo_cert/'
ball = open(d + 'lo_ball_certify.py').read()
sph = open(d + 'cert_spheres.py').read()
meas = json.load(open(d + sys.argv[1]))
meas.pop('source', None)

# parameters become fixed constants
ball = ball.replace("RHO = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-3", "RHO = 2.5e-4")
ball = ball.replace("print(f'{NCHK[0]} checks, {NCHK[1]} failed')", "print(f'part A: {NCHK[0]} checks, {NCHK[1]} failed')\nreturn NCHK")
sph = sph.replace("inp = json.load(open(sys.argv[1]))", "inp = MEASURE")
sph = sph.replace("RHO_HI = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-3", "RHO_HI = 1e-3")
sph = sph.replace("RHO_LO = float(sys.argv[3]) if len(sys.argv) > 3 else 3e-3", "RHO_LO = 3e-3")
sph = sph.replace("RTAIL = int(sys.argv[4]) if len(sys.argv) > 4 else 200", "RTAIL = 600")
sph = sph.replace("print(f'{NCHK[0]} checks, {NCHK[1]} failed')", "print(f'part B: {NCHK[0]} checks, {NCHK[1]} failed')\nreturn NCHK")
for name, src in (('ball', ball), ('spheres', sph)):
    assert 'sys.argv' not in src.split('"""', 2)[2], name


def body(src):
    # drop the module docstring, indent the rest into a function body
    doc_end = src.index('"""', src.index('"""') + 3) + 3
    return src[:doc_end], textwrap.indent(src[doc_end:], '    ')


docA, bodyA = body(ball)
docB, bodyB = body(sph)
out = ['#!/usr/bin/env python3',
       '"""Certified upper bounds for the Lieb-Oxford constant (one self-contained verification script).',
       '',
       'Part A: mu = nu = uniform unit ball, f = g (LLS 2022, Lemma 10)  ->  c_LO <= 1.604335.',
       'Part B: mu = mixture of 48 spheres, nu = one sphere (exact rational data below), f = g  ->  c_LO <= 1.597856.',
       'Run: python3 lo_certify.py      (about 7 minutes, < 1 GB)',
       'Docstrings of the two parts follow.',
       '"""',
       'MEASURE = ' + json.dumps(meas),
       '',
       '',
       'def part_A():',
       textwrap.indent(docA, '    '),
       bodyA,
       '',
       'def part_B():',
       textwrap.indent(docB, '    '),
       bodyB,
       '',
       "if __name__ == '__main__':",
       '    nA = part_A()',
       '    nB = part_B()',
       "    print(f'{nA[0] + nB[0]} checks, {nA[1] + nB[1]} failed')",
       '']
open(d + 'lo_certify.py', 'w').write('\n'.join(out))
print('written')
