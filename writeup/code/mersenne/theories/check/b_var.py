# variants of B at K=100 to see what reproduces the 'roughly' absolute values
import numpy as np, importlib.util, sys
src = open('b_model.py').read().split('for p_size, K in')[0]
exec(src)
for p_size in (2.6e5, 2.828e5, 3e5):
    print(p_size, {c: round(model(p_size, 100, c)[0], 4) for c in (1, 5, 7, 11)})
