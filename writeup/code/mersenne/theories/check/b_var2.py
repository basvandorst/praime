import numpy as np
src = open('b_model.py').read().split('for p_size, K in')[0]
exec(src)
for p_size in (3.5e5, 4e5, 4.5e5):
    print(p_size, {c: round(model(p_size, 100, c)[0], 4) for c in (1, 5, 7, 11)})
