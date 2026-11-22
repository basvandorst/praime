"""Compare own sieve with the study's common.primes_mask on [1e15, 1e15+1e7) (the study function is only the object under test)."""
import sys, numpy as np
sys.path.insert(0, '/home/bas/praime/writeup/code/study'); import common
sys.path.insert(0, '.'); from sieve import primes_in
L = 10**15; W = 10**7
A = primes_in(L, W); B = np.nonzero(common.primes_mask(L, W))[0].astype(np.int64) + L
print(len(A), len(B), np.array_equal(A, B), 'twins', int(np.sum(np.diff(A) == 2)))
