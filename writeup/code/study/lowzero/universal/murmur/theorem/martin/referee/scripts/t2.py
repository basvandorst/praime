import cypari2
pari = cypari2.Pari()
mf = pari.mfinit([15,2],0); B = pari.mfeigenbasis(mf)
print(len(B), pari.mffields(mf), [pari.mfcoef(F,2) for F in B], pari.mfatkineigenvalues(mf,3), pari.mfatkineigenvalues(mf,5))
mf = pari.mfinit([57,2],0); B = pari.mfeigenbasis(mf)
print(len(B), pari.mffields(mf), [pari.mfcoef(F,2) for F in B], pari.mfatkineigenvalues(mf,3), pari.mfatkineigenvalues(mf,19))
