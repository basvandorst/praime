import numpy as np, os
H = os.path.dirname(os.path.abspath(__file__))
for k in (2, 4):
    a = np.load(H + "/r159_exact_k%d.npz" % k); w = a["w"]
    m = {q: a[q] / w for q in ("V", "sm", "cnx", "cnp", "R", "Rp")}
    # empirical bin se of V and of R (about bin mean), and of the CN error, approx via second moments
    seV = np.sqrt(a["V2"] - w ** 2 * 0) / w  # raw second moment (mean small vs scatter)
    seR = np.sqrt(a["R2"]) / w; sed = np.sqrt(a["dcn2"]) / w
    lo = slice(0, 5); hi = slice(5, 60)
    print("k=%d: identity V = sm + cnx + R exact by construction" % k)
    print("  noise per bin (median over bins): se(V) %.4f  se(R) %.4f  rms(cnx - cnp) %.4f" % (np.median(seV), np.median(seR), np.median(sed)))
    zc = (m["cnx"] - m["cnp"]) / seV; zr = (m["R"] - m["Rp"]) / seV
    print("  chi2 (in data errors) of exact CN vs nu-average, bins y>1/4: %.1f / 55;  of R vs registered A2chi term: all %.1f / 60, y<1/4 %.2f / 5" % (
        np.sum(zc[hi] ** 2), np.sum(zr ** 2), np.sum(zr[lo] ** 2)))
    # class-number average vs nu at its own (much smaller) noise
    zc2 = (m["cnx"] - m["cnp"]) / np.where(sed > 0, sed, np.inf)
    print("  exact CN / nu-average, pooled y>1/4 (weighted): %.5f ; chi2 in own CN noise %.1f / 55" % (
        np.sum(m["cnx"][hi] * m["cnp"][hi] / sed[hi] ** 2) / np.sum(m["cnp"][hi] ** 2 / sed[hi] ** 2), np.sum(zc2[hi] ** 2)))
    # amplitude of R against its registered smooth term (k=2), in data errors
    if k == 2:
        A = np.sum(m["R"] * m["Rp"] / seV ** 2) / np.sum(m["Rp"] ** 2 / seV ** 2); sA = 1 / np.sqrt(np.sum(m["Rp"] ** 2 / seV ** 2))
        print("  R = a * 24 rho y^1.5: a = %.4f +- %.4f" % (A, sA))
    sel = [1, 4, 6, 10, 20, 30, 40, 50, 59]
    for i in sel:
        print("   y=%.3f  V %.3f  sm %.3f  cnx %.3f cnp %.3f  R %.3f Rp %.3f  se(V) %.3f" % ((i + .5) * .05, m["V"][i], m["sm"][i], m["cnx"][i], m["cnp"][i], m["R"][i], m["Rp"][i], seV[i]))
