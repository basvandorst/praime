"""The phase rule for bands of zeros (RIEMANN.md section 12).

After demodulation, prime power m carries the phase e^{-ic ln m}. Around n,
  c ln(m/n) = c (m - n)/n - c (m - n)^2 / (2 n^2) + ...
For candidates s apart (s = 2 for odd numbers) the step between neighbours is c s / n (mod 2 pi), and
across a window of K candidates the chirp is c (K s)^2 / (2 n^2). When the chirp is small and the
step is a multiple of pi, the phases of a window are real up to one common factor: the real and the
imaginary part of the samples carry the same real signal, and the band gives half its information.
The bad points are n_j = s c / (pi j), j = 1, 2, ...

usage: python3 band_phase.py LOG   (LOG: output lines of band_decoder.py; single bands with c > 0, odd candidates)
The bad points are listed when the chirp over the window is below 100 rad; above that the phases are random.
"""
import math
import re
import sys

S, K = 2, 12
fracs = []

for line in open(sys.argv[1]):
    m = re.match(r'1 band\(s\) at centre ([\d.e+]+), total width (\d+).*wheel \(2,\):.*first wrong \[(\d+)', line)
    if not m or float(m.group(1)) == 0:
        continue
    c, W, n = float(m.group(1)), int(m.group(2)), int(m.group(3))
    chirp = c * (K * S) ** 2 / (2 * n * n)
    head = f"centre {c:.6g}, width {W}: reach {n} = {n / W:.2f} W; "
    if chirp < 100:
        j = S * c / (math.pi * n)
        jn = round(j)
        bad = S * c / (math.pi * jn)
        fracs.append(abs(n - bad) / (math.pi * n * n / (S * c) / 2))
        print(head + f"step {j:.2f} pi; nearest bad point (j = {jn}) {bad:.0f}, reach {n - bad:+.0f} from it; "
              f"bad points {math.pi * n * n / (S * c):.0f} apart; chirp over the window {chirp:.2f} rad")
    else:
        print(head + f"chirp over the window {chirp:.1e} rad: random phases")
if fracs:
    from scipy.stats import chi2
    stat = -2 * sum(math.log(f) for f in fracs)
    print(f"{len(fracs)} bands with a small chirp; distance to the nearest bad point in half spacings: "
          + ', '.join(f'{f:.3f}' for f in fracs)
          + f". If the failures did not care about bad points these would be uniform on [0, 1]; "
            f"Fisher's combined test: chi2 = {stat:.1f} with {2 * len(fracs)} degrees of freedom, p = {chi2.sf(stat, 2 * len(fracs)):.1e}")
