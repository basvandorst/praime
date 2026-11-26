"""kerlib.py: kernels of the S3 cubic-field murmuration and their Poisson-summed primitives.

Kernels (variable u = 4 pi sqrt(rho y)):  D < 0 (Gamma_C(s), weight 1):  K_-(u) = J_0(u);
  D > 0 (Gamma_R(s)^2, even Maass form, spectral parameter 0):  K_+(u) = -Y_0(u) + (2/pi) K_0(u).
  Both are fixed by  int_0^oo u^s K(u) du = (4 pi)^s gamma(1/2 + s/2) / gamma(1/2 - s/2),  gamma the gamma factor:
  J_0: 2^s G((1+s)/2)/G((1-s)/2);  K_+: 4^s G((1+s)/4)^2/G((1-s)/4)^2 = 2^s G((1+s)/2)^2 (1 - sin(pi s/2))/pi.
Squares summed by Poisson (S(Z) = sum_{n>=1} K(4 pi n sqrt Z)) and integrated, C(Z) = int_0^Z S:
  C_-(Z) = sqrt(Z)/(2 pi) - Z/2 + (1/2pi) sum_{1<=k<2 sqrt Z} sqrt(4Z - k^2)
  C_+(Z) = sqrt(Z)/(2 pi) + (1/pi) [Z ln Z - Z + L(Z)],
           L(Z) = lim_N [ sum_{k<=N} (sqrt(k^2+4Z) - sqrt((k^2-4Z)_+))/2 - 2 Z ln N ]
  large Z:  C_-(Z) ~ (1/pi) Z^(1/4) zeta(-1/2, th) - zeta(-3/2, th)/(8 pi Z^(1/4)),  th = frac(2 sqrt Z)
            C_+(Z) ~ -(1/pi) Z^(1/4) zeta(-1/2, 1-th) - zeta(-3/2, 1-th)/(8 pi Z^(1/4))
"""
import math, numpy as np
import mpmath as mp

EULER = 0.57721566490153286061
ZETA5 = 1.0369277551433699263
ZETA9 = 1.0020083928260822144

_TH = np.linspace(0.0, 1.0, 4001)
def _hz_table(s):  # zeta(s, th + 1) on [0, 1]: smooth
    return np.array([float(mp.zeta(s, 1 + float(t))) for t in _TH])
_H12 = None; _H32 = None
def _tables():
    global _H12, _H32
    if _H12 is None:
        _H12 = _hz_table(-0.5); _H32 = _hz_table(-1.5)
    return _H12, _H32

def hurwitz(s, th):
    """zeta(s, th) for th in (0, 1], s in {-1/2, -3/2}: th^(-s) + zeta(s, th + 1) (table, cubic-accurate interpolation)."""
    H12, H32 = _tables(); tab = H12 if s == -0.5 else H32
    return th ** (-s) + np.interp(th, _TH, tab)

def C_minus_exact(Z):
    Z = np.asarray(Z, float); out = np.sqrt(Z) / (2 * math.pi) - Z / 2
    kmax = int(2 * math.sqrt(Z.max())) + 1
    acc = np.zeros_like(Z)
    for k in range(1, kmax + 1):
        acc += np.sqrt(np.maximum(4 * Z - k * k, 0.0))
    return out + acc / (2 * math.pi)

def L_exact(Z):
    Z = np.asarray(Z, float); N = int(10 * math.sqrt(max(Z.max(), 1.0))) + 40
    acc = 2 * EULER * Z
    for k in range(1, N + 1):
        acc += (np.sqrt(k * k + 4 * Z) - np.sqrt(np.maximum(k * k - 4 * Z, 0.0))) / 2 - 2 * Z / k
    # tail k > N: (sqrt(k^2+4Z) - sqrt(k^2-4Z))/2 - 2Z/k = 4 Z^3/k^5 + 28 Z^5/k^9 + ...
    acc += 4 * Z**3 * float(mp.zeta(5, N + 1)) + 28 * Z**5 * float(mp.zeta(9, N + 1))
    return acc

def C_plus_exact(Z):
    Z = np.asarray(Z, float); Zs = np.maximum(Z, 1e-300)
    return np.sqrt(Z) / (2 * math.pi) + (Z * np.log(Zs) - Z + L_exact(Z)) / math.pi

def C_asym(Z, sgn):
    Z = np.asarray(Z, float); th = np.mod(2 * np.sqrt(Z), 1.0); q = Z ** 0.25
    if sgn < 0:
        th = np.where(th == 0.0, 1.0, th)
        return q * hurwitz(-0.5, th) / math.pi - hurwitz(-1.5, th) / (8 * math.pi * q)
    t2 = 1.0 - th
    return -q * hurwitz(-0.5, t2) / math.pi - hurwitz(-1.5, t2) / (8 * math.pi * q)

ZSWITCH = 10.0
def C(Z, sgn):
    """C_sgn(Z) for an array Z >= 0: exact form below ZSWITCH, asymptotic above."""
    Z = np.asarray(Z, float); out = np.empty_like(Z)
    lo = Z < ZSWITCH; hi = ~lo
    if lo.any(): out[lo] = (C_minus_exact if sgn < 0 else C_plus_exact)(Z[lo])
    if hi.any(): out[hi] = C_asym(Z[hi], sgn)
    return out

def C_small(Z, sgn):
    """Expansion for Z < 1/4 (exact for sgn < 0; through Z^5 for sgn > 0, k = 1 term exact)."""
    Z = np.asarray(Z, float)
    if sgn < 0: return np.sqrt(Z) / (2 * math.pi) - Z / 2
    Lz = 2 * EULER * Z + (np.sqrt(1 + 4 * Z) - np.sqrt(1 - 4 * Z)) / 2 - 2 * Z + 4 * Z**3 * (ZETA5 - 1) + 28 * Z**5 * (ZETA9 - 1)
    return np.sqrt(Z) / (2 * math.pi) + (Z * np.log(Z) - Z + Lz) / math.pi
