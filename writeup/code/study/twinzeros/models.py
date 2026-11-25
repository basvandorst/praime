"""models.py -- regressors and predicted coefficients for each pattern (shared by the discovery scripts and the test)."""
import numpy as np
import tzlib, theory

QS = (3, 4, 5, 7, 8, 9, 11, 12, 13, 15, 17, 19, 20, 21, 23, 24, 28, 29)
QS2 = (3, 4, 5, 7, 11, 12, 13, 15, 20, 21)
PAT_KEY = {"h2": "P2", "h4": "P4", "h6": "P6", "h8": "P8", "h10": "P10", "h12": "P12", "h30": "P30",
           "SG": "SG", "T026": "T026", "T046": "T046"}
# characters predicted ABSENT for the translation patterns, used for the null tests (one per conjugate pair; Re and Im)
NULL_SET = ["3.1", "4.1", "5.2", "7.1", "7.4", "8.1", "8.2", "11.1", "12.1", "9.1"]


class Data:
    def __init__(self, paths):
        k, X, cells = tzlib.load_cells(paths)
        self.k, self.X = k, X
        self.cum = tzlib.cumulative(cells)
        self.Xf = X.astype(float)
        self.u = np.log(self.Xf)
        self.sq = np.sqrt(self.Xf)
        self.chars = theory.primitive_characters(30)
        self.byl = {c["label"]: c for c in self.chars}
        self._cache = {}

    def J(self, x, right=False):
        return int(np.searchsorted(self.X, x, side="right" if right else "left"))

    def Yzeta(self, a=1):
        if a == 1:
            return (self.cum["th"] - self.Xf) / self.sq
        return (self.cum["th_2"] - 2 * self.Xf) / self.sq

    def Ychi(self, lab, a=1):
        key = (lab, a)
        if key not in self._cache:
            c = self.byl[lab]
            self._cache[key] = tzlib.char_series(self.cum, c["v"], c["q"], doubled=(a == 2)) / self.sq
        return self._cache[key]

    def Ypat(self, name):
        S = theory.singular_const(name)
        return (self.cum[PAT_KEY[name]] - S * self.Xf) / self.sq, S


def chi_combination(D, rows, labels_ok, a):
    """Real series sum over present characters (one per conjugate pair) of the predicted contribution at scale a."""
    out = 0.0
    seen = set()
    for r in rows:
        if r["label"] in seen or not r["present"] or abs(r["beta"].get(a, 0)) < 1e-12:
            continue
        if r["q"] not in labels_ok:
            continue
        y = D.Ychi(r["label"], a)
        if r["real"]:
            out = out + r["beta"][a].real * y.real
        else:
            out = out + 2 * np.real(r["beta"][a] * y)
            seen.add(r["conj"])
        seen.add(r["label"])
    return out


def build(D, name, null_set=NULL_SET, per_conductor=False):
    """Regressors (list of (label, series, predicted coefficient)) for the pattern."""
    S, rows = theory.prediction_table(name, D.chars)
    forms = theory.PATTERNS[name]
    scales = theory.scales(name)
    regs = []
    for a in scales:
        n_a = sum(1 for a2, b in forms if a2 == a)
        regs.append(("zeta(a=%d)" % a, D.Yzeta(a), S * n_a / a))
    qs_ok = {a: (QS if a == 1 else QS2) for a in scales}
    if per_conductor:
        groups = sorted(set(r["q"] for r in rows if r["present"]))
        for q in groups:
            ser = 0.0
            for a in scales:
                if q in qs_ok[a]:
                    ser = ser + chi_combination(D, [r for r in rows if r["q"] == q], qs_ok[a], a)
            if not np.isscalar(ser):
                regs.append(("chi q=%d (loudness)" % q, ser, 1.0))
    else:
        ser = 0.0
        for a in scales:
            ser = ser + chi_combination(D, rows, qs_ok[a], a)
        regs.append(("all present chi (loudness)", ser, 1.0))
    byl = {r["label"]: r for r in rows}
    for lab in null_set:
        r = byl[lab]
        if r["present"]:
            continue
        y = D.Ychi(lab, 1)
        if r["real"]:
            regs.append(("null %s" % lab, y.real, 0.0))
        else:
            regs.append(("null Re %s" % lab, y.real, 0.0))
            regs.append(("null Im %s" % lab, y.imag, 0.0))
    return S, regs


def kforms(name):
    return len(theory.PATTERNS[name])
