"""Class-response law (derivation.md section 6), general version.

Short-range model M_Q (Q = 30030): every residue coprime to Q is 'prime' independently with probability
p = (Q/phi(Q))/lam. In M_Q the class counts N_c (c mod Q) are independent, Var = (1-p) E N_c, and by periodicity the
linear part of any window statistic is its regression on the N_c. The law replaces that class covariance by the
Hardy-Littlewood one (closed form, section 4):

  E chi^2(primes) = E chi^2(M_Q) - FO_sieve + FO_HL,   FO_X = sum_cells (1/E_cell) r_cell^T Cov_X r_cell,

r_cell[c] = response of the cell residual to one more prime at a site of class c, measured by insertion.
E chi^2(M_Q) = block mean of chi^2 (block against the pooled table) x NB/(NB-1).

Statistics (cells), all row-conditioned except 'band':
  band:k        D_k = p_{n+k}-p_n in bands [lo, lo+10], lo = 2,14,..,230 (D <= 240), expected N_block * pooled share
  trans:q:k     classes (p_n mod q, p_{n+k} mod q), expected N_a * pooled P(b|a), pair kept if D <= DMAX
  resgap        rows p_n mod 210, columns min(g/2, 20) (gaps 2..38, >= 40), cells with pooled share 0 dropped
usage: simlaw.py x H_for_v B_blocks_in_design DMAX nchunks chunk_span block_span seed out.json stat [stat ...]
 H_for_v: window length entering v_Q; B_blocks_in_design: number of blocks of the registered design (adds
 ln B/((B-1) lam) to v_Q, the block-minus-pooled correction), 1 for a single window against a fixed law."""
import sys, math, json, time
import numpy as np

X = float(sys.argv[1]); HV = float(sys.argv[2]); BDES = int(sys.argv[3]); DMAX = int(sys.argv[4])
NCH = int(sys.argv[5]); CHUNK = int(float(sys.argv[6])); BSPAN = int(float(sys.argv[7])); SEED = int(sys.argv[8])
OUT = sys.argv[9]; SPECS = sys.argv[10:]
LAM = math.log(X)
Q = 30030; PQ = (2, 3, 5, 7, 11, 13)
adm = np.array([a for a in range(Q) if math.gcd(a, Q) == 1], np.int64); NC = len(adm)
cidx = np.full(Q, -1, np.int64); cidx[adm] = np.arange(NC)
p = (Q / NC) / LAM
EG = 0.5772156649015329; BMS = 1 - EG - math.log(2 * math.pi)
SQ = sum(math.log(r) / (r - 1) for r in PQ)
vQ = 1 - (math.log(HV / Q) + SQ - BMS) / LAM + (math.log(BDES) / ((BDES - 1) * LAM) if BDES > 1 else 0.0)
rng = np.random.default_rng(SEED)


class Stat:
    def __init__(self, spec):
        f = spec.split(':'); self.spec = spec; self.kind = f[0]
        if self.kind == 'band':
            self.k = int(f[1]); self.ncell = 20; self.rows = None
        elif self.kind == 'trans':
            self.q = int(f[1]); self.k = int(f[2])
            red = [a for a in range(self.q) if math.gcd(a, self.q) == 1]
            self.lut = np.full(self.q, -1); self.lut[red] = np.arange(len(red)); self.R = len(red)
            self.ncell = self.R * self.R
        elif self.kind == 'resgap':
            self.k = 1; red = [a for a in range(210) if math.gcd(a, 210) == 1]
            self.lut = np.full(210, -1); self.lut[red] = np.arange(48); self.R = 48; self.ncell = 48 * 20
        self.counts = []
        self.acc = [np.zeros((NC, self.ncell)), np.zeros((NC, self.ncell))]

    def cell(self, left, right):
        d = right - left
        if self.kind == 'band':
            c = (d - 2) // 12
            return np.where((d <= 240) & (c >= 0) & (c < 20), c, -1)
        if self.kind == 'trans':
            return np.where(d <= DMAX, self.lut[left % self.q] * self.R + self.lut[right % self.q], -1)
        g = np.minimum(d // 2, 20) - 1
        return self.lut[left % 210] * 20 + g

    def row_of_cell(self):
        return np.arange(self.ncell) // (self.R if self.kind == 'trans' else 20)


def insertion(st, P, M, h):
    k = st.k
    j = np.searchsorted(P, M) - 1
    ok = (j - k + 1 >= 0) & (j + k < len(P))
    M = M[ok]; j = j[ok]; cl = cidx[M % Q]
    def put(left, right, w):
        c = st.cell(left, right); sel = c >= 0
        st.acc[h] += w * np.bincount(cl[sel] * st.ncell + c[sel], minlength=NC * st.ncell).reshape(NC, st.ncell)
    for t in range(k):
        i = j - t
        put(P[i], P[i + k], -1.0)
    put(P[j - k + 1], M, 1.0); put(M, P[j + k], 1.0)
    for off in range(k - 1):
        i = j - k + 2 + off
        put(P[i], P[i + k - 1], 1.0)


def main():
    t0 = time.time()
    stats = [Stat(s) for s in SPECS]
    nsite = [np.zeros(NC), np.zeros(NC)]; periods = 0; nblk = 0
    for ch in range(NCH):
        base = (np.arange(CHUNK // Q, dtype=np.int64)[:, None] * Q + adm[None, :]).ravel() + Q
        isp = rng.random(len(base)) < p
        P = base[isp]; M = base[~isp]; del base, isp
        periods += CHUNK // Q; h = ch % 2
        nsite[h] += np.bincount(cidx[M % Q], minlength=NC)
        edges = np.arange(Q, Q + (CHUNK // Q) * Q + 1, BSPAN)
        for st in stats:
            c = st.cell(P[:-st.k], P[st.k:])
            blk = np.searchsorted(edges, P[:-st.k], side='right') - 1
            sel = (c >= 0) & (blk < len(edges) - 1)
            nb = len(edges) - 1
            st.counts.append(np.bincount(blk[sel] * st.ncell + c[sel], minlength=nb * st.ncell).reshape(nb, st.ncell))
            insertion(st, P, M, h)
        nblk += len(edges) - 1
        print('chunk', ch, '%.0f s' % (time.time() - t0), flush=True)
    fib = {r: (adm % (Q // r)) for r in PQ}
    res = {'x': X, 'lam': LAM, 'p': p, 'H_for_v': HV, 'B_design': BDES, 'v_Q': vQ, 'blocks_sim': nblk,
           'block_span': BSPAN, 'stats': {}}
    for st in stats:
        C = np.concatenate(st.counts).astype(float); NBk = C.shape[0]
        pooled = C.sum(0)
        if st.kind == 'band':
            Pc = pooled / pooled.sum(); E = C.sum(1, keepdims=True) * Pc[None, :]
            keep = (Pc * 4.8e6 >= 50) if st.k == 1 else np.ones(20, bool)
            nrows = 1
        else:
            rows = st.row_of_cell(); nr = st.R
            rowtot = np.zeros((NBk, nr)); np.add.at(rowtot.T, rows, C.T)
            prow = np.bincount(rows, weights=pooled, minlength=nr)
            Pc = pooled / prow[rows]
            E = rowtot[:, rows] * Pc[None, :]
            keep = Pc > 0
            nrows = nr
        chi = (((C - E) ** 2 / np.where(E > 0, E, 1))[:, keep]).sum(1)
        dshort = float(chi.mean() * NBk / (NBk - 1)); dsd = float(chi.std(ddof=1))
        r = []
        for h in (0, 1):
            psi = st.acc[h] / nsite[h][:, None]
            if st.kind == 'band':
                rr = psi - psi.sum(1, keepdims=True) * Pc[None, :]
            else:
                rt = np.zeros((NC, nrows)); np.add.at(rt.T, rows, psi.T)
                rr = psi - rt[:, rows] * Pc[None, :]
            r.append(rr[:, keep])
        Ecell = (pooled / periods)[keep]
        FO_s = float(np.sum(p * (1 - p) * (r[0] * r[1]).sum(0) / Ecell))
        cross = np.zeros_like(r[1])
        for q_ in PQ:
            fs = np.zeros((Q // q_, r[1].shape[1])); np.add.at(fs, fib[q_], r[1])
            cross += (math.log(q_) / ((q_ - 1) * LAM)) * (fs[fib[q_]] - r[1])
        FO_hl = float(np.sum(p * (r[0] * (vQ * r[1] - cross)).sum(0) / Ecell))
        dof = int(keep.sum()) - nrows
        colpart = None
        if st.kind == 'trans':                                   # realized column part per block (sieve model)
            R_ = st.R; cols = []
            Pm = Pc.reshape(R_, R_)
            for bi in range(NBk):
                eps = (C[bi] - E[bi]).reshape(R_, R_); cvec = eps.sum(0)
                rt_ = rowtot[bi]
                Mb = sum(rt_[a] * (np.diag(Pm[a]) - np.outer(Pm[a], Pm[a])) for a in range(R_))
                cols.append(float(cvec @ np.linalg.pinv(Mb) @ cvec))
            colpart = float(np.mean(cols) * NBk / (NBk - 1))
        res['stats'][st.spec] = {'chi_blocks': chi.tolist(), 'col_part_sieve_per_block': colpart, 'P_pooled': Pc.tolist(),'dof': dof, 'chi2_sieve_model': dshort, 'chi2_sd_per_block': dsd,
                                 'chi2_sieve_err': dsd / math.sqrt(NBk), 'FO_sieve': FO_s, 'FO_HL': FO_hl,
                                 'chi2_law': dshort - FO_s + FO_hl, 'ratio_law': (dshort - FO_s + FO_hl) / dof,
                                 'ratio_sieve': dshort / dof}
        print(st.spec, {a: round(v, 4) if isinstance(v, float) else v for a, v in res['stats'][st.spec].items() if a not in ('chi_blocks', 'P_pooled')}, flush=True)
    json.dump(res, open(OUT, 'w'), indent=1)
    print('done %.0f s' % (time.time() - t0))


if __name__ == '__main__':
    main()
