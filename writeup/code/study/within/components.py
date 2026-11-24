"""Stage 1 of law.py: the Hardy-Littlewood components of the within part at one height lam (cached as JSON).
usage: components.py LAM [LAM ...]  ->  data/comp_<lam>.json"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, math, json, time
from chainlaw import Model2
from triple import D3_model, tau1_from_D3
from tailmod import autocorr, tail_values, responses_all
from hlcore import S_tail_table
HS = [10 ** 4, 10 ** 5, 10 ** 6]; QS = (3, 5, 7); J = 60
for lam in [float(a) for a in sys.argv[1:]]:
    t0 = time.time(); out = {'lam': lam}
    # (1) lag-1 term from the HL triple law, y = 13 and y = 11 (truncation check)
    for y in (13, 11):
        D3s, n = D3_model(y, lam, QS, T1=360, T2=560)
        out['tau1_triple_y%d' % y] = {q: tau1_from_D3(D3s[q])[0] for q in QS}
    # (2) fine-class chain (HL pair factor to the predecessor), k = 1 and k = 2, y = 13; tails (w2) at y = 13, 11
    M = Model2(lam, y=13)
    for q in QS:
        tau1, s1 = M.chain_tau(q, J=J); tau2, s2 = M.chain_tau_k2(q, J=J)
        A1 = autocorr(M.Q, M.U, M.responses(q, s1)); A2 = autocorr(M.Q, M.U, M.responses_k2(q, s2))
        C = s1['C']
        t1 = tail_values(M.Q, M.U, A1, C, lam, 13, HS, Stab=M.Stab)
        t1w3 = tail_values(M.Q, M.U, A1, C, lam, 13, HS, Stab=M.Stab, weight='w3')
        t2 = tail_values(M.Q, M.U, A2, C, lam, 13, HS, Stab=M.Stab)
        Pm = s1['Pm']; Bq = C * np.mean(np.diag(Pm)) - 1
        out[q] = {'chain_k1': tau1.tolist(), 'chain_k2': tau2.tolist(), 'tail_k1': t1['W_tail'], 'tail_k1_w3': t1w3['W_tail'],
                  'tail_k2': t2['W_tail'], 'V_k1': t1['V_per_dir'], 'V_k2': t2['V_per_dir'], 'Bq': float(Bq)}
    Q11, U11, R11 = responses_all(11, lam, QS); St11 = S_tail_table(10 ** 6, 11)
    for q in QS:
        A = autocorr(Q11, U11, R11[q]['r'])
        out[q]['tail_k1_y11'] = tail_values(Q11, U11, A, R11[q]['C'], lam, 11, HS, Stab=St11)['W_tail']
    # (3) the sieve chain M_Q (S_{>y} == 1), y = 13
    M0 = Model2(lam, y=13, sieve_only=True)
    for q in QS:
        out[q]['MQ_k1'] = M0.chain_tau(q, J=J)[0].tolist(); out[q]['MQ_k2'] = M0.chain_tau_k2(q, J=J)[0].tolist()
    json.dump(out, open('data/comp_%.4f.json' % lam, 'w'))
    print('lam %.4f done in %.0fs' % (lam, time.time() - t0), flush=True)
