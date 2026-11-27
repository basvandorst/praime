"""Check the Gamma factor of Sym^2 (and Sym^3, Sym^4) of Delta by the functional-equation test of PARI.
Analytic normalisation: Sym^2: Gamma_R(s+1) Gamma_C(s+11); Sym^3: Gamma_C(s+33/2) Gamma_C(s+11/2);
Sym^4: Gamma_R(s) Gamma_C(s+11) Gamma_C(s+22). PARI uses the arithmetic variable s_a = s + w/2, w = 11k,
weight k_PARI = w + 1, Vga = shifts of Gamma_R in s_a (Gamma_C(s) = Gamma_R(s) Gamma_R(s+1)).
lfuncheckfeq returns log2 of the relative error of the functional equation (very negative = holds)."""
import cypari2
pari = cypari2.Pari()
pari.allocatemem(300 * 10**6)
pari.set_real_precision(int(__import__("sys").argv[1]) if len(__import__("sys").argv) > 1 else 38)

euler = {
    2: "1/((1-p^11*X)*(1-(ramanujantau(p)^2-2*p^11)*X+p^22*X^2))",
    3: "1/((1-ramanujantau(p)*p^11*X+p^33*X^2)*(1-(ramanujantau(p)^3-3*p^11*ramanujantau(p))*X+p^33*X^2))",
    4: "1/((1-p^22*X)*(1-(ramanujantau(p)^2-2*p^11)*p^11*X+p^44*X^2)*(1-((ramanujantau(p)^2-2*p^11)^2-2*p^22)*X+p^44*X^2))",
}
# Vga candidates in the arithmetic variable (w = 11k): correct one first
cands = {
    2: {"Gamma_R(s+1)Gamma_C(s+11) [note]": "[-10,0,1]", "Gamma_R(s)Gamma_C(s+11)": "[-11,0,1]"},
    3: {"Gamma_C(s+33/2)Gamma_C(s+11/2) [note]": "[0,1,-11,-10]", "Gamma_C(s+33/2)Gamma_R(s+11/2)Gamma_R(s+11/2)": "[0,1,-11,-11]"},
    4: {"Gamma_R(s)Gamma_C(s+11)Gamma_C(s+22) [note]": "[-22,-11,-10,0,1]", "Gamma_R(s+1)Gamma_C(s+11)Gamma_C(s+22)": "[-21,-11,-10,0,1]"},
}
for k in (2, 3, 4):
    w = 11 * k
    for name, vga in cands[k].items():
        for eps in (1, -1):
            try:
                L = pari("lfuncreate([n->direuler(p=2,n,%s), 0, %s, %d, 1, %d])" % (euler[k], vga, w + 1, eps))
                e = pari.lfuncheckfeq(L)
                print("Sym^%d  %-48s eps=%+d  lfuncheckfeq = %s" % (k, name, eps, e), flush=True)
            except Exception as ex:
                print("Sym^%d  %-48s eps=%+d  error %s" % (k, name, eps, str(ex)[:80]), flush=True)
