import sympy as smp


l,E = smp.symbols('l E')
r = smp.Function('r')
n = smp.Symbol('n',integer=True)
r_1 = -2*E
r_0 = 1

def hydrogen(n,rn_2,rn_3):
    return 1/(8*E*n)*((1-n)*(n*(n-2)-4*l*(l+1))*rn_3-4*(2*n-1)*rn_2) # return rn_1


rn_2 = r_0
rn_3 = r_1
for i in range(2,10):
    rp = hydrogen(i, rn_2, rn_3)
    rn_2,rn_3 = rp,rn_2
    smp.print_latex(rp.simplify())

