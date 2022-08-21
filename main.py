# Ksenia Fedosova and Kim Klinger-Logan
#
# 26 July 2022 

from bib_fkl import *
from sympy import symbols

# mode; we choose the signs of n1 and sn2
n1_sgn = 1
n2_sgn = 1

#calculateNminus3(n1_sgn, n2_sgn)
#calculateNminus2(n1_sgn, n2_sgn)

alpha = 3/2
beta = 7/2

r = 9

coeff = find_solutions(alpha, beta, r, n1_sgn, n2_sgn)
print_nice_output(coeff, n1_sgn, n2_sgn, "", r, is_solution=True)
