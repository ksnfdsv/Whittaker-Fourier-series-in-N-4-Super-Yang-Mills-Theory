from sympy import Matrix, symbols
from sympy import linsolve, simplify, fraction, factor, expand, cancel
import datetime
from sympy import mathematica_code as mcode

n1, n2, absn1, absn2, y, EulerGamma, C, l1, l2, pi, x, b0, b1, p0, p1 = symbols(
    'n1, n2, absn1, absn2, y, EulerGamma, C, l1, l2, pi, x, b0, b1, p0, p1')
alpha3, alpha5, alpha7, alpha9, beta5, beta7, beta9, gamma5, gamma7, gamma9 = symbols(
    'alpha3, alpha5, alpha7, alpha9, beta5, beta7, beta9, gamma5, gamma7, gamma9')
d4, d6 = symbols(
    'd4, d6')

# We will use b0 as a symbol to denote BesselK(0, n_1 x)
#             b1                       BesselK(1, n_1 x)
# We will use p0 as a symbol to denote BesselK(0, n_2 x)
#             p1                       BesselK(1, n_2 x)

# We recursively define BesselK(k, n_1 x) from BesselK(0, n_1 x) and BesselK(1, n_1 x) 
# (this is not the best solution time-wise, but other calculations are way more 
# time-consuming; roughly speaking, if we need too much time for this recursion, we would
# have to wait for eternity to solve corresponding linear equation)

def bk1(k, t):
    if k == 1:
        return b1
    if k == 0:
        return b0
    return bk1(k - 2, t) + 2 * (k - 1) * bk1(k - 1, t) / t


# We recursively define BesselK(k, n_2 x) from BesselK(0, n_2 x) and BesselK(1, n_2 x)
def bk2(k, t):
    if k == 1:
        return p1
    if k == 0:
        return p0
    return bk2(k - 2, t) + 2 * (k - 1) * bk2(k - 1, t) / t


# we initialize the array that represents the differential operator 
def init_arr(n1_sign, n2_sign, min_n, max_n, lambd):
    arr = [[0 for _ in range(4 * (max_n - min_n + 3))] for _ in range(4 * (max_n - min_n + 1))]
    for j in range(max_n - min_n + 1):
        p = 4 * j
        degree = j + min_n

        arr[p][p] = degree * degree - degree - lambd
        arr[p][p + 5] = absn2 - 2 * degree * absn2
        arr[p][p + 6] = absn1 - 2 * degree * absn1
        arr[p][p + 8] = - n1_sign * n2_sign * 2 * absn1 * absn2
        arr[p][p + 11] = 2 * absn1 * absn2

        arr[p + 1][p + 1] = 2 + degree * degree - 3 * degree - lambd
        arr[p + 1][p + 4] = absn2 - 2 * degree * absn2
        arr[p + 1][p + 7] = 3 * absn1 - 2 * degree * absn1
        arr[p + 1][p + 9] = - n1_sign * n2_sign * 2 * absn1 * absn2
        arr[p + 1][p + 10] = 2 * absn1 * absn2

        arr[p + 2][p + 2] = 2 - 3 * degree + degree * degree - lambd
        arr[p + 2][p + 4] = absn1 - 2 * degree * absn1
        arr[p + 2][p + 7] = 3 * absn2 - 2 * degree * absn2
        arr[p + 2][p + 9] = 2 * absn1 * absn2
        arr[p + 2][p + 10] = -2 * n1_sign * n2_sign * absn1 * absn2

        arr[p + 3][p + 3] = 6 - 5 * degree + degree * degree - lambd
        arr[p + 3][p + 5] = 3 * absn1 - 2 * degree * absn1
        arr[p + 3][p + 6] = 3 * absn2 - 2 * degree * absn2
        arr[p + 3][p + 8] = 2 * absn1 * absn2
        arr[p + 3][p + 11] = - 2 * n1_sign * n2_sign * absn1 * absn2

    arr_inv = Matrix([[arr[j][i] for j in range(4 * (max_n - min_n + 1))] for i in range(4 * (max_n - min_n + 3))])
    return arr_inv


def init_rhs(min_n, max_n, bessels_rhs):
    rhs = [[0] for _ in range(4 * (max_n - min_n + 3))]
    begin_of_array = (-min_n + 1) * 4

    coeff00 = bessels_rhs.coeff(b0).coeff(p0)
    for j in range(min_n, max_n + 1):
        rhs[begin_of_array + 4 * j] = [coeff00.coeff(x ** j)]

    coeff01 = bessels_rhs.coeff(b0).coeff(p1)
    for j in range(min_n, max_n + 1):
        rhs[begin_of_array + 4 * j + 1] = [coeff01.coeff(x ** j)]

    coeff10 = bessels_rhs.coeff(b1).coeff(p0)
    for j in range(min_n, max_n + 1):
        rhs[begin_of_array + 4 * j + 2] = [coeff10.coeff(x ** j)]

    coeff11 = bessels_rhs.coeff(b1).coeff(p1)
    for j in range(min_n, max_n + 1):
        rhs[begin_of_array + 4 * j + 3] = [coeff11.coeff(x ** j)]

    return rhs


# Prints out the solution in Mathematica-friendly form
def print_nice_output(coeff, n1_sgn, n2_sgn, prfx, r, is_solution=True):
    str_bessel = ["" for _ in range(4)]
    str_bessel[0] = "eta00" + prfx + "[y_] :="
    str_bessel[1] = "eta01" + prfx + "[y_] :="
    str_bessel[2] = "eta10" + prfx + "[y_] :="
    str_bessel[3] = "eta11" + prfx + "[y_] :="

    for _ in range(0, 4):
        print(str_bessel[_])
        print(str(mcode(coeff[_])) + ";")

    sn1 = str(n1_sgn * n1)
    sn2 = str(n2_sgn * n2)

    math_output = "f[y_] := eta00" + prfx + "[y] BesselK[0, " + sn1 + " y] BesselK[0," + sn2 + "  y]"
    math_output += " + eta01" + prfx + "[y] BesselK[0, " + sn1 + " y] BesselK[1," + sn2 + "  y]"
    math_output += " + eta10" + prfx + "[y] BesselK[1, " + sn1 + " y] BesselK[0," + sn2 + "  y]"
    math_output += " + eta11" + prfx + "[y] BesselK[1, " + sn1 + " y] BesselK[1," + sn2 + "  y]"

    print(math_output)
    if is_solution:
        print("FullSimplify[y^2 D[D[f[y], y], y] - " + str(r * (r + 1)) + "f[y] - (n1 + n2)^2 y^2 f[y]]")


# Main function; finds a solution of the system of linear equations
def find_solutions(alpha, beta, r, n1_sgn, n2_sgn):
    # the eigenvalue
    lambd = r * (r + 1)

    # minN is the lowest degree of p_{i,j}(y)
    min_n = -4 * r - 4 - 4 * int(alpha) - 4 * int(beta)

    # maxN is the highest degree of p_{i,j}(y)
    max_n = 2

    # E_{\alpha} E_{\beta} has Bessel_{s1} and Bessel_{s2} in its right hand side
    s1 = int(alpha - 1 / 2)
    s2 = int(beta - 1 / 2)

    # This is the full expression of Bessel(s1, absn1 x) Bessel(s2, absn2 x)
    #t1 = datetime.datetime.now()
    # noinspection SpellCheckingInspection
    bessels_rhs = expand(bk1(s1, absn1 * x) * bk2(s2, absn2 * x))
    #t2 = datetime.datetime.now()
    #print(str(t2 - t1) + " needed for the Bessel functions recursion")
    t1 = datetime.datetime.now()

    arr_inv = init_arr(n1_sgn, n2_sgn, min_n, max_n, lambd)

    rhs = init_rhs(min_n, max_n, bessels_rhs)

    # all matrices shall be in the sympy form
    rhs_matrix = Matrix(rhs)
    lhs_matrix = Matrix(arr_inv)

    # solving a system of linear equations
    solution = linsolve((lhs_matrix, rhs_matrix))

    if not solution:  # this weird line actually checks if the set is empty!
        print("No solutions found")

    t2 = datetime.datetime.now()
    print(str(t2 - t1) + " seconds needed to find solutions for r = "+str(r)+", alpha = "+str(alpha)+", beta = "+str(beta))


    # noinspection SpellCheckingInspection
    coeff = [0 for _ in range(4)]

    for a_solution in solution:

        for i in range(4):

            for j in range(max_n - min_n + 1):
                p = 4 * j
                degree = j + min_n

                sol = a_solution[p + i].subs(absn1, n1_sgn * l1).subs(absn2, n2_sgn * l2)

                if sol != 0:
                    m, n = fraction(cancel(sol))
                    nom = simplify(m).subs(l1, n1).subs(l2, n2)
                    den = factor(n).subs(l1, n1).subs(l2, n2)
                    coeff[i] += (y ** degree) * nom / den

    return [coeff[_] for _ in range(0, 4)]


def calculateNminus3(n1_sgn, n2_sgn):
    coeffa3 = find_solutions(3 / 2, 3 / 2, 5, n1_sgn, n2_sgn)
    coeffa3 = [alpha3 * _ for _ in coeffa3]

    coeffa5 = find_solutions(3 / 2, 3 / 2, 5, n1_sgn, n2_sgn)
    coeffa7 = find_solutions(3 / 2, 3 / 2, 7, n1_sgn, n2_sgn)
    coeffa9 = find_solutions(3 / 2, 3 / 2, 9, n1_sgn, n2_sgn)

    coeffa5 = [alpha5 * _ for _ in coeffa5]
    coeffa7 = [alpha7 * _ for _ in coeffa7]
    coeffa9 = [alpha9 * _ for _ in coeffa9]

    coeffb5 = find_solutions(5 / 2, 5 / 2, 5, n1_sgn, n2_sgn)
    coeffb7 = find_solutions(5 / 2, 5 / 2, 7, n1_sgn, n2_sgn)
    coeffb9 = find_solutions(5 / 2, 5 / 2, 9, n1_sgn, n2_sgn)

    coeffb5 = [beta5 * _ for _ in coeffb5]
    coeffb7 = [beta7 * _ for _ in coeffb7]
    coeffb9 = [beta9 * _ for _ in coeffb9]

    coeffg5 = find_solutions(3 / 2, 7 / 2, 5, n1_sgn, n2_sgn)
    coeffg7 = find_solutions(3 / 2, 7 / 2, 7, n1_sgn, n2_sgn)
    coeffg9 = find_solutions(3 / 2, 7 / 2, 9, n1_sgn, n2_sgn)

    coeffg5 = [gamma5 * _ for _ in coeffg5]
    coeffg7 = [gamma7 * _ for _ in coeffg7]
    coeffg9 = [gamma9 * _ for _ in coeffg9]

    coeff = [
        coeffa3[_] + coeffa5[_] + coeffa7[_] + coeffa9[_] + coeffb5[_] + coeffb7[_] + coeffb9[_] + coeffg5[_] + coeffg7[
            _] + coeffg9[_] for _ in range(0, 4)]

    print_nice_output(coeff, n1_sgn, n2_sgn, "", 0, False)

def calculateNminus2(n1_sgn, n2_sgn):
    coeffd4 = find_solutions(5 / 2, 3 / 2, 4, n1_sgn, n2_sgn)
    coeffd4 = [d4 * _ for _ in coeffd4]

    coeffd6 = find_solutions(5 / 2, 3 / 2, 6, n1_sgn, n2_sgn)
    coeffd6 = [d6 * _ for _ in coeffd6]

    coeff = [
        coeffd4[_]+coeffd6[_] for _ in range(0, 4)]

    print_nice_output(coeff, n1_sgn, n2_sgn, "", 0, False)
