import math

var = """y = arcctg(x)
a) Xi = -3, -1, 1, 3
b) Xi = -3,  0, 1, 3

X* = -0.5
"""


def polynom_printer(P):
    print("P(x) = ", end='')
    plus = False
    for member in P:
        if plus:
            print(" + ", end='')
        else:
            plus = True

        print("{:5.3f}".format(member[0]), end='')
        for i in range(1, len(member)):
            print("(x - {:3.1f})".format(member[i]), end='')
    print()


def calc(P, x):
    f = 0
    for member in P:
        f += member[0] * mult(
            x - member[i] for i in range(1, len(member)))
    return f


def mult(gen):
    r = 1
    for e in gen:
        r *= e
    return r


def func(x):
    return math.atan(x)


def getl(i, x):
    r = [1.0]
    for j in range(len(x)):
        if j != i:
            r.append(x[j])
            r[0] /= x[i] - x[j]
    return r


def lagrange(x, f):
    """ return
        L : [ [0.86, -3, 5], [0.42, 4, -2] ]
        ~ : 0.86(x - 3)(x + 5) + 0.42(x + 4)(x - 2)
    """
    n = len(x)
    y = [f(e) for e in x]

    L = []
    for i in range(n):
        li = getl(i, x)
        li[0] *= y[i]
        L.append(li)

    return L


def newton(x, f):
    n = len(x)

    D = [[f(e) for e in x]]

    N = [ [D[0][0]] ]  # f(x0)
    K = 1
    for i in range(n - 1):
        D.append([])
        for j in range(len(D[i]) - 1):
            D[i + 1].append(
                (D[i][j] - D[i][j + 1]) / (x[j] - x[j + K]))
        N.append(
            [D[i + 1][0]] + [x[k] for k in range(i + 1)])
        K += 1
    """
    debug output calculate heap

    for i in range(n):
        print("{:6.3f}".format(x[i]), end=' ')
    print('\n')
    for i in range(n):
        for j in range(len(D[i])):
            print("{:6.3f}".format(D[i][j]), end=' ')
        print()
    """
    return N


def test(x):
    return math.sin(math.pi*x/6)


def main():
    x_star = -0.5

    x_a = [-3, -1, 1, 3]
    x_b = [-3,  0, 1, 3]
    x = 1.5
    x_test = [0, 1, 2, 3]
    print(var)

    print("LAGRANGE\n")
    L = lagrange(x_b, func)
    polynom_printer(L)

    print()
    Lx = calc(L, x_star)
    print("P({}) = {:6.3f}".format(x_star, Lx))

    print("NEWTON\n")
    N = newton(x_a, func)
    polynom_printer(N)
    print()
    Nx = calc(N, x_star)
    print("P({}) = {:6.3f}".format(x_star, Nx))

    print("FUNCTION\n")
    Fx = func(x_star)
    print("F({:6.3f}) = {:6.3f}".format(x_star, Fx))
    print("|F(x) - L(x)| == {:7.4}".format(abs(Lx - Fx)))
    print("|F(x) - N(x)| == {:7.4}".format(abs(Nx - Fx)))

if __name__ == '__main__':
    main()
