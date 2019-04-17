import math

"""
System:
{
	x1 - cos(x2) = 2,
	x2 - sin(x1) = 2
}
"""

# first function ---


def f1(x1, x2):
    return x1 - math.cos(x2) - 2


def df1_dx1(x1, x2):
    return 1


def df1_dx2(x1, x2):
    return math.sin(x2)


# second function ---

def f2(x1, x2):
    return x2 - math.sin(x1) - 2


def df2_dx1(x1, x2):
    return -math.cos(x1)


def df2_dx2(x1, x2):
    return 1


func = {
    "f1": f1,
    "f2": f2,

    "df1_dx1": df1_dx1,
    "df1_dx2": df1_dx2,
    "df2_dx1": df2_dx1,
    "df2_dx2": df2_dx2
}

# shape = 2 only


def det(A):
    return A[0][0] * A[1][1] - A[0][1] * A[1][0]


def getx(x1, x2, fc):
    A1 = [
        [fc["f1"](x1, x2), fc["df1_dx2"](x1, x2)],
        [fc["f2"](x1, x2), fc["df2_dx2"](x1, x2)]
    ]
    A2 = [
        [fc["df1_dx1"](x1, x2), fc["f1"](x1, x2)],
        [fc["df2_dx1"](x1, x2), fc["f2"](x1, x2)]
    ]

    J = [
        [fc["df1_dx1"](x1, x2), fc["df1_dx2"](x1, x2)],
        [fc["df2_dx1"](x1, x2), fc["df2_dx2"](x1, x2)]
    ]
    detJ = det(J)
    return x1 - det(A1) / detJ, x2 - det(A2) / detJ


def norm(x, x_prev):
    return max(x[0] - x_prev[0], x[1] - x_prev[1])


def newton(x1, x2, precision=0.01):
    x1_prev, x2_prev = x1_cur, x2_cur = x1, x2
    k = 0
    while True:
        k += 1
        x1_cur, x2_cur = getx(x1_prev, x2_prev, func)
        
        epsilon = norm(
            (x1_cur,  x2_cur), (x1_prev, x2_prev))
        print(
            "x1: {:7.4f} | x2: {:7.4f} | eps: {:7.4f}"
            .format(x1_cur, x2_cur, epsilon))
        
        if epsilon <= precision:
            break
        x1_prev, x2_prev = x1_cur, x2_cur
    return x1_cur, x2_cur


def simple_iteration():
    pass


if __name__ == '__main__':
    newton(x1=0.5, x2=2, precision=0.001)
