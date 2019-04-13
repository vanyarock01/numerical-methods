import utils
import math

# TODO: перестановки
def LU_decomposition(size, A):
    L = [[0.0 for _ in range(size)] for _ in range(size)]
    U = [row.copy() for row in A]

    for i in range(size):
        L[i][i] = 1
        if U[i][i] == 0:
            max_elem = abs(U[i][i])
            row = i
            for k in range(i + 1, size):
                if abs(U[k][i]) > max_elem:
                    row = k
                    max_elem = abs(U[k][i])

            for k in range(i, size):
                U[k][row], U[k][i] = U[k][i], U[k][row]

        for j in range(i + 1, size):
            L[j][i] = U[j][i] / U[i][i]
            for k in range(i + 1, size):
                U[j][k] -= L[j][i] * U[i][k]

    for i in range(1, size):
        for j in range(i):
            U[i][j] = 0.0

    return L, U


def LU_solve(size, L, U, B):
    Z = [0 for _ in range(size)]
    X = [0 for _ in range(size)]

    for i in range(size):
        Z[i] = B[i] - sum(L[i][j] * Z[j] for j in range(i))

    for i in reversed(range(size)):
        z = sum(U[i][j] * X[j] for j in range(i + 1, size))
        X[i] = (Z[i] - z) / U[i][i]
    return X


def tridiagonal(size, A, B):
    a = [A[i + 1][i] for i in range(size - 1)]
    b = [A[i][i] for i in range(size)]
    c = [A[i][i + 1] for i in range(size - 1)]
    d = B.copy()

    for i in range(1, size):
        m = a[i - 1] / b[i - 1]
        b[i] = b[i] - m * c[i - 1]
        d[i] = d[i] - m * d[i - 1]

    x = b.copy()
    x[size - 1] = d[size - 1] / b[size - 1]

    for i in reversed(range(0, size - 1)):
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i]

    return x


def simple_iteration(size, A, B, precision=0.01):
    A = [row.copy() for row in A]

    alpha = [[0.0 for _ in range(size)] for _ in range(size)]
    beta = [0.0 for _ in range(size)]

    for i in range(size):
        beta[i] = float(B[i]) / A[i][i]
        for j in range(size):
            alpha[i][j] = -A[i][j] / A[i][i]
        alpha[i][i] = 0

    q = utils.m_norm(size, alpha)
    cur = beta.copy()
    last = []
    while True:
        last = cur
        mult = utils.mv_mult(size, alpha, cur)
        cur = utils.vv_add(size, beta, mult)
        norma = utils.v_norm(
            size, utils.vv_substr(size, cur, last))
        if norma <= precision * (1 - q) / q:
            break

    return cur


def zeidel_method(size, A, B, precision=0.01):
    x = [0.0 for i in range(size)]

    cov = False
    while not cov:
        x_next = x.copy()
        for i in range(size):
            S1 = sum(A[i][j] * x_next[j] for j in range(i))
            S2 = sum(A[i][j] * x[j] for j in range(i + 1, size))
            x_next[i] = (B[i] - S1 - S2) / A[i][i]

        cov = math.sqrt(sum((x_next[i] - x[i]) ** 2 for i in range(size))) <= precision
        x = x_next

    return x

