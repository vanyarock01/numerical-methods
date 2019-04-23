import utils
import math
import numpy as np
from numpy.linalg import norm, solve, inv


def swap_lines(matrix, i, k):
    matrix[i], matrix[k] = matrix[k].copy(), matrix[i].copy()


def permute(A):
    n = len(A)
    p = []
    m = [row.copy() for row in A]
    for j in range(n - 1):
        tmp = [(m[i][j], i) for i in range(j, n)]
        idx = max(tmp, key=lambda x: abs(x[0]))[1]
        if idx != j:
            p.append((j, idx))
            swap_lines(m, j, idx)
    return p, m


def create_p_matrix(matrix, p):
    p_m = matrix.copy()

    for i, j in p:
        p_m[i], p_m[j] = p_m[j], p_m[i]
    return p_m


def LU_decomposition(size, A):
    L = [[0.0 for _ in range(size)] for _ in range(size)]

    p, U = permute(A)

    for i in range(size):
        L[i][i] = 1
        for j in range(i + 1, size):
            L[j][i] = U[j][i] / U[i][i]
            for k in range(i + 1, size):
                U[j][k] -= L[j][i] * U[i][k]

    for i in range(1, size):
        for j in range(i):
            U[i][j] = 0.0

    return L, U, p


def LU_solve(size, L, U, B, p):
    Z = [0 for _ in range(size)]
    X = [0 for _ in range(size)]
    Bp = create_p_matrix(B, p)
    for i in range(size):
        Z[i] = Bp[i] - sum(L[i][j] * Z[j] for j in range(i))

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


def get_equal_form(A, B):
    size = len(A)
    alpha = [[0.0 for _ in range(size)] for _ in range(size)]
    beta = [0.0 for _ in range(size)]

    for i in range(size):
        beta[i] = float(B[i]) / A[i][i]
        for j in range(size):
            alpha[i][j] = -A[i][j] / A[i][i]
        alpha[i][i] = 0
    return alpha, beta


def simple_iteration(size, A, B, precision=0.01):
    A = [row.copy() for row in A]
    alpha, beta = get_equal_form(A, B)

    q = utils.m_norm(size, alpha)
    cur = beta.copy()
    last = []
    while True:
        last = cur
        mult = utils.mv_mult(size, alpha, cur)
        cur = utils.vv_add(size, beta, mult)
        norma = utils.v_norm(
            size, utils.vv_substr(size, cur, last))
        if norma * q / (1 - q) <= precision:
            break

    return cur


def get_norm(A, alpha, S, C):
    return norm(C, np.inf) / (1. - norm(S, np.inf))


def zeidel_method(size, A, B, precision=0.01):
    alpha, beta = get_equal_form(A, B)

    alpha = np.array(alpha)
    beta = np.array(beta)

    b = np.tril(alpha, -1)

    K = alpha - b

    T1 = inv(
        np.eye(size, size) - b) @ K

    T2 = inv(
        np.eye(size, size) - b) @ beta

    X = T2
    C = get_norm(A, alpha, T1, K)

    cov = True
    while cov:
        X_next = T2 + T1 @  X
        if C * norm(X_next - X, np.inf) <= precision:
            cov = False
        X = X_next

    return X_next


def rotate_jacobi(size, A, precision=0.01):
    Ak = [row.copy() for row in A]

    idx = range(size)
    U = [[0. if i != j else 1. for i in idx] for j in idx]

    cov = False
    while not cov:
        ik, jk = 0, 1
        for i in range(size - 1):
            for j in range(i + 1, size):
                if abs(Ak[i][j]) > abs(Ak[ik][jk]):
                    ik, jk = i, j

        if Ak[ik][ik] == Ak[jk][jk]:
            phi = math.pi / 4
        else:
            phi = 0.5 * math.atan(
                2 * Ak[ik][jk] / (Ak[ik][ik] - Ak[jk][jk]))

        Uk = [[0. if i != j else 1. for i in idx] for j in idx]
        Uk[ik][jk] = math.sin(phi)
        Uk[jk][ik] = -Uk[ik][jk]

        Uk[ik][ik] = math.cos(phi)
        Uk[jk][jk] = math.cos(phi)

        tmp = utils.mm_mult(Uk, Ak)

        Uk[ik][jk], Uk[jk][ik] = Uk[jk][ik], Uk[ik][jk]

        Ak = utils.mm_mult(tmp, Uk)
        U = utils.mm_mult(U, Uk)

        accum = 0
        for i in range(size - 1):
            for j in range(i + 1, size):
                accum += Ak[i][j] ** 2

        avg = math.sqrt(accum)
        if avg < precision:
            cov = True

    return [Ak[i][i] for i in range(size)], U


def sign(x):
    if x > 0:
        return 1
    elif x < 0:
        return -1
    else:
        return 0


def QR_decomposition(size, A, precision=0.01):
    Ak = [row.copy() for row in A]

    E = utils.e(size)
    Q = utils.e(size)

    for i in range(size):
        V = [0 for _ in range(size)]

        V[i] = Ak[i][i] + sign(Ak[i][i]) * math.sqrt(
            sum(Ak[j][i] ** 2 for j in range(i, size)))

        for k in range(i + 1, size):
            V[k] = Ak[k][i]

        Vt = [V]
        V = [[V[i]] for i in range(size)]

        M = utils.mm_mult(V, Vt)
        C = utils.mm_mult(Vt, V)[0][0]
        for j in range(size):
            for k in range(size):
                M[j][k] /= C
                M[j][k] *= 2

        Hk = utils.mm_substr(E, M)
        Q = utils.mm_mult(Q, Hk)
        Ak = utils.mm_mult(Hk, Ak)

    return Q, Ak


def roots(a, types):
    n = len(a)
    A = np.array(a)
    solution = []
    k = 0
    for t in types:
        if t == 'real':
            solution.append(A[k, k])
        else:

            A11 = A[k, k]
            A12 = A21 = A22 = 0

            if k + 1 < n:
                A12 = A[k, k + 1]
                A21 = A[k + 1, k]
                A22 = A[k + 1, k + 1]

            solution.extend(np.roots(
                (1, -A11 - A22, A11 * A22 - A12 * A21)))
            k += 1
        k += 1
    return solution





def check(matrix, precision=0.01):
    n = len(matrix)
    check = []
    k = 0

    def square_norm(x):
        return sum(e ** 2 for e in x) ** 0.5

    def get_column(a, k):
        return [a[i][k] for i in range(k+1, len(a))]

    while k < n:
        if square_norm(get_column(matrix, k)) <= precision:
            check.append('real')
        elif square_norm(get_column(matrix, k + 1)) <= precision:
            check.append('img')
            k += 1
        else:
            check.append(None)
        k += 1
    return check


def QR_method(size, A, precision=0.01):
    Ak = [row.copy() for row in A]
    it, max_it = 0, 100

    step = True

    while it < max_it:
        Q, R = QR_decomposition(size, Ak)
        Ak = utils.mm_mult(R, Q)
        utils.matrix_print(Q, header=f"A{it}")

        types = check(Ak, precision)
        if all(types):
            if step:
                step = False
            else:
                return roots(Ak, types)
        it += 1

    return None
