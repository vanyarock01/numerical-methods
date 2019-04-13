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


def mv_mult(size, matrix, vector):
    R = [0.0 for _ in range(size)]
    for i in range(size):
        for j in range(size):
            R[i] += matrix[i][j] * vector[j]
    return R


def vv_substr(size, vec1, vec2):
    return [vec1[i] - vec2[i] for i in range(size)]


def vv_add(size, vec1, vec2):
    return [vec1[i] + vec2[i] for i in range(size)]


def v_norm(size, vec):
    return sum(abs(vec[i]) for i in range(size))


def m_norm(size, matrix):
    v = [sum(matrix[i][j] for i in range(size)) for j in range(size)]
    return max(v)


def simple_iteration(size, A, B, precision=0.01):
    A = [row.copy() for row in A]
    # X = [0.0 for _ in range(size)]

    alpha = [[0.0 for _ in range(size)] for _ in range(size)]
    beta = [0.0 for _ in range(size)]

    for i in range(size):
        beta[i] = float(B[i]) / A[i][i]
        for j in range(size):
            alpha[i][j] = -A[i][j] / A[i][i]
        alpha[i][i] = 0

    q = m_norm(size, alpha)
    print(q)
    cur = beta.copy()
    last = []
    while True:
        last = cur
        mult = mv_mult(size, alpha, cur)
        cur = vv_add(size, beta, mult)
        print(v_norm(size, vv_substr(size, cur, last)))
        if v_norm(size, vv_substr(size, cur, last)) <= precision * (1 - q) / q:
            break

    return cur
