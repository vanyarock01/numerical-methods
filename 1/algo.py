def LU_decomposition(size, A):
    L = [[0.0 for _ in range(size)] for _ in range(size)]
    U = [[0.0 for _ in range(size)] for _ in range(size)]

    for i in range(size):
        L[i][i] = 1
        max_elem = abs(U[i][i])
        row = i
        for k in range(i + 1, size):
            if abs(U[k][i]) > max_elem:
                row = k
                max_elem = abs(U[k][i])

        for k in range(i, size):
            U[row][k], U[i][k] = U[i][k], U[row][k]

        for j in range(i + 1):
            s = sum(U[k][i] * L[j][k] for k in range(j))
            U[j][i] = A[j][i] - s

        for j in range(i, size):
            s = sum(U[k][i] * L[j][k] for k in range(i))
            L[j][i] = (A[j][i] - s) / U[i][i]

    return L, U


def LU_solve(size, L, U, B):
    Z = [0 for _ in range(size)]
    for i in range(0, size):
        Z[i] = B[i] - sum(L[i][j] * Z[j] for j in range(i - 1))

    X = [0 for _ in range(size)]
    for i in reversed(range(size)):
        X[i] = (1 / U[i][i]) * (Z[i] - sum(U[i][j] * X[j]
                                           for j in range(i + 1, size)))
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
