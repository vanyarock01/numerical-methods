#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import numpy as np
import utils


def LU_decomposition(size, matrix):
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
            U[j][i] = matrix[j][i] - s

        for j in range(i, size):
            s = sum(U[k][i] * L[j][k] for k in range(i))
            L[j][i] = (matrix[j][i] - s) / U[i][i]

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


if __name__ == '__main__':
    with open('matrix.json', 'r') as json_data:
        data = json.load(json_data)

        size = data[0]['size']
        A = data[0]['A']
        B = data[0]['B']

        L, U = LU_decomposition(size, A)

        print("matrix")
        utils.matrix_print(A)
        print("L")
        utils.matrix_print(L)
        print("L")
        utils.matrix_print(U)
        print('L*U')

        X = LU_solve(size, L, U, B)
        a = np.array(L)
        b = np.array(U)
        r = np.matmul(a, b)
        utils.matrix_print(r)
        utils.vector_print(X)
