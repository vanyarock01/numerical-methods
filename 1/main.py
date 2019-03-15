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


if __name__ == '__main__':
    with open('matrix.json', 'r') as json_data:
        data = json.load(json_data)

        size = data[0]['size']
        matrix = data[0]['A']
        L, U = LU_decomposition(size, matrix)

        print("matrix")
        utils.matrix_print(matrix)
        print("L")
        utils.matrix_print(L)
        print("L")
        utils.matrix_print(U)
        print('L*U')
        a = np.array(L)
        b = np.array(U)
        r = np.matmul(a, b)
        utils.matrix_print(r)
