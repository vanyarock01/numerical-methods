#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import numpy as np
import utils
import algo


def show_LU(data):
    size = data[0]['size']
    A = data[0]['A']
    B = data[0]['B']

    L, U = algo.LU_decomposition(size, A)

    print("matrix")
    utils.matrix_print(A)
    print("L")
    utils.matrix_print(L)
    print("U")
    utils.matrix_print(U)
    print('L * U')

    X = algo.LU_solve(size, L, U, B)
    a = np.array(L)
    b = np.array(U)
    r = np.matmul(a, b)
    utils.matrix_print(r)
    print('X')
    utils.vector_print(X)
    print('NP')
    utils.vector_print(np.linalg.solve(A, B))


def show_tridiagonal(data):
    size = data[1]['size']
    A = data[1]['A']
    B = data[1]['B']
    x = algo.tridiagonal(size, A, B)
    utils.vector_print(x)
    utils.vector_print(np.linalg.solve(A, B))

if __name__ == '__main__':
    with open('matrix.json', 'r') as json_data:
        data = json.load(json_data)
        show_LU(data)
        #show_tridiagonal(data)
