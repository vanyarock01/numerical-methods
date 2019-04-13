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
    X = algo.LU_solve(size, L, U, B)

    npL = np.array(L)
    npU = np.array(U)
    np_matrix = np.matmul(npL, npU)

    utils.matrix_print(A, header="matrix")
    utils.matrix_print(L, header="L")
    utils.matrix_print(U, header="U")

    utils.matrix_print(np_matrix, header="L * U")
    utils.vector_print(X, header="X")
    utils.vector_print(np.linalg.solve(A, B), header="LINALG SOLVE")


def show_tridiagonal(data):
    size = data[1]['size']
    A = data[1]['A']
    B = data[1]['B']
    utils.matrix_print(A, header="A")
    utils.vector_print(B, header="B")
    X = algo.tridiagonal(size, A, B)
    utils.vector_print(X, header="X")
    utils.vector_print(np.linalg.solve(A, B), header="LINALG SOLVE")


def show_simple_iteration(data):
    size = data[2]['size']
    A = data[2]['A']
    B = data[2]['B']
    utils.matrix_print(A, header="A")
    utils.vector_print(B, header="B")
    X = algo.simple_iteration(size, A, B)
    utils.vector_print(X, header="X")
    utils.vector_print(np.linalg.solve(A, B), header="LINALG SOLVE")

def show_zeidel_method(data):
    size = data[2]['size']
    A = data[2]['A']
    B = data[2]['B']
    utils.matrix_print(A, header="A")
    utils.vector_print(B, header="B")
    X = algo.zeidel_method(size, A, B)
    utils.vector_print(X, header="X")
    utils.vector_print(np.linalg.solve(A, B), header="LINALG SOLVE")

if __name__ == '__main__':
    with open('matrix.json', 'r') as json_data:
        data = json.load(json_data)
        print("zeidel")
        show_zeidel_method(data)