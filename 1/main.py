#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import numpy as np
import utils
import algo
from pprint import pprint


def show_LU(data, idx):
    print("LAB 1.1: LU decomposition", '\n')
    size = data[idx]['size']
    A = data[idx]['A']
    B = data[idx]['B']

    L, U, p = algo.LU_decomposition(size, A)
    X = algo.LU_solve(size, L, U, B, p)

    npL = np.array(L)
    npU = np.array(U)
    np_matrix = np.matmul(npL, npU)

    utils.matrix_print(A, header="matrix")
    utils.matrix_print(L, header="L")
    utils.matrix_print(U, header="U")

    utils.matrix_print(np_matrix, header="L * U")
    utils.vector_print(X, header="X")
    utils.vector_print(np.linalg.solve(A, B), header="LINALG SOLVE")


def show_tridiagonal(data, idx):
    print("LAB 1.2: tridiagonal", '\n')
    size = data[idx]['size']
    A = data[idx]['A']
    B = data[idx]['B']
    utils.matrix_print(A, header="A")
    utils.vector_print(B, header="B")
    X = algo.tridiagonal(size, A, B)
    utils.vector_print(X, header="X")
    utils.vector_print(np.linalg.solve(A, B), header="LINALG SOLVE")


def show_simple_iteration(data, idx):
    print("LAB 1.3(a): Simple iteration", '\n')

    size = data[idx]['size']
    A = data[idx]['A']
    B = data[idx]['B']
    epsilon = data[idx]['precision']
    utils.matrix_print(A, header="A")
    utils.vector_print(B, header="B")
    X = algo.simple_iteration(size, A, B, epsilon)
    utils.vector_print(X, header="X")
    utils.vector_print(np.linalg.solve(A, B), header="LINALG SOLVE")


#        seidel
def show_zeidel_method(data, idx):
    print("LAB 1.3(b): Seidel method", '\n')

    size = data[idx]['size']
    A = data[idx]['A']
    B = data[idx]['B']
    epsilon = data[idx]['precision']
    utils.matrix_print(A, header="A")
    utils.vector_print(B, header="B")
    X = algo.zeidel_method(size, A, B, epsilon)
    utils.vector_print(X, header="X")
    utils.vector_print(np.linalg.solve(A, B), header="LINALG SOLVE")


def show_iteration(data, idx):
    show_simple_iteration(data, idx)
    show_zeidel_method(data, idx)


def show_jacobi_rotate(data, idx):
    print("LAB 1.4: Jacobi rotate", '\n')
    size = data[idx]['size']
    A = data[idx]['A']
    epsilon = data[idx]['precision']
    utils.matrix_print(A, header="A")
    X, U = algo.rotate_jacobi(size, A, epsilon)

    utils.vector_print(X, header="X")
    utils.matrix_print(U, header="U")

    X, U = np.linalg.eig(A)
    utils.vector_print(X, header="LINALG SOLVE: X")
    utils.matrix_print(U, header="LINALG SOLVE: U")


def show_QR_method(data, idx):
    print("LAB 1.5: QR decompose", '\n')
    size = data[idx]['size']
    A = data[idx]['A']
    epsilon = data[idx]['precision']
    utils.matrix_print(A, header="A")
    X = algo.QR_method(size, A, epsilon)
    pprint(X)
    X, U = np.linalg.eig(A)
    pprint(X)


task = [
    show_LU,
    show_tridiagonal,
    show_iteration,
    show_jacobi_rotate,
    show_QR_method
]

if __name__ == '__main__':
    with open('matrix.json', 'r') as json_data:
        data = json.load(json_data)
        task[-1](data, 5)
        """
        for idx, method in zip(range(len(task)), task):
            print('\n', '-'*20)
            method(data, idx)
        """