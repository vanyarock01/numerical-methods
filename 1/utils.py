def matrix_print(matrix):
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            print('{:6.2f}'.format(float(matrix[i][j])), end=' ')
        print()
    print()


def vector_print(vec):
    for i in range(len(vec)):
        print('{:6.2f}'.format(float(vec[i])), end=' ')
    print()
