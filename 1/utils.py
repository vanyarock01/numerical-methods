def matrix_print(matrix):
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            print('{:6.2f}'.format(float(matrix[i][j])), end=' ')
        print()
    print()
