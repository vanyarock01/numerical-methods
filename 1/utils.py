def matrix_print(matrix, header=None):
    if header:
        print(header)
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            print('{:6.2f}'.format(float(matrix[i][j])), end=' ')
        print()


def vector_print(vec, header=None):
    if header:
        print(header)
    for i in range(len(vec)):
        print('{:6.2f}'.format(float(vec[i])), end=' ')
    print()


def mv_mult(size, matrix, vector):
    R = [0.0 for _ in range(size)]
    for i in range(size):
        for j in range(size):
            R[i] += matrix[i][j] * vector[j]
    return R


def mm_mult(size, matrix1, matrix2):
    R = [[0.0 for _ in range(size)] for _ in range(size)]

    for i in range(size):
        for j in range(size):
            for k in range(size):
                R[i][j] += matrix1[i][k] * matrix2[k][j]
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
