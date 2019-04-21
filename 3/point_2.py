import math


var = """X* = -0.5
i :      0       1       2        3        4
Xi: -    3      -1       1        3        5
Fi: 2.8198, 2.3562, 0.7854, 0.32175, 0.19740.
"""

raw_data = {
    "X*": -0.5,
    "n": 4,
    "x": [-3.0, -1.0, 1.0, 3.0, 5.0],
    "f": [2.8198, 2.3562, 0.7854, 0.32175, 0.19740],

    "x_": [0.0, 1.0, 2.0, 3.0, 4.0],
    "f_": [0.0, 1.8415, 2.9093, 3.1411, 3.2432]
}


def spline_printer(S, x):
    spline_template = "S(x)={:5.2f} + {:5.2f}(x - {:5.2f}) + {:5.2f}(x - {:5.2f})^2 + {:5.2f}(x - {:5.2f})^3"
    for i in range(1, len(x)):
        print("{} | [{:2.0f}, {:2.0f}]"
              .format(i, x[i-1], x[i]), end=' | ')
        print(spline_template
              .format(
                  S[i-1][1],
                  S[i-1][2],
                  S[i-1][0],
                  S[i-1][3],
                  S[i-1][0],
                  S[i-1][4],
                  S[i-1][0]))

    # simple S line example: [x[i-1], A[i], B[i], C[i], D[i]]


def tridiagonal(a, b, c, d):
    n = len(d)

    a = a.copy()
    b = b.copy()
    c = c.copy()
    d = d.copy()
    for i in range(1, n):
        m = a[i - 1] / b[i - 1]
        b[i] = b[i] - m * c[i - 1]
        d[i] = d[i] - m * d[i - 1]

    x = b.copy()
    x[n - 1] = d[n - 1] / b[n - 1]

    for i in reversed(range(0, n - 1)):
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i]
    return x


def cubic_spline(x, f):
    n = len(x)

    def h(i): return x[i] - x[i-1]

    a, b, c, d = [], [], [], []

    a = [h(i - 1) for i in range(3, n)]
    b = [2 * (h(i-1) + h(i)) for i in range(2, n)]
    c = [h(i) for i in range(2, n-1)]
    d = [3 * ((f[i] - f[i-1]) / h(i) - ((f[i-1] - f[i-2]) / h(i-1)))
         for i in range(2, n)]

    C = [0, 0] + tridiagonal(a, b, c, d)
    A = [0] + [f[i] for i in range(n-1)]
    B = [0]
    for i in range(1, n-1):
        B.append((f[i] - f[i - 1]) / h(i) - 1/3 * h(i) * (C[i + 1] + 2 * C[i]))
    B.append((f[n-1] - f[n-2]) / h(n-1) - 2/3 * h(n-1) * C[n-1])
    D = [0] + [(C[i+1] - C[i]) / (3 * h(i)) for i in range(1, n-1)]
    D.append(-C[n-1] / 3 * h(n-1))

    S = []
    for i in range(1, n):
        S.append(
            [x[i-1], A[i], B[i], C[i], D[i]])
    return S


def main():
    print(var)

    S = cubic_spline(raw_data["x"], raw_data["f"])
    spline_printer(S, raw_data["x"])


if __name__ == '__main__':
    main()
