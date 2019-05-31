import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def draw_plot(X, Y, labels, save_file="plot.png"):
    fig, ax = plt.subplots()
    for i in range(len(X)):
        ax.plot(X[i], Y[i], label=f"{labels[i]}")

    ax.legend(loc='upper right')
    ax.grid()
    if save_file:
        fig.savefig(save_file)
        print(f"Saved in {save_file} succesfuly")
        plt.close(fig)


params = {
    "h":    0.1,
    "xmin": 0.0,
    "xmax": 1.0,
    "y0":   2.0,
    "z0":   0.0
}


var = """+++++++++++++++++++++++++++++++
KOSHI
+++++++++++++++++++++++++++++++
+ y'' + y'tg(x) - ycos^2x = 0 
+ y(0)  = 2                   
+ y'(0) = 0                
+ x in [0, 1], h = 0.1
+++++++++++++++++++++++++++++++
+ SOLVE
+++++++++++++++++++++++++++++++
+ y = e^(sin(x)) + e^(-sin(x))
+++++++++++++++++++++++++++++++
"""


def precise(x):
    return math.e ** math.sin(x) + math.e ** (-math.sin(x))


def g(x, y, z):
    return -z * math.tan(x) + y * (math.cos(x) ** 2)


def f(x, y, z):
    return z


def getn(a, b, step=0.1):
    return int((b - a) / step)


def solve(f, p):
    X = [x for x in
         np.arange(p["xmin"], p["xmax"] + p["h"], p["h"])]
    Y = [f(x) for x in X]
    return X, Y


def euler(g, f, p):
    x0, y0, z0, h = p["xmin"], p["y0"], p["z0"], p["h"]
    n = getn(p["xmin"], p["xmax"], h)

    X, Y, Z = [x0], [y0], [z0]
    for k in range(n):
        Y.append(Y[k] + h * f(X[k], Y[k], Z[k]))
        Z.append(Z[k] + h * g(X[k], Y[k], Z[k]))
        X.append(X[k] + h)

    return X, Y, Z


def runge_kutta(g, f, p):
    degree = 4
    a = [0, 1/2, 1/2,   1]
    c = [1/6, 1/3, 1/3, 1/6]
    b = [[],
         [1/2, ],
         [0,  1/2],
         [0,    0, 1/2]]

    err = []
    x0, y0, z0, h = p["xmin"], p["y0"], p["z0"], p["h"]
    n = getn(p["xmin"], p["xmax"], h)
    X, Y, Z = [x0], [y0], [z0]
    K, L = [0.] * degree, [0.] * degree
    for k in range(n - 1):
        y_, z_ = 0, 0
        for i in range(degree):
            sum_K = sum( b[i][j] * K[j] for j in range(i) )
            sum_L = sum( b[i][j] * L[j] for j in range(i) )

            K[i] = h * \
                f(X[k] + a[i] * h, Y[k] + sum_K, Z[k] + sum_L)
            L[i] = h * \
                g(X[k] + a[i] * h, Y[k] + sum_K, Z[k] + sum_L)
            
            y_ += c[i] * K[i]
            z_ += c[i] * L[i]
        Y.append(Y[k] + y_)
        Z.append(Z[k] + z_)
        X.append(X[k] + h)
        if abs(K[0] - K[1]) > 0:
            err.append((K[1] - K[2]) / (K[0] - K[1]))
        else:
            err.append(0.0)
        # print(err)
    return X, Y, Z


def adams(g, f, p):
    x0, y0, z0, h = p["xmin"], p["y0"], p["z0"], p["h"]
    n = getn(p["xmin"], p["xmax"], h)

    X, Y, Z = runge_kutta(g, f, p)

    for k in range(3, n):
        y_ = 55 * f(X[k], Y[k], Z[k]) \
            - 59 * f(X[k - 1], Y[k - 1], Z[k - 1]) \
            + 37 * f(X[k - 2], Y[k - 2], Z[k - 2]) \
            - 9 * f(X[k - 3], Y[k - 3], Z[k - 3])
        z_ = 55 * g(X[k], Y[k], Z[k]) \
            - 59 * g(X[k - 1], Y[k - 1], Z[k - 1]) \
            + 37 * g(X[k - 2], Y[k - 2], Z[k - 2]) \
            - 9 * g(X[k - 3], Y[k - 3], Z[k - 3])
        Y.append(Y[k] + h / 24 * y_)
        Z.append(Z[k] + h / 24 * z_)
        X.append(X[k] + h)
    return X, Y, Z


methods = {
    "euler":       euler,
    "runge_kutta": runge_kutta,
    "adams":       adams
}


def main():
    print(var)

    s = solve(precise, params)
    X, Y, labels = [s[0]], [s[1]], ["precision"]

    for label, bar in methods.items():
        x_, y_, _ = bar(g, f, params)
        X.append(x_)
        Y.append(y_)
        labels.append(label)

    draw_plot(X, Y, labels, "point_1.png")


if __name__ == '__main__':
    main()
