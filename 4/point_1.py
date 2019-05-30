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
    return z * math.tan(x) - y * (math.cos(x) ** 2)


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
    x0, y0, z, h = p["xmin"], p["y0"], p["z0"], p["h"]
    n = getn(p["xmin"], p["xmax"], h)

    X, Y = [x0], [y0]
    for k in range(n):
        z += h * g(X[k], Y[k], z)

        X.append(X[k] + h)
        Y.append(Y[k] + h * f(X[k], Y[k], z))

    return X, Y


methods = {
    "euler": euler
}



def main():
    print(var)

    s = solve(precise, params)
    X, Y, labels = [s[0]], [s[1]], ["precision"]

    for label, bar in methods.items():
        x_, y_ = bar(g, f, params)
        X.append(x_)
        Y.append(y_)
        labels.append(label)

    draw_plot(X, Y, labels, "point_1.png")


if __name__ == '__main__':
    main()
