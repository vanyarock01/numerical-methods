import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


var = """Solve MNK system
i:      0      1      2      3       4      5
x:     -5     -3      1      1       3      5
y: 2.9442 2.8198 2.3562 0.7854 0.32175 0.1974
"""

raw_data = {
    "x": [-5, -3, 1, 1, 3, 5],
    "y": [2.9442, 2.8198, 2.3562, 0.7854, 0.32175, 0.1974],
    "x_": [0.0, 1.7, 3.4, 5.1, 6.8, 8.5],
    "y_": [0.0, 1.3038, 1.8439, 2.2583, 2.6077, 2.9155] 
}


def n_degree_polynom(x, y, n):
    N = len(x)
    max_pow = 0
    for _ in range(n):
        max_pow += 2

    xpow = [N] + [0 for _ in range(max_pow)]

    b = [0 for _ in range(n + 1)]

    for i in range(N):
        for k in range(1, len(xpow)):
            xpow[k] += x[i] ** (k)
        for k in range(len(b)):
            b[k] += y[i] * (x[i] ** k)

    shift = 0
    system = []
    for _ in range(n + 1):
        system.append([xpow[i + shift] for i in range(n + 1)])
        shift += 1

    a = np.linalg.solve(system, b)    
    
    def f(x):
        return a[0] + sum(a[i] * (x ** i) for i in range(1, len(a)))

    return a, f


    for i in range(n):
        xsum += x[i]
        ysum += y[i]
        xysum += x[i] * y[i]
        xpow += x[i] ** 2



def sum_square_error(x, y, f):
    return sum(
        (f(x[i]) - y[i]) ** 2 for i in range(len(x)))


def func_printer(x, f):
    n = len(x)
    def print_line(x):
        for i in range(n):
            print("{:7.4f}".format(x[i]), end=" | ")
        print()
    print_line(x)
    print_line([f(x[i]) for i in range(n)])


def show_plot(f, x, y, save_file, step=0.1):
    X = np.arange(x[0], x[-1], step)
    Y = []
    for i in range(len(f)):
        Y.append([f[i](val) for val in X])
    
    fig, ax = plt.subplots()
    for i in range(len(Y)):
        ax.plot(X, Y[i], label=f'degree={i}')

    ax.plot(x,  y, label='original')
    ax.legend(loc='upper right') 

    ax.grid()

    if save_file:
        fig.savefig(save_file)
        print(f"Saved in {save_file} succesfuly")
        plt.close(fig) 

    plt.show()


def main():
    x = raw_data["x"]
    y = raw_data["y"]
    f = []
    for degree in [1, 2, 3, 4, 5, 6]:
        print(f"degree {degree}")
        a, ft = n_degree_polynom(x, y, degree)
        f.append(ft)       
        func_printer(x, ft)
        print("sum square error: {:6.3f}\n"
            .format(sum_square_error(x, y, ft)))
    show_plot(f, x, y, "plot.png")


if __name__ == '__main__':
    main()
