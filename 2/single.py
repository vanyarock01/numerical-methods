import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def func(x):
    return math.exp(x) - x ** 3 + 3 * (x ** 2) - 2 * x - 3


def dfunc(x):
    return math.exp(x) - 3 * (x ** 2) + 6 * x - 2


def ddfunc(x):
    return math.exp(x) - 6 * x + 6


def phi(x):
    return math.log(x ** 3 - 3 * (x ** 2) + 2 * x + 3)


def dphi(x):
    return (3 * (x ** 2) - 6 * x + 2) / \
        (x ** 3 - 3 * (x ** 2) + 2 * x + 3)


def newton(f, df, x0, precision=0.001):
    x = x_next = x0
    k = 0
    cov = True
    while cov:
        k += 1
        x_next = x - f(x) / df(x)
        print(f"x: {x_next} k: {k} |x_next - x|: {abs(x_next - x)}")
        if abs(x_next - x) <= precision:
            cov = False
        x = x_next


def simple_iteration(phi, dphi, a, b, precision=0.01):
    q = max(abs(dphi(a)), abs(dphi(b)))
    x = (b - a) / 2
    x_next = x
    k = 0
    cov = True
    while cov:
        k += 1
        x_next = phi(x)

        print(f"x: {x_next} k: {k} |x_next - x|: {abs(x_next - x)}")
        if abs(x_next - x) <= precision:
            cov = False
        x = x_next
        if k == 10:
            break


def show_plot(f, df,  x, save_file, step=0.5,  ddf=None):
    X = np.arange(x[0], x[-1], step)
    Y = [f(val) for val in X]
    dY = [df(val) for val in X]
    
    if ddf:
        ddY = [ddf(val) for val in X]
    fig, ax = plt.subplots()
    ax.plot(X, Y, label='f')
    ax.plot(X, dY, label='df')
    
    if ddf:
        ax.plot(X, ddY, label='ddf')
    ax.legend(loc='upper right') 
    ax.grid()

    if save_file:
        fig.savefig(save_file)
        print(f"Saved in {save_file} succesfuly")
        plt.close(fig) 

    plt.show()


if __name__ == '__main__':
    show_plot(func, dfunc,  [-2, 2], "newton.png", step=0.1, ddf=ddfunc)
    show_plot(phi, dphi, [0, 2], "iteration.png", step=0.1)
    print("newton")
    newton(func, dfunc, x0=0.5)
    print("simple_iteration")
    simple_iteration(phi, dphi, 0, 2)
