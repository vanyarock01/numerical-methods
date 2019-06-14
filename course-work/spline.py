import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt


def draw_plot(X, Y, dot_x, dot_y, labels, save_file="plot.png"):
    fig, ax = plt.subplots()
    
    for i in range(len(X)):
        ax.plot(X[i], Y[i], label=f"{labels[i]}")
    ax.scatter(dot_x, dot_y, s=20)
    ax.legend(loc='upper right')
    ax.grid()
    if save_file:
        fig.savefig(save_file)
        print(f"Saved in {save_file} succesfuly")
        plt.close(fig)


def vector_print(vec, header=None):
    if header:
        print(header, '\n')
    for i in range(len(vec)):
        print('{:6.2f}'.format(float(vec[i])), end=' ')
    print('\n')


def log_info(title, func, msg):
    print(f"=={title} from {func} <==\n=={msg}")


def checkout(state, msg, u):
    if state:
        log_info("ERR", u, msg)


def exp_coeff(n, x, y, p, phi):
    """computation of the coefficients of the exponential spline"""
    u = "exp_coeff"

    # output
    d = [0.0] * (n + 1)
    dq = [0.0] * (n)
    h = [0.0] * (n)
    hp = [0.0] * (n)
    ph = [0.0] * (n)

    # computation of the elements of the tridiagonal system
    i1, n1 = 0, n - 1
    c, c1, c2, u, v, w = 0.0, 0.0, 0.0, y[0], 0.0, 0.0
    checkout(n < 3, "the value n must be > or == 3", u)
    q = [0.0] * (n + 1)
    r = [0.0] * (n + 1)
    for i in range(1, n):
        i1 = i - 1
        h[i1] = x[i] - x[i1]
        checkout(h[i1] < 0.0, "the value h[i] must be > or == 0", u)
        v = y[i]
        hp[i1] = abs(h[i1] * p[i1])
        if h[i1] == 0.0:
            d[i] = v
        else:
            d[i] = (v - u) / h[i1]
            u = v
        if hp[i1] > 0.5:
            ph[i1] = math.exp(-hp[i1])
            c = ph[i1] ** 2
            c1 = 1.0 - c
            c2 = c1 / hp[i1]
            c1 *= hp[i1] / h[i1]
            q[i] = (1.0 - c2 + c) / c1
            r[i] = (c2 - 2.0 * ph[i1]) / c1
        else:  # auxiliary function phi
            c = hp[i1] * hp[i1]
            ph[i1] = phi(c)
            w = h[i1] / (1.0 + c * ph[i1])
            c *= 0.25
            c2 = 1.0 + c * phi(c)
            q[i] = (0.5 * c2 ** 2 - ph[i1]) * w
            r[i] = ph[i1] * w

    """
    solution of the tridiagonal system with
    diagonal:        q[i] + q[i+1] i = 1 , n1
    off-diagonal:    r[i]          i = 2,  n1
    right hand side: d[i+1] - d[i] i= 1,   n1
    second difference quotient: dq[i]: = d[i+1] - d[i] i= 1, n1
    """

    d[0], u = 0.0, 0.0
    for i in range(n1):
        q[i] = q[i] + q[i + 1] - u * r[i]
        dq[i] = d[i + 1] - d[i]
        d[i] = dq[i] - u * d[i - 1]
        u = r[i + 1] / q[i]
    d[n] = 0.0
    for i in range(n1, 1, -1):
        d[i] = (d[i] - r[i + 1] * d[i + 1]) / q[i]

    return d, dq, h, hp, ph


def exp_spl(n, x, y, p, d, h, hp, ph, xx, i, phi):
    """evaluation of the exponential spline at the abscissa xx"""
    u = "exp_spl"
    # output
    ex = 0.0

    i1 = i + 1
    checkout(i > n - 1, "incorrect params", u)
    t = (xx - x[i]) / h[i]
    t1 = 1.0 - t
    checkout(t > 1.0 or t < 0.0, "incorrect params", u)
    if hp[i] > 0.5:
        e = math.exp(-t * hp[i])
        e1 = math.exp(-t1 * hp[i])
        c = 1.0 - ph[i] ** 2
        ex = y[i1] * t + y[i] * t1 + \
            (d[i1] * (e1 * (1.0 - e ** 2) / c - t) +
             d[i] * (e * (1.0 - e1 ** 2) / c - t1)) / (p[i] ** 2)
    else:
        e = t * hp[i]
        e1 = t1 * hp[i]
        c = h[i] ** 2 / (1.0 + hp[i] ** 2 * ph[i])
        ex = t * (y[i1] + d[i1] * c * (t ** 2 * phi(e ** 2) - ph[i])) + \
            t1 * (y[i] + d[i] * c * (t1 ** 2 * phi(e1 ** 2) - ph[i]))
    return ex


def phi(a):
    """auxiliary function phi with best approximation for sinh"""
    return ((0.27713991169e1 - 5 * a + 0.19840927713e1 - 3)*a +
            0.83333336379e1 - 2)*a + 0.16666666666


def generator(n, x, y, d, dq, h):
    # output
    p = [0.0] * n
    # determination of the tension parameters p
    if h[0] > 0.0:
        y11 = y[1]
    else:
        y11 = y[0]
    y1 = abs(y11)
    d1 = d[1] * dq[1]
    for i in range(n-2):  # !!!
        i1 = i + 1
        d2 = d[i1] * dq[i1]
        y2 = abs(y[i1])
        if y2 > y1:
            y1 = y2
        if d1 * d2 > 0.0:
            p[i] = 0.0
        else:
            if y1 == 0.0:
                p[i] = 15 / h[i]
            else:
                p[i] = (4 + 1 / (0.1 + abs(y[i1] - y11) / y1)) / h[i]
        y11 = y[i1]
        y1 = y2
        d1 = d2
    p[0] = p[1]
    p[n - 1] = p[n - 2]
    return p

def main():
    def get_spline(x, y, n, step, p=None):
        if not p:
            p = [0.0] * (n-1)
        d, dq, h, hp, ph = exp_coeff(n, x, y, p, phi)
        #p = generator(n, x, y, d, dq, h)
        X = []
        Y = []
        for i, point in enumerate(x[:-1]):
            X.append(x[i])
            Y.append(y[i])
            for xx in np.arange(point + step, x[i+1], step):
                X.append(xx)
                ex = exp_spl(n, x, y, p, d, h, hp, ph, xx, i, phi)
                Y.append(ex)
        return X, Y


    tests = [
        {   
            'n': 8,
            'x': [-6, 1,   3,   6,   8, 10, 11, 12],
            'y': [-2, 2, 3.5, 3.5, 2.8, -4, 2.8, 5],
            'p': [0, 0, 1, 3.6, 0, 0, 0],
            'step': 0.1
        },
        {
            'n': 12,
            'x': [-11, -10, -9, -6, -5, -2.5, 0, 5, 8, 9, 10, 11],
            'y': [0, 1.4, 1.6, 1.7, 2, 4, 4, 4, 2, 1.7, 1.4, 0],
            'p': [8.4, 8.4, 3.4, 0, 0, 5.6, 2.8, 0, 8, 7.6, 7.6],
            'step': 0.1
        },
        {
            'n': 5,
            'x': [-5, -2, -0.5, 3, 8],
            'y': [2, 4, -1, 5, 1],
            'p': [0, 1, 1, 0],
            'step': 0.01
        },
        {
            'n': 7,
            'x': [0, 1, 1.2, 2.5, 3.8, 4, 5],
            'y': [-2, 0, 5, 0, 5, 0, -2],
            'p': [0, 5, 25, 25,25, 25],
            'step': 0.01
        }
    ]


    for idx, test in enumerate(tests):
        print(test)
        X, Y, labels = [], [], []
        cs = CubicSpline(np.array(test['x']), np.array(test['y']))
        xs = np.arange(test['x'][0], test['x'][-1], test['step'])

        x, y = get_spline(test['x'], test['y'], test['n'], test['step'], test['p'])

        X.append(xs)
        Y.append(cs(xs))
        labels.append('cubic_spline ')

        X.append(x)
        Y.append(y)
        labels.append('exp_spline')

        draw_plot(X, Y, test['x'], test['y'], labels, f'plot{idx}.png')


if __name__ == '__main__':
    main()