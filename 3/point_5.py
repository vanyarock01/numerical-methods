var = """y = 1 / (x^3 + 64)
x0 =  -2
xk =   2
h1 = 1.0
h2 = 0.5
"""

def runge_romberg(h1, h2, y1, y2, pow=2):
    return abs((y1 - y2) / ((h2 / h1) ** pow - 1.0))


def f(x):
    return 1 / (x ** 3 + 64)


def f_(x):
    return x / ((3*x + 4) ** 2)


def arange(x0, xk, step=0.1):
    x = x0
    r = []
    while x <= xk:
        r.append(x)
        x += step
    return r


raw_data = {
    "x0": -2,
    "xk": 2,
    "h1": 1.0,
    "h2": 0.5,
    "f": f,
    "val": 0,
    "x0_": -1,
    "xk_": 1,
    "h1_": 0.5,
    "h2_": 0.25,
    "f_": f_,
    "val_": -0.16474
}


def rectangle_method(h, x, f):
    return h * sum(
        f((x[i-1] + x[i]) / 2) for i in range(1, len(x)))


def trapezium_method(h, x, f):
    return h * (
        f(x[0]) / 2 + sum(
            f(x[i]) for i in range(1, len(x) - 1)) + f(x[len(x) - 1]))


def simpson_method(h, x, f):
    return (h / 3) * (f(x[0]) +
                      sum(4 * f(x[i]) for i in range(1, len(x)-1, 2)) +
                      sum(2 * f(x[i]) for i in range(2, len(x)-1, 2)) +
                      f(x[len(x) - 1]))


def main():
    print(var)

    h1 = raw_data['h1_']
    h2 = raw_data['h2_']
    x0 = raw_data['x0_']
    xk = raw_data['xk_']
    f = raw_data['f_']
    val = raw_data['val_']

    methods = {
        "rectangle": rectangle_method,
        "trapezium": trapezium_method,
        "simpson": simpson_method
    }

    for title, func in methods.items():
        print('\n', title, '\n')

        first = func(h1, arange(x0, xk, h1), f)
        second = func(h2, arange(x0, xk, h2), f)

        ptrn = "h{} ~~ F(x) = {:8.5f}"
        print(ptrn.format(1, first))
        print(ptrn.format(2, second))
        err = runge_romberg(h1, h2, first, second, pow=2)
        real_err = abs(first - val)
        print("error: runge = {:7.3f}, real = {:7.3f}"
            .format(err, real_err))



if __name__ == '__main__':
    main()
