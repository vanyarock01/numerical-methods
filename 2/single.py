import math


def func(x):
    return math.exp(x) - x ** 3 + 3 * (x ** 2) - 2 * x - 3


def dfunc(x):
    return math.exp(x) - 3 * (x ** 2) + 6 * x - 2


def phi(x, f, df):
    return math.log(x ** 3 - 3 * (x ** 2) + 2 * x + 3)


def dphi(x, df):
    return (3 * (x ** 2) - 6 * x + 2) / \
        (x ** 3 - 3 * (x ** 2) + 2 * x + 3)


def newton(f, df, x0, precision=0.0001):
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


def simple_iteration(phi, dphi, a, b, precision=0.000   1):
    q = max(abs(dphi(a, dfunc)), abs(dphi(b, dfunc)))
    x = (b - a) / 2
    x_next = x
    k = 0
    cov = True
    while cov:
        k += 1
        x_next = phi(x, func, dfunc)

        print(f"x: {x_next} k: {k} |x_next - x|: {abs(x_next - x)}")
        if abs(x_next - x) <= precision:
            cov = False
        x = x_next
        if k == 10:
            break


if __name__ == '__main__':
    print("newton")
    newton(func, dfunc, x0=0.5)
    print("simple_iteration")
    simple_iteration(phi, dphi, -2, 2)
