

var = """x* = 1
x: 0.0 0.5 1.0 1.5 2.0
y: 1.0 1.3776 1.5403 1.5707 1.5839
"""

raw_data = {
    "x*": 1.0,
    "x": [0.0, 0.5, 1.0, 1.5, 2.0],
    "y": [1.0, 1.3776, 1.5403, 1.5707, 1.5839],
    "x_": [0.0, 0.1, 0.2, 0.3, 0.4],
    "y_": [1.0, 1.1052, 1.2214, 1.3499, 1.4918],
    "x*_": 0.2
}


def derivative(x, y, val):
    n = len(x)
    k = 0
    try:
        while x[k + 1] < val:
            k += 1
    except IndexError:
        print("value does not fall within the interval")
        exit(1)

    if k + 1 > n + 1:
        k -= 1
    first = (y[k+1] - y[k])/(x[k+1] - x[k])
    second = first + \
        ((y[k+2] - y[k+1])/(x[k+2] - x[k+1]) - (y[k+1] - y[k])/(x[k+1] - x[k])) * \
        (2 * val - x[k] - x[k+1]) / \
        (x[k+2] - x[k])

    return first, second


def main():
    x, y, val = raw_data['x'], raw_data['y'], raw_data['x*']
    print(derivative(x, y, val))


if __name__ == '__main__':
    main()
