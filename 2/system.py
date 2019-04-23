import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
"""
System:
{
	x1 - cos(x2) = 2,
	x2 - sin(x1) = 2
}
"""

# first function ---


def printer(data):
    for title, value in data.items():
        print("| {} == {:7.4f} ".format(title, value), end='')
    print('|')


def f1(x1, x2):
    return x1 - math.cos(x2) - 2

def x1(x2):
    return math.cos(x2) + 2

def df1_dx1(x1, x2):
    return 1


def df1_dx2(x1, x2):
    return math.sin(x2)


# second function ---

def f2(x1, x2):
    return x2 - math.sin(x1) - 2

def x2(x1):
    return math.sin(x1) + 2

def df2_dx1(x1, x2):
    return -math.cos(x1)


def df2_dx2(x1, x2):
    return 1


func = {
    "f1": f1,
    "f2": f2,

    "df1_dx1": df1_dx1,
    "df1_dx2": df1_dx2,
    "df2_dx1": df2_dx1,
    "df2_dx2": df2_dx2
}

# shape = 2 only

def det(A):
    return A[0][0] * A[1][1] - A[0][1] * A[1][0]


def getx(x1, x2, fc):
    A1 = [
        [fc["f1"](x1, x2), fc["df1_dx2"](x1, x2)],
        [fc["f2"](x1, x2), fc["df2_dx2"](x1, x2)]
    ]
    A2 = [
        [fc["df1_dx1"](x1, x2), fc["f1"](x1, x2)],
        [fc["df2_dx1"](x1, x2), fc["f2"](x1, x2)]
    ]

    J = [
        [fc["df1_dx1"](x1, x2), fc["df1_dx2"](x1, x2)],
        [fc["df2_dx1"](x1, x2), fc["df2_dx2"](x1, x2)]
    ]
    detJ = det(J)
    return x1 - det(A1) / detJ, x2 - det(A2) / detJ


def norm(x, x_prev):
    return max(x[0] - x_prev[0], x[1] - x_prev[1])


def newton(x1, x2, precision=0.01):
    x1_prev, x2_prev = x1_cur, x2_cur = x1, x2
    k = 0
    while True:
        k += 1
        x1_cur, x2_cur = getx(x1_prev, x2_prev, func)

        epsilon = norm(
            (x1_cur,  x2_cur), (x1_prev, x2_prev))
        printer(
            {"x1": x1_cur, "x2": x2_cur, "eps": epsilon})

        if epsilon <= precision:
            break
        x1_prev, x2_prev = x1_cur, x2_cur
    return x1_cur, x2_cur


def ph1(x1, x2):
    return math.cos(x2) + 2


def ph2(x1, x2):
    return math.sin(x1) + 2


def dph1_dx1(x1, x2):
    return 0


def dph1_dx2(x2):
    return -math.sin(x2)


def dph2_dx1(x1):
    return math.cos(x1)


def dph2_dx2(x1, x2):
    return 0


phi = {
    "ph1": ph1,
    "ph2": ph2,

    "dph1_dx1": dph1_dx1,
    "dph1_dx2": dph1_dx2,
    "dph2_dx1": dph2_dx1,
    "dph2_dx2": dph2_dx2
}


def getq(x1, x2, phi):
	return max(
		abs(phi["dph1_dx1"](x1, x2)) + abs(phi["dph1_dx2"](x2)),
		abs(phi["dph2_dx1"](x1)) + abs(phi["dph2_dx2"](x1, x2)))

#TODO fix jumping
def simple_iteration(x1, x2, precision=0.01):
    x1_prev, x2_prev = x1_cur, x2_cur = x1, x2
    k = 0
    
    q = getq(x1, x2, phi)
    print(q)
    if q >= 1:
    	print("leave G field")
    	printer({"q": q})
    while True:
        k += 1
        x1_cur, x2_cur = ph1(x1_prev, x2_prev), ph2(x1_prev, x2_prev) 
        
        epsilon = abs((q / (1 - q)) * norm(
            (x1_cur,  x2_cur), (x1_prev, x2_prev)))
        printer(
            {"x1": x1_cur, "x2": x2_cur, "eps": epsilon})
        if epsilon <= precision:
        	break
        x1_prev, x2_prev = x1_cur, x2_cur
    return x1_cur, x2_cur


def show_plot(f1, f2,  x, save_file, step=0.1):
    X = np.arange(x[0], x[-1], step)
    Y1 = [f1(val) for val in X]
    Y2 = [f2(val) for val in X]

    fig, ax = plt.subplots()
    ax.plot(X, Y1, label='f1')
    ax.plot(Y2, X, label='f2')
    
    ax.legend(loc='upper right') 
    ax.grid()

    if save_file:
        fig.savefig(save_file)
        print(f"Saved in {save_file} succesfuly")
        plt.close(fig) 

    plt.show()



if __name__ == '__main__':
    show_plot(x1, x2, [-3, 5], 'newton.png')
    show_plot(dph1_dx2, dph2_dx1, [-3, 5], 'iteration.png')
    print("newton")
    newton(x1=0.5, x2=2, precision=0.001)
    print("iteration")
    simple_iteration(x1=0.5, x2=2, precision=0.0001)
