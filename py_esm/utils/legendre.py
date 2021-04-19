import numpy as np

# generates legendre polynomial
# n -> order
# x -> x coordinate
def l(n, x):
    if n == 0:
        return 1.0
    elif n == 1:
        return x
    else:
        return ((2.0 * n - 1.0) * x * l(n - 1, x) - (n - 1) * l(n - 2, x)) / n

# returns the derivative of the legendre polynomial
# n -> order
# x -> xcoordinate
def grad_l(n, x):
    if n == 0:
        return 0
    elif n == 1:
        return 1.0
    else:
        return (n / (x ** 2 - 1.0)) * (x * l(n, x) - l(n - 1, x))

# generates the roots of the legendre polynomial, which are used as the nodes for integration in gaussian quadrature
# n -> order
# tolerance -> tolerance
def l_roots(n, tolerance=1e-12):
    if n < 2:
        raise ValueError(" order should be above 2")
    else:
        roots = []

        for i in range(1, int(n / 2) + 1):
            # estimate of the roots
            x = np.cos(np.pi * (i - 0.25) / (n + 0.5))
            # set the movement well above the tolerance
            movement = 10 * tolerance

            while (movement > tolerance):
                # newton-rhapson method
                dx = -l(n, x) / grad_l(n, x)
                x = x + dx
                movement = abs(dx)

            roots.append(x)

        # putting the roots into an array and mirroring them around the center
        roots = np.array(roots)
        if n % 2 == 0:
            roots = np.concatenate((-1.0 * roots, roots[::-1]))
        else:
            roots = np.concatenate((-1.0 * roots, [0.0], roots[::-1]))

        return roots

# calculates the gaussian legendre quadrature
# f -> function to quad
# a -> lower bound
# b -> higher bound
# n -> number of nodes
def gauss_legendre_quadrature(func, a, b, n):
    roots = l_roots(n)
    weights = 2.0 / ((1.0 - roots**2) * (grad_l(n, roots) ** 2))
    c1 = (b + a) / 2.0
    c2 = (b - a) / 2.0
    return c2 * np.sum(weights * func(c1 + c2 * roots))

def leggauss(n):
    # roots = l_roots(n)
    # weights = 2.0 / ((1.0 - roots**2) * (grad_l(n, roots) ** 2))

    from numpy.polynomial.legendre import leggauss
    roots, weights = leggauss(n)

    grid = []
    for i in range(len(roots)):
        grid.append((roots[i], weights[i]))

    return grid
