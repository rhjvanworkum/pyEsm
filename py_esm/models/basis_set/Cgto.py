import numpy as np
from scipy.special import factorial2 as fact2
from scipy.special import hyp1f1


class ContractedGaussianOrbital:
    """
    Data object of a contracted gaussian type orbital(CGTO)
    """

    def __init__(self, origin, shell, ang_mom_numbers, exponents, coeffs):
        self.origin = origin
        self.shell = shell
        self.ang_mom = ang_mom_numbers
        self.exponents = exponents
        self.coeffs = coeffs
        self.n_gaussians = len(exponents)
        self.norm = None

        self.normalize()


    def __call__(self, x, y, z):
        i, j, k = self.ang_mom
        dx, dy, dz = x - self.origin[0], y - self.origin[1], z - self.origin[2]
        d2 = dx ** 2 + dy ** 2 + dz ** 2

        value = 0.0
        for idx, coef in enumerate(self.coeffs):
            value += self.norm[idx] * coef * dx ** i * dy ** j * dz ** k * np.exp(-self.exponents[idx] * d2)

        return value


    def normalize(self):
        """
        Function to normalize all the basis_set functions in the object
        - Formula's of the PGTO / CGTO normalization can be found at:
          https://arxiv.org/ftp/arxiv/papers/2007/2007.12057.pdf
        """

        L = np.sum(self.ang_mom)

        # Normalization of the gaussian primitives
        assert isinstance(self.ang_mom, tuple)
        self.norm = np.sqrt(np.power(2, 2 * L + 1.5) *
                        np.power(self.exponents, L + 1.5) /
                        fact2(2 * self.ang_mom[0] - 1) /
                        fact2(2 * self.ang_mom[1] - 1) /
                        fact2(2 * self.ang_mom[2] - 1) /
                        np.power(np.pi, 1.5))

        # Normalization fo the contracted gaussians
        prefactor = np.power(np.pi, 1.5) * fact2(2 * self.ang_mom[0] - 1) * \
                                           fact2(2 * self.ang_mom[1] - 1) * \
                                           fact2(2 * self.ang_mom[2] - 1) / np.power(2.0, L)

        N = 0.0
        for a in range(self.n_gaussians):
            for b in range(self.n_gaussians):
                N += self.norm[a] * self.norm[b] * self.coeffs[a] * self.coeffs[b] / \
                     np.power(self.exponents[a] + self.exponents[b], L + 1.5)

        N *= prefactor
        N = np.power(N, -0.5)

        self.norm *= N
        # for a in range(self.n_gaussians):
        #     self.coeffs[a] *= N


def E(a, b, i, j, t, Xab):
    """
    Recursive definition of the Expansion coeffiecents for Hermite Gaussians
    :param a: exponent of gaussian 1
    :param b: exponent of gaussian 2
    :param i: ang_mom_number of gaussian 1
    :param j: ang_mom_number of gaussian 2
    :param t: index of hermite gaussian
    :param Xab: distance between origins of gaussian 1 & 2
    :return: E_t(ij)
    """
    p = a + b
    mu = a * b / p

    # out of bounds for t
    if (t > (i + j)) or (t < 0):
        return 0.0
    # starting case
    elif i == j == t == 0:
        return np.exp(-mu * Xab ** 2)  # K_AB
    # j index is equivalent to starting case -> decrement index i
    elif j == 0:
        return (1 / (2 * p)) * E(a, b, i - 1, j, t - 1, Xab) + \
               (mu * Xab / a) * E(a, b, i - 1, j, t, Xab) + \
               (t + 1) * E(a, b, i - 1, j, t + 1, Xab)
    # decrement index j
    else:
        return (1 / (2 * p)) * E(a, b, i, j - 1, t - 1, Xab) + \
               (mu * Xab / b) * E(a, b, i, j - 1, t, Xab) + \
               (t + 1) * E(a, b, i, j - 1, t + 1, Xab)


def overlap(a, ijk, A, b, lmn, B):
    """
    Calculates the Overlap Integral between two gaussians
    :param a: exponent of gaussian 1
    :param ijk: tuple containing the ang_mom_numbers of gaussian 1
    :param A: center of gaussian 1
    :param b: exponent of gaussian 2
    :param lmn: tuple containing the ang_mom_numbers of gaussian 2
    :param B: center of gaussian 2
    :return: overlap integral
    """
    i, j, k = ijk
    l, m, n = lmn

    S_x = E(a, b, i, l, 0, A[0] - B[0])
    S_y = E(a, b, j, m, 0, A[1] - B[1])
    S_z = E(a, b, k, n, 0, A[2] - B[2])

    return S_x * S_y * S_z * (np.pi / (a + b)) ** 1.5


def S(Ga : ContractedGaussianOrbital, Gb : ContractedGaussianOrbital) -> float:
    """
    Evalues the overlap integral between two CGTO functions
    :param Ga: ContractedGaussianOrbital 1
    :param Gb: ContractedGaussianOrbital 2
    :return: float
    """

    s = 0.0
    for index1, coeffs1 in enumerate(Ga.coeffs):
        for index2, coeffs2 in enumerate(Gb.coeffs):
            # apply the PGTO normalization constant here
            s += Ga.norm[index1] * Gb.norm[index2] * coeffs1 * coeffs2 * \
                 overlap(Ga.exponents[index1], Ga.ang_mom, Ga.origin, Gb.exponents[index2], Gb.ang_mom, Gb.origin)

    return s


def kinetic(a, ijk, A, b, lmn, B):
    """
    Calculates the Kinetic Integral between two gaussians
    :param a: exponent of gaussian 1
    :param ijk: tuple containing the ang_mom_numbers of gaussian 1
    :param A: center of gaussian 1
    :param b: exponent of gaussian 2
    :param lmn: tuple containing the ang_mom_numbers of gaussian 2
    :param B: center of gaussian 2
    :return: kinetic integral
    """

    l, m, n = lmn

    term0 = b * (2 * (l + m + n) + 3) * overlap(a, ijk, A, b, lmn, B)
    term1 = -2 * b ** 2 * (overlap(a, ijk, A, b, (l + 2, m, n), B) +
                           overlap(a, ijk, A, b, (l, m + 2, n), B) +
                           overlap(a, ijk, A, b, (l, m, n + 2), B))
    term2 = -0.5 * (l * (l + 1) * overlap(a, ijk, A, b, (l - 2, m, n), B) +
                    m * (m + 1) * overlap(a, ijk, A, b, (l, m - 2, n), B) +
                    n * (n + 1) * overlap(a, ijk, A, b, (l, m, n - 2), B))

    return term0 + term1 + term2

def T(Ga : ContractedGaussianOrbital, Gb : ContractedGaussianOrbital) -> float:
    """
    Evalues the kinetic integral between two CGTO functions
    :param Ga: ContractedGaussianOrbital 1
    :param Gb: ContractedGaussianOrbital 2
    :return: float
    """

    t = 0.0
    for index1, coeffs1 in enumerate(Ga.coeffs):
        for index2, coeffs2 in enumerate(Gb.coeffs):
            t += Ga.norm[index1] * Gb.norm[index2] * coeffs1 * coeffs2 * \
                 kinetic(Ga.exponents[index1], Ga.ang_mom, Ga.origin, Gb.exponents[index2], Gb.ang_mom, Gb.origin)

    return t


def boys_function(n, x):
    """
    the Boys Function
    :param n: order of the boys function
    :param x: argument
    :return:
    """
    return hyp1f1(n + 0.5, n + 1.5, -x) / (2.0 * n + 1.0)


def R(t, u, v, n, p, PCx, PCy, PCz, Xpc) -> float:
    """
    Hermite Coulomb Integrals through the recursive definition
    :param t: order of Hermite derivative in x
    :param u: order of Hermite derivative in y
    :param v: order of Hermite derivative in z
    :param n: order of the boys function
    :param p: exponent of the gaussian composite
    :param PCx: cartesian x distance between P & C
    :param PCy: cartesian y distance between P & C
    :param PCz: cartesian z distance between P & C
    :param Xpc: distance between Gaussian composite center P and nuclear center C
    :return:
    """

    x = p * Xpc ** 2

    value = 0.0
    # TODO : Why do we keep adding to the value here while recursing at the same time

    # starting case
    if t == u == v == 0:
        value += (-2 * p) ** n * boys_function(n, x)
    # t, u index is equivalent to starting case -> decrement index v
    elif t == u == 0:
        if v > 1:
            value += (v - 1) * R(t, u, v - 2, n + 1, p, PCx, PCy, PCz, Xpc)
        value += PCz * R(t, u, v - 1, n + 1, p, PCx, PCy, PCz, Xpc)
    # t index is equivalent to starting case -> decrement index u
    elif t == 0:
        if u > 1:
            value += (u - 1) * R(t, u - 2, v, n + 1, p, PCx, PCy, PCz, Xpc)
        value += PCy * R(t, u - 1, v, n + 1, p, PCx, PCy, PCz, Xpc)
    # decrement index t
    else:
        if t > 1:
            value += (t - 1) * R(t - 2, u, v, n + 1, p, PCx, PCy, PCz, Xpc)
        value += PCx * R(t - 1, u, v, n + 1, p, PCx, PCy, PCz, Xpc)

    return value


def nuclear_attraction(a, ijk, A, b, lmn, B, C):
    """
    Calculates the Nuclear Attraction Integral between two gaussians
    :param a: exponent of gaussian 1
    :param ijk: tuple containing the ang_mom_numbers of gaussian 1
    :param A: center of gaussian 1
    :param b: exponent of gaussian 2
    :param lmn: tuple containing the ang_mom_numbers of gaussian 2
    :param B: center of gaussian 2
    :param C: center of the nucleus
    :return: nuclear attraction integral
    """

    i, j, k = ijk
    l, m, n = lmn

    p = a + b
    P = (a * A + b * B) / p
    Xpc = np.linalg.norm(P - C)

    AB = A - B
    PC = P - C

    value = 0.0
    for t in range(i + l + 1):
        for u in range(j + m + 1):
            for v in range(k + n + 1):
                value += E(a, b, i, l, t, AB[0]) * E(a, b, j, m, u, AB[1]) * E(a, b, k, n, v, AB[2]) * \
                         R(t, u, v, 0, p, PC[0], PC[1], PC[2], Xpc)

    value *= 2 * np.pi / p
    return value


def V(Ga : ContractedGaussianOrbital, Gb : ContractedGaussianOrbital, C : np.ndarray) -> float:
    """
    Evalues the nuclear attraction integral between two CGTO functions
    :param Ga: ContractedGaussianOrbital 1
    :param Gb: ContractedGaussianOrbital 2
    :param C: array containing the nucleus coordinates
    :return: float
    """

    v = 0.0
    for index1, coeffs1 in enumerate(Ga.coeffs):
        for index2, coeffs2 in enumerate(Gb.coeffs):
            v += Ga.norm[index1] * Gb.norm[index2] * coeffs1 * coeffs2 * \
                 nuclear_attraction(Ga.exponents[index1], Ga.ang_mom, Ga.origin, Gb.exponents[index2], Gb.ang_mom, Gb.origin, C)

    return v


def electron_repulsion(a, ijk, A, b, lmn, B, c, xyz, C, d, rst, D):
    """
    Electron repulsion integral
    """

    i, j, k = ijk
    l, m, n = lmn
    x, y, z = xyz
    r, s, t2 = rst

    p = a + b
    q = c + d
    alpha = p * q / (p + q)

    P = (a * A + b * B) / p
    Q = (c * C + d * D) / q
    RPQ = np.linalg.norm(P - Q)

    AB = A - B
    CD = C - D
    PQ = P - Q

    value = 0.0
    for t in range(i + l + 1):
        for u in range(j + m + 1):
            for v in range(k + n + 1):
                for tau in range(x + r + 1):
                    for nu in range(y + s + 1):
                        for phi in range(z + t2 + 1):
                            value += E(a, b, i, l, t, AB[0]) * \
                                     E(a, b, j, m, u, AB[1]) * \
                                     E(a, b, k, n, v, AB[2]) * \
                                     E(c, d, x, r, tau, CD[0]) * \
                                     E(c, d, y, s, nu, CD[1]) * \
                                     E(c, d, z, t2, phi, CD[2]) * \
                                     (-1) ** (tau + nu + phi) * \
                                     R(t + tau, u + nu, v + phi, 0, alpha, PQ[0], PQ[1], PQ[2], RPQ)

    value *= (2 * np.pi ** 2.5 / (p * q * np.sqrt(p + q)))
    return value


def ERI(Ga : ContractedGaussianOrbital, Gb : ContractedGaussianOrbital,
        Gc : ContractedGaussianOrbital, Gd : ContractedGaussianOrbital) -> float:
    """
    Returns multi electron integral
    :param Ga:
    :param Gb:
    :param Gc:
    :param Gd:
    :return:
    """

    eri = 0.0
    for index1, coeffs1 in enumerate(Ga.coeffs):
        for index2, coeffs2 in enumerate(Gb.coeffs):
            for index3, coeffs3 in enumerate(Gc.coeffs):
                for index4, coeffs4 in enumerate(Gd.coeffs):
                    eri += Ga.norm[index1] * Gb.norm[index2] * Gc.norm[index3] * Gd.norm[index4] * \
                           coeffs1 * coeffs2 * coeffs3 * coeffs4 * \
                           electron_repulsion(Ga.exponents[index1], Ga.ang_mom, Ga.origin,
                                              Gb.exponents[index2], Gb.ang_mom, Gb.origin,
                                              Gc.exponents[index3], Gc.ang_mom, Gc.origin,
                                              Gd.exponents[index4], Gd.ang_mom, Gd.origin)
    return eri
