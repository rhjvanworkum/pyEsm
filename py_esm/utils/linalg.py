import numpy as np


def get_2d_matrix_index(a, b):
    if a > b:
        ab = a*(a+1)/2 + b
    else:
        ab = b*(b+1)/2 + a

    return ab


def get_4d_matrix_index(a, b, c, d):
    ab = get_2d_matrix_index(a, b)
    cd = get_2d_matrix_index(c, d)

    if ab > cd:
        abcd = ab*(ab+1)/2 + cd
    else:
        abcd = cd*(cd+1)/2 + ab

    return abcd


def list_to_array(list):
    """
    :param list:
    :return:
    """

    arr = np.zeros(len(list))
    for i in range(len(list)):
        arr[i] = list[i]

    return arr


def get_magnitude(vec):
    """
    :param vec: Vector
    :return: Returns the magnitude of a vector
    """
    return vec[0] ** 2 + vec[1] ** 2 + vec[2] ** 2


def get_rotational_matrix(theta, r):
    """
    :param theta: angle in radians to rotate
    :param r: axis around which to rotate
    :return: Rotational matrix of theta radians around axis r
    """
    cos = np.cos(theta)
    icos = 1 - np.cos(theta)
    sin = np.sin(theta)

    return np.array([
        [cos + r[0] ** 2 * icos, r[0] * r[1] * icos - r[2] * sin, r[0] * r[2] * icos + r[1] * sin],
        [r[1] * r[0] * icos + r[2] * sin, cos + r[1] ** 2 * icos, r[1] * r[2] * icos - r[0] * sin],
        [r[2] * r[0] * icos - r[1] * sin, r[2] * r[1] * icos + r[0] * sin, cos + r[2] ** 2 * icos]
    ])
