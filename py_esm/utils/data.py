import numpy as np
from itertools import permutations, product

from py_esm.utils.linalg import get_2d_matrix_index, get_4d_matrix_index

def xyz_reader(file_name):
    """
    :param file_name: absolute path to the xyz file
    :return: natoms, list of atomic symbols, list of atomic coordinates
    """
    file = open(file_name, 'r')

    number_of_atoms = 0
    atom_types = []
    atom_coordinates = []

    for index, line in enumerate(file):
        if index == 0:
            try:
                number_of_atoms = int(line.split()[0])
            except:
                print('your file does not seem compatible')

        if index >= 1:
            split = line.split()
            atom_types.append(split[0])
            atom_coordinates.append([float(split[1]), float(split[2]), float(split[3])])

    file.close()

    return number_of_atoms, atom_types, atom_coordinates

def read_2d_matrix(file_name):
    """
    :param file_name:
    :return:
    """

    file = open(file_name, 'r')

    for index, line in enumerate(file):
        if index == 0:
            try:
                indices = int(line.split()[0])
                array = np.zeros((indices, indices), dtype=np.float64)
            except:
                print('your file does not seem compatible')

        else:
            split = line.split()

            idx1 = int(split[0]) - 1
            idx2 = int(split[1]) - 1

            if len(split) == 3:
                array[idx1][idx2] = float(split[2])
                array[idx2][idx1] = float(split[2])
            else:
                idx3 = int(split[2]) - 1
                idx4 = int(split[3]) - 1

                array[get_4d_matrix_index(idx1, idx2, idx3, idx4)] = float(split[4])

    return array

def read_4d_matrix(file_name):
    file = open(file_name, 'r')

    for index, line in enumerate(file):
        if index == 0:
            try:
                indices = int(line.split()[0])
                array = {}
            except:
                print('your file does not seem compatible')

        else:
            split = line.split()

            idx1 = int(split[0]) - 1
            idx2 = int(split[1]) - 1
            idx3 = int(split[2]) - 1
            idx4 = int(split[3]) - 1

            array[get_4d_matrix_index(idx1, idx2, idx3, idx4)] = float(split[4])

    for p in product(np.arange(0, 7), repeat=4):
        if get_4d_matrix_index(p[0], p[1], p[2], p[3]) not in array.keys():
            array[get_4d_matrix_index(p[0], p[1], p[2], p[3])] = 0.0

    return array

def print_sym_matrix(mat):
    for i in range(len(mat)):
        for j in range(0, i + 1):
            print(i + 1, j + 1, mat[i][j])
