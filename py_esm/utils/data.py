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