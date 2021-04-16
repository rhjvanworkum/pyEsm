from collections import Counter


def sort_bend_angle(atoms):
    """
    :param atoms: list of two bonds shared by a central atom
    :return: list of 3 atoms, defining the bend angle between the 2 bonds
    """
    center = 0
    center_idx = []
    atom_list = list(atoms)
    counter = Counter(atom_list)
    for i in range(len(atom_list)):
        if counter[atom_list[i]] == 2:
            center = atom_list[i]
            center_idx.append(i)

    atom_list.pop(center_idx[1])
    atom_list.pop(center_idx[0])

    return [atom_list[0], center, atom_list[1]]


def sort_torsion_bonds(bond1, bond2, bond3):
    """
    sorts the atoms in the 3 bonds of a torsional angle
    in geometrical order
    :return: list of atoms in torsional angle
    """
    t = []
    if bond1[0] in bond2:
        t.append(bond1[1])
        t.append(bond1[0])
    else:
        t.append(bond1[0])
        t.append(bond1[1])

    if bond2[0] in t:
        t.append(bond2[1])
    else:
        t.append(bond2[0])

    if bond3[0] in t:
        t.append(bond3[1])
    else:
        t.append(bond3[0])

    return t


def sort_bend_angle_atoms(angle, bonds_array):
    """
    Sorts an angle in the order of outside atom - middle atom - outside atom
    :param angle: the angle to be sorted
    :param bonds_array: array containing the bonds of the corresponding angle atoms
    :return: sorted angle
    """
    center_id = 0

    for index, atom1 in enumerate(angle):
        bonds = bonds_array[index]

        match = [False, False]

        i = 0
        for atom2 in angle:
            if atom1 != atom2:
                for bond in bonds:
                    if atom2 in bond:
                        match[i] = True
                        i += 1

        if False not in match:
            center_id = index
            break

    if center_id != 1:
        temp = angle[1]
        angle[1] = angle[center_id]
        angle[center_id] = temp

    return angle
