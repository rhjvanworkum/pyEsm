import numpy as np

from py_esm.helpers.forcefield import get_parameter


def ff_mol(mol):
    # check bonded forces
    E_bonds = bond_forces(mol)
    E_angles = angle_forces(mol)
    E_torsions_p = torsion_p_forces(mol)
    # E_torsions_im = torsion_im_forces(mol)

    # check non-bonded forces
    E_vdw = vdw_forces(mol)
    E_el = el_forces(mol)

    return E_bonds + E_angles + E_torsions_p + E_vdw + E_el


def bond_forces(mol):
    E = 0
    for bond in mol.bonds:
        p = get_parameter(mol, 'Bonds', bond)
        r = mol.get_bond_length(bond)
        E += 0.5 * p[0] * (r - p[1])**2
    return E


def angle_forces(mol):
    E = 0
    for angle in mol.angles:
        p = get_parameter(mol, 'Angles', angle)
        E += 0.5 * p[0] * (mol.get_bend_angle(angle) - p[1])**2
    return E


def torsion_p_forces(mol):
    E = 0
    for torsion in mol.torsions:
        p = get_parameter(mol, 'ProperTorsions', torsion)
        tau = mol.get_torsion_angle(torsion)
        for i in range(len(p[1])):
            E += p[0][i].value_in_unit(p[0][i].unit) * (1 + np.cos(p[1][i] * tau -
                                                                   (p[2][i].value_in_unit(p[2][i].unit)) * np.pi / 180))
    return E


#
# # def torsion_im_forces(torsions):
# #     return 0


def vdw_forces(mol):
    E = 0
    for i in range(0, mol.natoms - 1):
        for j in range(i + 1, mol.natoms):
            p = get_parameter(mol, 'vdW', [i, j])
            if mol.get_bond_length([i, j]) < p[2]:
                rij_6 = (1 / mol.get_bond_length([i, j])) ** 6
                E += p[0] * ((p[1] / rij_6) + (p[1] / rij_6**2))
    return E


def el_forces(mol):
    E = 0
    p = get_parameter(mol, 'el', [])
    for a in mol.atoms:
        for b in mol.atoms:
            if a != b and np.linalg.norm(a.coordinates - b.coordinates) < p[0]:
                E += (a.charge * b.charge) / np.linalg.norm(a.coordinates - b.coordinates)

    return 0.5 * E
