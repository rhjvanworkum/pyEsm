from py_esm.models.Molecule import Molecule
from py_esm.models.basis_set.CgtoBasisSet import CgtoBasisSet
# from py_esm.models.methods.ccsd import ccsd_mol
from py_esm.models.methods.ff import ff_mol
from py_esm.models.methods.dft import dft_mol
from py_esm.models.methods.hartree_fock import hf_mol

# import numpy as np
#
# a = np.array([
#     [1, 2],
#     [3, 4]
# ])
#
# b = np.array([
#     [4, 4],
#     [5, 5]
# ])
#
# P = 0
#
# for i in range(2):
#     for j in range(2):
#         P += a[i, j] * b[i, j]
#
# print(P)
# print(np.einsum('pq,qp', b, a))

mol = Molecule('O')
mol.set_atom_positions(
    [0, 1, 2],
    [
        [0.000000, -0.075791844, 0.000000],
        [0.866811829, 0.601435779, 0.000000],
        [-0.866811829, 0.601435779, 0.000000]
    ]
)

basis = CgtoBasisSet(mol)

e = dft_mol(mol, basis)

print(e)

