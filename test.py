from py_esm.helpers.basis_set import get_basis_functions
from py_esm.helpers.gaussian_integrals import s_integral
from py_esm.models.Molecule import Molecule
from py_esm.models.methods.hartree_fock import hf_mol

import numpy as np

mol = Molecule('O')

mol.set_geometry(1.1, [8, 1])
mol.set_geometry(104.0 / 180 * np.pi, [1, 8, 1])

evals, evecs, e = hf_mol(mol)

print(e)


#
# Ra = mol.atoms[orbital['atom_index']].coordinates
# coeffs_vec_m = orbital['coeffs']
# zeta = get_zeta_value(orbital['atom'])
# alpha_vec_m = orbital['alphas']
#
# s = 0
#
# for p in range(len(orbital['coeffs'])):
#
#     Rb = mol.atoms[orbital2['atom_index']].coordinates
#     coeffs_vec_n = orbital2['coeffs']
#     zeta = get_zeta_value(orbital2['atom'])
#     alpha_vec_n = orbital2['alphas']
#
#     for q in range(len(orbital2['coeffs'])):
#         print(coeffs_vec_m[p] * coeffs_vec_n[q] * \
#                                  s_integral(alpha_vec_m[p], Ra, alpha_vec_n[q], Rb))
#         s += coeffs_vec_m[p] * coeffs_vec_n[q] * \
#                                  s_integral(alpha_vec_m[p], Ra, alpha_vec_n[q], Rb)
#
# print(s)