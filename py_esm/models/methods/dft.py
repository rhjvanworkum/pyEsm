import numpy as np
import scipy

from py_esm.helpers.dft import get_xc, LSDA
from py_esm.models.basis_set.CgtoBasisSet import CgtoBasisSet
from py_esm.models.grid.Grid import Grid
from py_esm.models.methods.scf.KohnShamProcedure import KohnShamProcedure


def dft_mol(mol, basis):

    scf = KohnShamProcedure(mol, basis)
    scf.run_hf(1e-6)

    return scf.E, scf.C, scf.energy


# def run_dft(mol, N, basis, e_nuc, n_scf):
#
#     # intialize variables
#     evals = np.zeros(basis.n_basis)
#     evecs = np.zeros((basis.n_basis, basis.n_basis))
#     e_tot = 0
#
#     B = basis.n_basis
#
#     # T and V are not affected by the iterative process
#     H_core = basis.kinetic + basis.extern
#
#     grid = Grid(mol.atoms)
#     grid.set_func(basis.basis_functions)
#
#     P = np.zeros((B, B))
#
#     for n in range(n_scf):
#
#         print('cycle: ', n)
#
#         J = np.zeros((B, B), dtype=np.float64)
#         for i in range(B):
#             for j in range(B):
#                 for x in range(B):
#                     for y in range(B):
#                         J[i, j] += P[x, y] * basis.coulomb([i, j, x, y])
#
#         e_xc, v_xc = get_xc(grid, 0.5 * P)
#
#         Fock = H_core + 2 * J + v_xc
#
#         # print(P)
#
#         evals, evecs = scipy.linalg.eigh(Fock, basis.S)
#
#         P = np.zeros((B, B))
#         for i in range(B):
#             for j in range(B):
#                 for a in range(int(N / 2)):
#                     P[i, j] += evecs[i, a] * evecs[j, a]
#
#         # calculate total energy
#         # nucleur repulsion energy
#         e_tot = e_nuc + e_xc
#         for i in range(B):
#             for j in range(i, B):
#                 e_tot += 2 * P[i, j] * (H_core[i, j] + J[i, j])
#
#         print(e_tot)
#
#     return evals, evecs, e_tot


# def build_rho(N, evecs, grid, basis):
#     rho = 0.0
#
#     occups = [[1], [1], [2], [2], [4]]
#
#     for index in range(len(evecs)):
#
#         orbitals_grid = np.dot(evecs[index], basis.func) / (np.sqrt(4 * np.pi)) / grid.points
#
#         # plt.plot(np.linspace(0, 5, len(grid.points)), orbitals_grid ** 2)
#         # plt.plot(np.linspace(0, 5, len(grid.points)), basis.func[index])
#         # plt.show()
#
#         rho += occups[index] * orbitals_grid ** 2
#
#     return rho

#
# def solve_poisson_spectral(grid, rho):
#     """Solve the radial poisson equation for a spherically symmetric density."""
#     norm = grid.integrate(4 * np.pi * grid.points**2 * rho)
#
#     int1 = grid.antiderivative(grid.points * rho)
#     int1 -= int1[0]
#     int2 = grid.antiderivative(int1)
#     int2 -= int2[0]
#
#     v_h = -(4 * np.pi) * int2
#
#     alpha = (norm - v_h[-1]) / grid.points[-1]
#     v_h += alpha * grid.points
#     v_h /= grid.points
#
#     e_h = 0.5 * grid.integrate(v_h * rho * 4 * np.pi * grid.points ** 2)
#
#     return e_h, v_h


# def run_dft(natoms, N, basis, grid, e_nuc, n_scf):
#     B = basis.nbasis
#
#     # T and V are not affected by the iterative process
#     Hcore = basis.kinetic + basis.extern
#
#     v_xc = np.zeros(len(grid.points))
#     v_hartree = np.zeros(len(grid.points))
#
#     for n in range(n_scf):
#
#         Fock = Hcore + basis.potential(v_xc + v_hartree)
#
#         evals, evecs = scipy.linalg.eigh(Fock, basis.S)
#
#         rho = build_rho(N, evecs, grid, basis)
#
#         exc, v_xc = xcfunctional(rho, LSDA)
#
#         # e_hartree, v_hartree = hartreefunctional(rho, grid)
#         e_hartree, v_hartree = solve_poisson_spectral(grid, rho)
#
#         P = np.zeros((B, B))
#         for i in range(B):
#             for j in range(B):
#                 for a in range(int(N / 2)):
#                     P[i, j] += 2 * evecs[i, a] * evecs[j, a]
#
#         # calculate total energy
#         # nucleur repulsion energy
#         e_tot = e_nuc
#         for i in range(B):
#             for j in range(i, B):
#                 e_tot += 0.5 * P[i, j] * (Hcore[i, j] + Fock[i, j])
#
#         # volume_element = 4 * np.pi * grid.points ** 2
#         # print(e_tot)
#         # print(np.sum(evals) - e_hartree - grid.integrate(rho * exc * volume_element) - grid.integrate(v_xc * rho * volume_element))
#
#     return (evals, evecs, e_tot)