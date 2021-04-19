import numpy as np

from py_esm.models.methods.scf.HartreeFockProcedure import HartreeFockProcedure


def mp2_mol(mol, basis):
    """
    MP2 calculation on a molecule
    :param mol: Molecule Object
    :return: eigenvalues, eigenvectors and total energy of the system
    """

    scf = HartreeFockProcedure(mol, basis)
    scf.run_hf(25, False)

    calculator = MP2EnergyCalculator(mol.n_electrons, scf.basis.n_basis, scf.basis.met)

    return scf.E, scf.C, (scf.energy + calculator.calculate(scf.E))


class MP2EnergyCalculator():

    def __init__(self, N, dim, multi_electron_tensor):
        self.N = N
        self.spatial_dim = dim
        self.met = multi_electron_tensor

    def calculate(self, e_values):
        return self.EMP2(e_values)

    # with spatial orbital basis functions
    # e_values: orbital energies = eigenvalues of the Fock matrix
    # returns e_mp2
    def EMP2_spatial(self, e_values):
        e_mp2 = 0

        for i in range(self.spatial_dim):
            for j in range(i, self.spatial_dim):
                for x in range(self.spatial_dim):
                    for y in range(x, self.spatial_dim):
                        denum = e_values[i] + e_values[j] - e_values[x] - e_values[y]
                        if (denum != 0):
                            e_mp2 += (self.met[i, j, x, y] - self.met[i, x, j, y]) ** 2 / denum

        return e_mp2

    # with spin-spatial orbital basis functions
    def EMP2(self, e_values):

        dim = 2 * self.spatial_dim

        spin_ints = np.zeros((dim, dim, dim, dim))
        for p in range(dim):
            for q in range(dim):
                for r in range(dim):
                    for s in range(dim):
                        value1 = self.met[p // 2, r // 2, q // 2, s // 2] * (p % 2 == r % 2) * (q % 2 == s % 2)
                        value2 = self.met[p // 2, s // 2, q // 2, r // 2] * (p % 2 == s % 2) * (q % 2 == r % 2)
                        spin_ints[p, q, r, s] = value1 - value2

        fs = np.zeros((dim))
        for i in range(dim):
            fs[i] = e_values[i // 2]

        e_mp2 = 0
        for i in range(self.N):
            for j in range(self.N):
                for a in range(self.N, dim):
                    for b in range(self.N, dim):
                        e_mp2 += 0.25 * spin_ints[i, j, a, b] * spin_ints[i, j, a, b] / (fs[i] + fs[j] - fs[a] - fs[b])

        return e_mp2
