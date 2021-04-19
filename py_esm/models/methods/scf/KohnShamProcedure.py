import numpy as np

from py_esm.helpers.dft import get_xc
from py_esm.models.grid.Grid import Grid


class KohnShamProcedure:

    def __init__(self, mol, basis):
        self.basis = basis
        self.mol = mol

        self.E = np.zeros(self.basis.n_basis)
        self.C = np.zeros((self.basis.n_basis, self.basis.n_basis))
        self.nuc_energy = mol.e_nuc_repulsion
        self.el_energy = 0
        self.energy = 0

        self.P_previous = np.zeros((self.basis.n_basis, self.basis.n_basis), dtype=np.float64)
        self.P = np.zeros((self.basis.n_basis, self.basis.n_basis), dtype=np.float64)

        e_s, U = np.linalg.eig(self.basis.overlap)
        diag_s = np.diag(e_s ** -0.5)
        self.X = np.dot(U, np.dot(diag_s, U.T))

    def run_hf(self, tol, max_steps=1000, DIIS=True):
        B = self.basis.n_basis
        # T and V are not affected by the iterative process
        H_core = self.basis.kinetic + self.basis.extern

        grid = Grid(self.mol.atoms)
        grid.set_func(self.basis.basis_functions)

        if DIIS:
            self.fock_set = []
            self.error_set = []

        for n in range(max_steps):

            e_xc, v_xc = get_xc(grid, 0.5 * self.P)

            J = np.zeros((B, B), dtype=np.float64)
            for i in range(B):
                for j in range(B):
                    for x in range(B):
                        for y in range(B):
                            J[i, j] += self.P[x, y] * self.basis.coulomb([i, j, x, y])

            Fock = H_core + 2 * J + v_xc

            if DIIS and n > 0:
                F_diis = self.updateDIIS(Fock, self.P)
                Fock_prime = self.toOrthonormal(F_diis)
            else:
                Fock_prime = self.toOrthonormal(Fock)

            evals_prime, C_prime = np.linalg.eig(Fock_prime)

            indices = evals_prime.argsort()

            evals_prime = evals_prime[indices]
            self.E = evals_prime

            C_prime = C_prime[:, indices]
            self.C = np.dot(self.X, C_prime)

            self.P_previous = self.P
            self.P = np.zeros((B, B))
            for i in range(B):
                for j in range(B):
                    for a in range(int(self.mol.n_electrons / 2)):
                        self.P[i, j] += self.C[i, a] * self.C[j, a]

            self.energy = 0.0
            self.el_energy = e_xc
            for i in range(B):
                for j in range(B):
                    self.el_energy += 2 * self.P[i, j] * (H_core[i, j] + J[i, j])
            self.energy = self.nuc_energy + self.el_energy

            if np.linalg.norm(self.P - self.P_previous) < tol:
                break

    def toOrthonormal(self, mat):
        return np.dot(self.X.T, np.dot(mat, self.X))

    # def toAO(self, mat):
    #     return np.dot

    def updateDIIS(self, Fock, P):
        FPS = np.dot(Fock, np.dot(P, self.basis.overlap))
        SPF = np.conjugate(FPS).T

        error = np.dot(self.X, np.dot(FPS - SPF, self.X))

        self.fock_set.append(Fock)
        self.error_set.append(error)
        n_fock = len(self.fock_set)

        if n_fock > 8:
            del self.fock_set[0]
            del self.error_set[0]
            n_fock -= 1

        B = np.zeros((n_fock + 1, n_fock + 1))
        B[-1, :] = B[:, -1] = 1.0
        B[-1, -1] = 0

        for i in range(n_fock):
            for j in range(i + 1):
                B[i, j] = B[j, i] = np.real(np.trace(np.dot(np.conjugate(self.error_set[i]).T, self.error_set[j])))

        residual = np.zeros(n_fock + 1)
        residual[-1] = 1.0
        weights = np.linalg.solve(B, residual)

        assert np.isclose(sum(weights[:-1]), 1.0)

        F = np.zeros((self.basis.n_basis, self.basis.n_basis))
        for i, fock in enumerate(self.fock_set):
            F += weights[i] * fock

        return F
