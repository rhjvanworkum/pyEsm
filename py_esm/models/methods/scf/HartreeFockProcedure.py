import numpy as np


class HartreeFockProcedure:

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

    def run_hf(self, tol, max_steps=100, DIIS=False, report_iterations=False):
        B = self.basis.n_basis
        # T and V are not affected by the iterative process
        H_core = self.basis.kinetic + self.basis.extern

        if DIIS:
            self.fock_set = []
            self.error_set = []

        if report_iterations:
            self.MO_list = []
            self.C_list = []

        for n in range(max_steps):

            G = np.zeros((B, B), dtype=np.float64)
            for i in range(B):
                for j in range(B):
                    for x in range(B):
                        for y in range(B):
                            G[i, j] += self.P[x, y] * (2 * self.basis.coulomb([i, j, x, y]) - self.basis.exchange([i, x, j, y]))

            Fock = H_core + G

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

            if report_iterations:
                self.MO_list.append(self.E)
                self.C_list.append(self.C)

            self.P_previous = self.P
            self.P = np.zeros((B, B))
            for i in range(B):
                for j in range(B):
                    for a in range(int(self.mol.n_electrons / 2)):
                        self.P[i, j] += self.C[i, a] * self.C[j, a]

            self.el_energy = 0.0
            self.energy = 0.0
            for i in range(B):
                for j in range(B):
                    self.el_energy += self.P[i, j] * (H_core[i, j] + Fock[i, j])
            self.energy = self.nuc_energy + self.el_energy

            if np.linalg.norm(self.P - self.P_previous) < tol:
                # calculate dipole moments
                for i in range(3):
                    self.basis.mu[i] = -2 * np.trace(np.dot(self.P, self.basis.M[i])) + \
                                     sum([atom.charge * (atom.coordinates[i] -
                                                         self.mol.center_of_charge[i]) for atom in self.mol.atoms])
                self.basis.mu *= 2.541765

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
