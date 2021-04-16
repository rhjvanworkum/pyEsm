import numpy as np

from py_esm.models.basis_set.CgtoBasisSet import CgtoBasisSet
from py_esm.models.methods.hartree_fock import run_hf


def ccsd_mol(mol):

    basis = CgtoBasisSet(mol, basis_set='sto-3g')
    calculator = CCSDEnergyCalculator(mol.n_electrons, basis.n_basis, basis.met)
    (evals, evecs, energy) = run_hf(mol.n_electrons, basis, mol.e_nuc_repulsion, 25, calculator)

    return evals, evecs, energy

# TODO: Test this class more carefully
class CCSDEnergyCalculator():

    def __init__(self, N, dim, multi_electron_tensor):
        self.N = N
        self.spatial_dim = dim
        self.met = multi_electron_tensor

        # spin spatial orbital dimensions
        self.dim = None

    def calculate(self, e_values):
        return self.ECCSD(e_values)

    def ECCSD(self, e_values):

        # spin dimensionality = 2 * spatial dimensionality, bc 2 electrons per orbital
        self.dim = 2 * self.spatial_dim

        spin_int = np.zeros((self.dim, self.dim, self.dim, self.dim))
        for p in range(0, self.dim):
            for q in range(0, self.dim):
                for r in range(0, self.dim):
                    for s in range(0, self.dim):
                        value1 = self.met[p // 2, r // 2, q // 2, s // 2] * (p % 2 == r % 2) * (q % 2 == s % 2)
                        value2 = self.met[p // 2, s // 2, q // 2, r // 2] * (p % 2 == s % 2) * (q % 2 == r % 2)
                        spin_int[p, q, r, s] = value1 - value2

        # orbital energy values from HF calculation
        fs = np.zeros((self.dim))
        for i in range(self.dim):
            fs[i] = e_values[i // 2]
        fs = np.diag(fs)

        # Init empty T1 (ts) and T2 (td) arrays
        ts = np.zeros((self.dim, self.dim))
        td = np.zeros((self.dim, self.dim, self.dim, self.dim))

        # Initial guess T2 --- from MP2 calculation!
        for a in range(self.N, self.dim):
            for b in range(self.N, self.dim):
                for i in range(0, self.N):
                    for j in range(0, self.N):
                        # this is zero, unless a, b, i, j are all different
                        # this makes sense, because of the orthogonality of the basis functions
                        # indeed if a, b, i, j are not all different (ab|ij) will just be o
                        td[a, b, i, j] += spin_int[i, j, a, b] / (fs[i, i] + fs[j, j] - fs[a, a] - fs[b, b])

        # equations to calculate the Coupled Cluster excitation coeffiecents from Stanton start here

        # Dai - Equation (12) of Stanton
        Dai = np.zeros((self.dim, self.dim))
        for a in range(self.N, self.dim):
            for i in range(0, self.N):
                Dai[a, i] = fs[i, i] - fs[a, a]

        # Dabij - Stanton eq (13)
        Dabij = np.zeros((self.dim, self.dim, self.dim, self.dim))
        for a in range(self.N, self.dim):
            for b in range(self.N, self.dim):
                for i in range(0, self.N):
                    for j in range(0, self.N):
                        Dabij[a, b, i, j] = fs[i, i] + fs[j, j] - fs[a, a] - fs[b, b]

        e_ccsd = 0
        delta_e = 1.0
        while delta_e > 1e-9:
            prev_e_ccsd = e_ccsd
            Fae, Fmi, Fme, Wmnij, Wabef, Wmbej = self.updateintermediates(fs, ts, td, spin_int)
            ts = self.makeT1(fs, ts, td, spin_int, Dai, Fae, Fmi, Fme)
            td = self.makeT2(ts, td, spin_int, Dabij, Fae, Fmi, Fme, Wmnij, Wabef, Wmbej)
            e_ccsd = self.ccsdenergy(fs, ts, td, spin_int)
            delta_e = abs(e_ccsd - prev_e_ccsd)

        return e_ccsd

    # Stanton eq (9)
    def taus(self, a, b, i, j, ts, td):
        taus = td[a, b, i, j] + 0.5 * (ts[a, i] * ts[b, j] - ts[b, i] * ts[a, j])
        return taus

    # Stanton eq (10)
    def tau(self, a, b, i, j, ts, td):
        tau = td[a, b, i, j] + ts[a, i] * ts[b, j] - ts[b, i] * ts[a, j]
        return tau


    # We need to update our intermediates at the beginning, and
    # at the end of each iteration. Each iteration provides a new
    # guess at the amplitudes T1 (ts) and T2 (td), that *hopefully*
    # converges to a stable, ground-state, solution.
    # returns Fae, Fmi, Fme, Wmnij, Wabef, Wmbej
    def updateintermediates(self, fs, ts, td, spin_int):
        # Stanton eq (3)
        Fae = np.zeros((self.dim, self.dim))
        for a in range(self.N, self.dim):
            for e in range(self.N, self.dim):
                Fae[a, e] = (1 - (a == e)) * fs[a, e]
                for m in range(0, self.N):
                    Fae[a, e] += -0.5 * fs[m, e] * ts[a, m]
                    for f in range(self.N, self.dim):
                        Fae[a, e] += ts[f, m] * spin_int[m, a, f, e]
                        for n in range(0, self.N):
                            Fae[a, e] += -0.5 * self.taus(a, f, m, n, ts, td) * spin_int[m, n, e, f]

        # Stanton eq (4)
        Fmi = np.zeros((self.dim, self.dim))
        for m in range(0, self.N):
            for i in range(0, self.N):
                Fmi[m, i] = (1 - (m == i)) * fs[m, i]
                for e in range(self.N, self.dim):
                    Fmi[m, i] += 0.5 * ts[e, i] * fs[m, e]
                    for n in range(0, self.N):
                        Fmi[m, i] += ts[e, n] * spin_int[m, n, i, e]
                        for f in range(self.N, self.dim):
                            Fmi[m, i] += 0.5 * self.taus(e, f, i, n, ts, td) * spin_int[m, n, e, f]

        # Stanton eq (5)
        Fme = np.zeros((self.dim, self.dim))
        for m in range(0, self.N):
            for e in range(self.N, self.dim):
                Fme[m, e] = fs[m, e]
                for n in range(0, self.N):
                    for f in range(self.N, self.dim):
                        Fme[m, e] += ts[f, n] * spin_int[m, n, e, f]

        # Stanton eq (6)
        Wmnij = np.zeros((self.dim, self.dim, self.dim, self.dim))
        for m in range(0, self.N):
            for n in range(0, self.N):
                for i in range(0, self.N):
                    for j in range(0, self.N):
                        Wmnij[m, n, i, j] = spin_int[m, n, i, j]
                        for e in range(self.N, self.dim):
                            Wmnij[m, n, i, j] += ts[e, j] * spin_int[m, n, i, e] - ts[e, i] * spin_int[m, n, j, e]
                            for f in range(self.N, self.dim):
                                Wmnij[m, n, i, j] += 0.25 * self.tau(e, f, i, j, ts, td) * spin_int[m, n, e, f]

        # Stanton eq (7)
        Wabef = np.zeros((self.dim, self.dim, self.dim, self.dim))
        for a in range(self.N, self.dim):
            for b in range(self.N, self.dim):
                for e in range(self.N, self.dim):
                    for f in range(self.N, self.dim):
                        Wabef[a, b, e, f] = spin_int[a, b, e, f]
                        for m in range(0, self.N):
                            Wabef[a, b, e, f] += ts[a, m] * spin_int[b, m, e, f] - ts[b, m] * spin_int[a, m, e, f]
                            for n in range(0, self.N):
                                Wabef[a, b, e, f] += 0.25 * self.tau(a, b, m, n, ts, td) * spin_int[m, n, e, f]

        # Stanton eq (8)
        Wmbej = np.zeros((self.dim, self.dim, self.dim, self.dim))
        for m in range(0, self.N):
            for b in range(self.N, self.dim):
                for e in range(self.N, self.dim):
                    for j in range(0, self.N):
                        Wmbej[m, b, e, j] = spin_int[m, b, e, j]
                        for f in range(self.N, self.dim):
                            Wmbej[m, b, e, j] += ts[f, j] * spin_int[m, b, e, f]
                        for n in range(0, self.N):
                            Wmbej[m, b, e, j] += -ts[b, n] * spin_int[m, n, e, j]
                            for f in range(self.N, self.dim):
                                Wmbej[m, b, e, j] += -(0.5 * td[f, b, j, n] + ts[f, j] * ts[b, n]) * spin_int[
                                    m, n, e, f]

        return Fae, Fmi, Fme, Wmnij, Wabef, Wmbej


    # Stanton eq (1)
    def makeT1(self, fs, ts, td, spin_int, Dai, Fae, Fmi, Fme):
        tsnew = np.zeros((self.dim, self.dim))
        for i in range(0, self.N):
            for a in range(self.N, self.dim):
                tsnew[a, i] = fs[i, a]
                for e in range(self.N, self.dim):
                    tsnew[a, i] += ts[e, i] * Fae[a, e]
                for m in range(0, self.N):
                    tsnew[a, i] += -ts[a, m] * Fmi[m, i]
                    for e in range(self.N, self.dim):
                        tsnew[a, i] += td[a, e, i, m] * Fme[m, e]
                        for f in range(self.N, self.dim):
                            tsnew[a, i] += -0.5 * td[e, f, i, m] * spin_int[m, a, e, f]
                        for n in range(0, self.N):
                            tsnew[a, i] += -0.5 * td[a, e, m, n] * spin_int[n, m, e, i]
                for n in range(0, self.N):
                    for f in range(self.N, self.dim):
                        tsnew[a, i] += -ts[f, n] * spin_int[n, a, i, f]
                tsnew[a, i] = tsnew[a, i] / Dai[a, i]
        return tsnew

    # Stanton eq (2)
    def makeT2(self, ts, td, spin_int, Dabij, Fae, Fmi, Fme, Wmnij, Wabef, Wmbej):
        tdnew = np.zeros((self.dim, self.dim, self.dim, self.dim))
        for a in range(self.N, self.dim):
            for b in range(self.N, self.dim):
                for i in range(0, self.N):
                    for j in range(0, self.N):
                        tdnew[a, b, i, j] += spin_int[i, j, a, b]
                        for e in range(self.N, self.dim):
                            tdnew[a, b, i, j] += td[a, e, i, j] * Fae[b, e] - td[b, e, i, j] * Fae[a, e]
                            for m in range(0, self.N):
                                tdnew[a, b, i, j] += -0.5 * td[a, e, i, j] * ts[b, m] * Fme[m, e] + 0.5 * td[
                                    b, e, i, j] * ts[a, m] * Fme[m, e]
                                continue
                        for m in range(0, self.N):
                            tdnew[a, b, i, j] += td[a, b, j, m] * Fmi[m, i] - td[a, b, i, m] * Fmi[m, j]
                            for e in range(self.N, self.dim):
                                tdnew[a, b, i, j] += -0.5 * td[a, b, i, m] * ts[e, j] * Fme[m, e] + 0.5 * td[
                                    a, b, j, m] * ts[e, i] * Fme[m, e]
                                continue
                        for e in range(self.N, self.dim):
                            tdnew[a, b, i, j] += ts[e, i] * spin_int[a, b, e, j] - ts[e, j] * spin_int[a, b, e, i]
                            for f in range(self.N, self.dim):
                                tdnew[a, b, i, j] += 0.5 * self.tau(e, f, i, j, ts, td) * Wabef[a, b, e, f]
                                continue
                        for m in range(0, self.N):
                            tdnew[a, b, i, j] += -ts[a, m] * spin_int[m, b, i, j] + ts[b, m] * spin_int[m, a, i, j]
                            for e in range(self.N, self.dim):
                                tdnew[a, b, i, j] += td[a, e, i, m] * Wmbej[m, b, e, j] - ts[e, i] * ts[a, m] * \
                                                     spin_int[m, b, e, j]
                                tdnew[a, b, i, j] += -td[a, e, j, m] * Wmbej[m, b, e, i] + ts[e, j] * ts[a, m] * \
                                                     spin_int[m, b, e, i]
                                tdnew[a, b, i, j] += -td[b, e, i, m] * Wmbej[m, a, e, j] + ts[e, i] * ts[b, m] * \
                                                     spin_int[m, a, e, j]
                                tdnew[a, b, i, j] += td[b, e, j, m] * Wmbej[m, a, e, i] - ts[e, j] * ts[b, m] * \
                                                     spin_int[m, a, e, i]
                                continue
                            for n in range(0, self.N):
                                tdnew[a, b, i, j] += 0.5 * self.tau(a, b, m, n, ts, td) * Wmnij[m, n, i, j]
                                continue
                        tdnew[a, b, i, j] = tdnew[a, b, i, j] / Dabij[a, b, i, j]
        return tdnew

    def ccsdenergy(self, fs, ts, td, spin_int):
        ECCSD = 0.0
        for i in range(0, self.N):
            for a in range(self.N, self.dim):
                ECCSD += fs[i, a] * ts[a, i]
                for j in range(0, self.N):
                    for b in range(self.N, self.dim):
                        ECCSD += 0.25 * spin_int[i, j, a, b] * td[a, b, i, j] + 0.5 * spin_int[i, j, a, b] * (
                        ts[a, i]) * (ts[b, j])
        return ECCSD