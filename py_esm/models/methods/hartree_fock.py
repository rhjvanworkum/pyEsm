import numpy as np
import scipy

from py_esm.models.basis_set.CgtoBasisSet import CgtoBasisSet


def hf_mol(mol):
    """
    Hartree Fock calculation on a molecule
    :param mol: Molecule Object
    :return: eigenvalues, eigenvectors and total energy of the system
    """

    basis = CgtoBasisSet(mol, basis_set='sto-3g')
    (evals, evecs, energy) = run_hf(mol.n_electrons, basis, mol.e_nuc_repulsion, 25)

    return evals, evecs, energy


def run_hf(N, basis, e_nuc, n_scf, calculator=None):
    """
    Runs the Hartree Fock SCF procedure
    :param N: number of electrons in the system
    :param basis: the basis_set set object used
    :param e_nuc: the Nuclear Repulsion energy of the system
    :param n_scf: number of iterative cycles to run
    :param calculator: an optional post-HF method calculator object
    :return:
    """
    # intialize variables
    evals = np.zeros(basis.n_basis)
    evecs = np.zeros((basis.n_basis, basis.n_basis))
    e_tot = 0

    B = basis.n_basis
    # T and V are not affected by the iterative process
    H_core = basis.kinetic + basis.extern

    P = np.zeros((B, B), dtype=np.float64)

    for n in range(n_scf):

        Fock = np.zeros((B, B), dtype=np.float64)
        for i in range(B):
            for j in range(B):
                Fock[i, j] = H_core[i, j]
                for x in range(B):
                    for y in range(B):
                        Fock[i, j] += P[x, y] * (2 * basis.coulomb([i, j, x, y]) - basis.exchange([i, x, j, y]))

        evals, evecs = scipy.linalg.eigh(Fock, basis.overlap)

        P = np.zeros((B, B))
        for i in range(B):
            for j in range(B):
                for a in range(int(N / 2)):
                    P[i, j] += evecs[i, a] * evecs[j, a]

        e_tot = e_nuc

        for i in range(B):
            for j in range(B):
                e_tot += P[i, j] * (H_core[i, j] + Fock[i, j])

        if calculator:
            e_tot += calculator.calculate(evals)

    return evals, evecs, e_tot
