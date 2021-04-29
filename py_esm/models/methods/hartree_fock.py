from py_esm.models.methods.scf.HartreeFockProcedure import HartreeFockProcedure


def hf_mol(mol, basis):
    """
    Hartree Fock calculation on a molecule
    :param mol: Molecule Object
    :return: eigenvalues, eigenvectors and total energy of the system
    """

    scf = HartreeFockProcedure(mol, basis)
    scf.run_hf(1e-12, DIIS=True)

    return scf.E, scf.C, scf.energy


def hf_mol_iterated(mol, basis):
    """
    Hartree Fock calculation on a molecule
    :param mol: Molecule Object
    :return: eigenvalues, eigenvectors and total energy of the system
    """

    scf = HartreeFockProcedure(mol, basis)
    scf.run_hf(1e-12, DIIS=True, report_iterations=True)

    return scf.MO_list, scf.C_list