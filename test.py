from py_esm.models.Molecule import Molecule
from py_esm.models.basis_set.CgtoBasisSet import CgtoBasisSet
from py_esm.models.methods.ccsd import ccsd_mol
from py_esm.models.methods.mp2 import mp2_mol

mol = Molecule('[HH]')
mol.set_atom_positions(
    [0, 1],
    [
        [0.000000, 0.000000, 0.373],
        [0.000000, 0.000000, -0.373],
    ]
)

evals, evecs, energy = ccsd_mol(mol)

print(energy)