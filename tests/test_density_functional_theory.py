import pytest

from py_esm.models.Molecule import Molecule
from py_esm.models.basis_set.CgtoBasisSet import CgtoBasisSet
from py_esm.models.methods.dft import dft_mol


def test_dft_procedure():
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

    evals, evecs, energy = dft_mol(mol, basis)

    assert energy == pytest.approx(-74.939072, abs=1e-1)