"""
Test Module for the Post Hartree Fock methods with a CGTO basis_set set in the py_esm module


"""
import pytest

from py_esm.models.Molecule import Molecule
from py_esm.models.methods.ccsd import ccsd_mol
from py_esm.models.methods.mp2 import mp2_mol


def test_mp2():
    mol = Molecule('O')
    mol.set_atom_positions(
        [0, 1, 2],
        [
            [0.000000, -0.075791844, 0.000000],
            [0.866811829, 0.601435779, 0.000000],
            [-0.866811829, 0.601435779, 0.000000]
        ]
    )

    evals, evecs, energy = mp2_mol(mol)

    assert energy == pytest.approx(-75.006040, abs=5e-2)


def test_ccsd():
    mol = Molecule('[HH]')
    mol.set_atom_positions(
        [0, 1],
        [
            [0.000000, 0.000000, 0.373],
            [0.000000, 0.000000, -0.373],
        ]
    )

    evals, evecs, energy = ccsd_mol(mol)

    assert energy == pytest.approx(-1.151522, abs=5e-2)