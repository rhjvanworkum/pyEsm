"""
Test Module for the Hartree Fock procedure with a CGTO basis_set set in the py_esm module


"""
import pytest

from py_esm.models.Molecule import Molecule
from py_esm.models.methods.hartree_fock import hf_mol
from py_esm.models.basis_set.CgtoBasisSet import CgtoBasisSet, PreComputedCgtoBasisSet


def test_scf_procedure():
    mol = Molecule('O')
    mol.set_atom_positions(
        [0, 1, 2],
        [
            [0.000000, -0.075791844, 0.000000],
            [0.866811829, 0.601435779, 0.000000],
            [-0.866811829, 0.601435779, 0.000000]
        ]
    )

    base_path = 'C:/Users/rhjva/PycharmProjects/pyESM/tests/data/'

    basis = PreComputedCgtoBasisSet(base_path + 's.dat', base_path + 't.dat', base_path + 'v.dat', base_path + 'eri.dat')

    evals, evecs, energy = hf_mol(mol, basis)

    assert energy == pytest.approx(-74.965901, abs=5e-2)


def test_hartree_fock_procedure():
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

    evals, evecs, energy = hf_mol(mol, basis)

    assert energy == pytest.approx(-74.965901, abs=5e-2)

# TODO: write test for DIIS scf procedure