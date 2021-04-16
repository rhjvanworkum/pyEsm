"""
Molecule Test module for the py_esm module

Data points for the geometrical tests are extracted for the acetaldehyde molecule using ChemDoodle 3D
"""
import pytest
import numpy as np

from py_esm.models.Molecule import Molecule


@pytest.mark.parametrize("bond, bond_length", [
    ([0, 1], 1.07),
    ([1, 3], 1.27),
    ([1, 2], 1.53),
    ([2, 4], 1.07),
    ([2, 5], 1.07),
    ([2, 6], 1.07),
])
def test_molecule_bonds(bond, bond_length):
    m = Molecule('[H]C(C)=O')
    assert m.get_bond_length(bond) == pytest.approx(bond_length, abs=5e-2)


@pytest.mark.parametrize("angle, bend_angle", [
    ([0, 1, 2], 120.46),
    ([2, 1, 3], 119.97),
    ([1, 2, 4], 109.58),
    ([4, 2, 5], 109.47),
    ([5, 2, 6], 109.47),
])
def test_molecule_angles(angle, bend_angle):
    m = Molecule('[H]C(C)=O')
    assert m.get_bend_angle(angle) == pytest.approx(bend_angle, abs=5e-1)


@pytest.mark.parametrize("torsion, torsional_angle", [
    ([3, 1, 2, 6], 120),
    ([4, 2, 1, 3], 10),
])
def test_molecule_torsions(torsion, torsional_angle):
    m = Molecule('[H]C(C)=O')
    assert m.get_torsion_angle(torsion) == pytest.approx(torsional_angle, abs=1)


@pytest.mark.parametrize("value, dof", [
    (1, [6, 6]),
    (3, [6, 8]),
    (5, [6, 1])
])
def test_molecule_set_geometry_bonds(value, dof):
    m = Molecule('[H]C(C)=O')
    m.set_geometry(value, dof)

    bonds = m.get_bonds_like(dof)
    for bond in bonds:
        assert m.get_bond_length(bond) == pytest.approx(value, abs=1e-1)


# TODO: Fix the accuracy of setting the angles in the molecular geometry
@pytest.mark.parametrize("value, dof", [
    (110 * (np.pi / 180), [6, 6, 1]),
    (120 * (np.pi / 180), [6, 6, 8]),
    (110 * (np.pi / 180), [1, 6, 1])
])
def test_molecule_set_geometry_angles(value, dof):
    m = Molecule('CC(=O)C')
    m.set_geometry(value, dof)

    angles = m.get_angles_like(dof)
    for angle in angles:
        a = m.get_bend_angle(angle) * (np.pi / 180)
        assert a == pytest.approx(value, abs=2e-1)
