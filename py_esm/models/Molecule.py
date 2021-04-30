from openbabel.pybel import readstring, ob
# from rdkit import Chem
# from rdkit.Chem import AllChem
# from rdkit.Chem.rdmolfiles import MolFragmentToSmarts
from pysmiles import read_smiles
import numpy as np
from scipy.spatial.transform import Rotation as rot
import functools

from py_esm.helpers.geometry import sort_bend_angle_atoms, sort_torsion_bonds, sort_bend_angle
from py_esm.models.Atom import Atom
from py_esm.utils.data import PySmilesCopy

ob.obErrorLog.StopLogging()

class Molecule:
    """
    Molecule class is a data structure containing 3 types of data:
    1. Cheminformatical: smiles, OpenBabel Mol class, RDkit Mol class
    2. Geometrical: Bonds, angles & torsions
    3. Compositional: Array of Atoms
    """

    def __init__(self, smiles):
        """
        Initalize Molecule Class
        :param smiles: smiles string
        """


        # Cheminformatics section
        self.smiles = smiles

        self.pybel_mol = readstring("smi", smiles)
        self.pybel_mol.make3D()

        self.mol_formula = self.pybel_mol.formula

        # self.rd_mol = Chem.MolFromSmiles(smiles)
        # self.rd_mol = Chem.AddHs(self.rd_mol)
        # AllChem.EmbedMolecule(self.rd_mol)
        # AllChem.MMFFOptimizeMolecule(self.rd_mol)

        try:
            self.pysmiles_mol = read_smiles(smiles, explicit_hydrogen=True)
        except ValueError:
            self.pysmiles_mol = PySmilesCopy([0, 1, 1])


        # compositional section
        self.natoms = len(self.pybel_mol.atoms)
        self.position = np.array([0.0, 0.0, 0.0], dtype=float)
        self.atoms = []

        for i in range(self.natoms):
            atom = self.pybel_mol.atoms[i]
            self.atoms.append(Atom(atom.atomicnum, atom.coords))

        # geometrical section
        self.bonds = []
        self.bond_orders = []
        for bond in self.pysmiles_mol.edges(data='order'):
            self.bonds.append(bond[:-1])
            self.bond_orders.append(bond[2])

        self.angles = []
        for i in range(self.natoms):
            bonds = self.get_bonds(i)

            for j in range(len(bonds) - 1):
                self.angles.append(sort_bend_angle(bonds[j] + bonds[j + 1]))

        self.torsions = []
        for index1 in range(self.natoms):
            bonds = self.get_bonds(index1)
            if len(bonds) >= 2:
                for i, middle_bond in enumerate(bonds):
                    index2 = [atom for atom in middle_bond if atom != index1][0]
                    bonds2 = self.get_bonds(index2)
                    if len(bonds2) >= 2:

                        if i == len(bonds) - 1:
                            bond1 = bonds[0]
                        else:
                            bond1 = bonds[i + 1]

                        for bond in bonds2:
                            if sorted(bond) != sorted(middle_bond):
                                bond2 = bond

                        self.torsions.append(sort_torsion_bonds(bond1, middle_bond, bond2))
            else:
                continue

    """ commented this out because of the rdkit dependecy """
    # @classmethod
    # def from_molfile(cls, file):
    #     """
    #     :param file: asbolute path to the Mol file
    #     :return: Molecule class
    #     """
    #     mol = Chem.MolFromMolFile(file)
    #     smiles = Chem.MolToSmiles(mol)
    #     return cls(smiles)

    def set_atom_positions(self, atoms, positions):
        """
        converts angstrom to bohr as well
        :param atoms: indices
        :param positions: positions
        :return:
        """
        # / 0.52917721092
        for atom in atoms:
            self.atoms[atom].coordinates = np.array([p * 1.889725989 for p in positions[atom]])

    @property
    def n_electrons(self):
        """
        :return: number of electrons in the molecule
        """
        n = 0
        for atom in self.atoms:
            n += atom.atom

        return n

    @property
    def e_nuc_repulsion(self):
        """
        Calculates the classical nuclear repulsion energy
        :return: nuclear repulsion energy
        """
        e_nuc = 0
        for a in self.atoms:
            for b in self.atoms:
                if a != b:
                    e_nuc += (a.charge * b.charge) / np.linalg.norm(a.coordinates - b.coordinates)
        return 0.5 * e_nuc

    @property
    def center_of_charge(self):
        return np.array([sum([atom.charge * atom.coordinates[0] for atom in self.atoms]),
                         sum([atom.charge * atom.coordinates[1] for atom in self.atoms]),
                         sum([atom.charge * atom.coordinates[2] for atom in self.atoms])]) * \
                         (1.0 / sum([atom.charge for atom in self.atoms]))

    @functools.lru_cache(maxsize=None)
    def get_nbonds(self, index):
        """
        :param index: index of atom in the atoms list
        :return: number of bonds of corresponding atom
        """
        n = 0
        for bond in self.bonds:
            if index in bond:
                n += 1

        return n

    @functools.lru_cache(maxsize=None)
    def get_bonds(self, index):
        """
        :param index: index of atom in the atoms list
        :return: array of bonds of corresponding atom
        """
        bonds = []
        for bond in self.bonds:
            if index in bond:
                bonds.append(bond)

        return bonds

    def get_bonds_like(self, bond):
        """
        :param bond: bond type to search for
        :return: array of alike bonds
        """
        bonds = []
        for b in self.bonds:
            atoms = [self.atoms[b[0]].atom, self.atoms[b[1]].atom]
            if bond in (atoms, [atoms[1], atoms[0]]):
                bonds.append(b)

        return bonds

    def get_angles_like(self, angle):
        """
        :param angle: angle type to search for
        :return: array of alike angles
        """
        angles = []
        for a in self.angles:
            if [self.atoms[i].atom for i in a] != angle:
                continue
            else:
                angles.append(a)

        return angles

    def get_bond_length(self, bond):
        """
        :param bond: bond to examine
        :return: distance(bond length) between atom 1 and atom 2
        """
        return np.sqrt(np.sum(np.power(self.atoms[bond[0]].coordinates - self.atoms[bond[1]].coordinates, 2)))

    def norm_dist_vec(self, bond):
        """
        :param bond: bond to examine
        :return: normalized unit vector between atom 1 and atom 2
        """
        return -(self.atoms[bond[0]].coordinates - self.atoms[bond[1]].coordinates) / self.get_bond_length(bond)

    def get_bend_angle(self, angle):
        """
        :param angle: angle to calculate
        :return: Bend angle between atom 12 and atom 23 in degrees
        """
        Eji = self.norm_dist_vec([angle[1], angle[0]])
        Ejk = self.norm_dist_vec([angle[1], angle[2]])

        return np.arccos(Eji.dot(Ejk)) / np.pi * 180

    def get_torsion_angle(self, torsion):
        """
        :param torsion: torsion to examine
        :return: Torsional angle between atom 123 and atom 234 in degrees
        """
        Eij = self.norm_dist_vec([torsion[0], torsion[1]])
        Ejk = self.norm_dist_vec([torsion[1], torsion[2]])
        Ekl = self.norm_dist_vec([torsion[2], torsion[3]])

        tau = (np.cross(Eij, Ejk).dot(np.cross(Ejk, Ekl))) / \
              (np.sin(self.get_bend_angle([torsion[0], torsion[1], torsion[2]])) *
               np.sin(self.get_bend_angle([torsion[1], torsion[2], torsion[3]])))

        while tau > 1 or tau < -1:
            if tau > 1:
                tau -= 2
            else:
                tau += 2

        return np.arccos(tau) / np.pi * 180

    """ ForceFields are no longer being used in this project """
    # def get_smart_string(self, atoms, isTorsion):
    #     """
    #     :param atoms: indices of the atoms
    #     :param isTorsion: Bool indicating a torsional term(4 atoms) or not
    #     :return: smart string of the molecule fragment, specified by the
    #     atoms
    #     """
    #     nbonds = []
    #     for atom in atoms:
    #         nbonds.append(len(self.get_bonds(atom)))
    #
    #     interbonds = MolFragmentToSmarts(self.rd_mol, atoms).replace('(', '').replace(')', '').split(']')
    #
    #     smarts = ''
    #     for i in range(len(atoms)):
    #         if i == len(atoms) - 1:
    #             end = ''
    #         else:
    #             end = interbonds[i + 1][0]
    #         if (nbonds[i] > 1 and self.atoms[atoms[i]].atom == 6) or (self.atoms[atoms[i]].atom > 1 and isTorsion):
    #             smarts += '[#' + str(self.atoms[atoms[i]].atom) + 'X' + str(nbonds[i]) + ':' + str(i + 1) + ']' + end
    #         else:
    #             smarts += '[#' + str(self.atoms[atoms[i]].atom) + ':' + str(i + 1) + ']' + end
    #     return smarts.replace('.', '-')

    def iterate_atoms(self, atom, bonds):
        """
        Iterater function over the molecule, return all the atoms that are on the graph
        of the molecule, except for the ones excluded by the given bonds
        :param atom: atom index of which node to start on the graph
        :param bonds: list of bonds to exclude from the iterator
        :return: list of atoms
        """
        atoms = [atom]
        curr_atoms = [atom]
        prev_bonds = bonds

        looping = True

        while looping:
            new_curr_atoms = []

            for a in curr_atoms:
                if self.get_nbonds(a) > 1:
                    for bond in self.get_bonds(a):
                        if bond not in prev_bonds:
                            # at bond to walked down bonds
                            prev_bonds.append(bond)
                            # add atoms to new current atoms
                            if bond[0] == a:
                                new_curr_atoms.append(bond[1])
                            else:
                                new_curr_atoms.append(bond[0])

            for i in new_curr_atoms:
                atoms.append(i)

            curr_atoms = new_curr_atoms

            if len(curr_atoms) == 0:
                looping = False

        return atoms

    def set_geometry(self, value, dof):
        """
        Sets a certain geometry on the molecule object
        :param value: value the degree of freedom should take
        :param dof: the degree of freedom to set it's value to
        -> when settings an angle enter the value in radians
        """

        if len(dof) > 3:
            print("Currently only DOF's of order 2 or 3 are supported")

        # bond length
        elif len(dof) == 2:
            # gather all affected bonds
            bonds = self.get_bonds_like(dof)
            # calculate needed update
            for b in bonds:
                R = self.get_bond_length(b)
                offset = 0.5 * (value - R)

                for atom in self.iterate_atoms(b[0], [b]):
                    rij = self.norm_dist_vec(b)
                    self.atoms[atom].coordinates -= offset * rij

                for atom in self.iterate_atoms(b[1], [b]):
                    rij = self.norm_dist_vec([b[1], b[0]])
                    self.atoms[atom].coordinates -= offset * rij

        # bend angles
        elif len(dof) == 3:
            # gather all affected angles
            angles = self.get_angles_like(dof)
            # calculate needed updates
            for a in angles:

                # sort the angle array
                a = sort_bend_angle_atoms(a, [self.get_bonds(i) for i in a])

                theta = self.get_bend_angle(a) * np.pi / 180
                diff = (value - theta) / 2

                r1 = (self.atoms[a[0]].coordinates - self.atoms[a[1]].coordinates) / self.get_bond_length([a[0], a[1]])
                r2 = (self.atoms[a[2]].coordinates - self.atoms[a[1]].coordinates) / self.get_bond_length([a[2], a[1]])
                norm_vec = np.cross(r1, r2)
                norm_vec /= np.linalg.norm(norm_vec)

                for atom in self.iterate_atoms(a[0], [(a[0], a[1])]):
                    rotation_vector = - diff * norm_vec
                    rotation = rot.from_rotvec(rotation_vector)
                    self.atoms[atom].coordinates = rotation.apply(self.atoms[atom].coordinates)

                for atom in self.iterate_atoms(a[2], [(a[1], a[2])]):
                    rotation_vector = diff * norm_vec
                    rotation = rot.from_rotvec(rotation_vector)
                    self.atoms[atom].coordinates = rotation.apply(self.atoms[atom].coordinates)