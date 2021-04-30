import numpy as np
import requests
import json

from py_esm.utils.data import read_2d_matrix, read_4d_matrix
from py_esm.utils.linalg import get_4d_matrix_index
from py_esm.utils.orbitals import ang_mom_dict
import py_esm.models.basis_set.Cgto as cgto


class CgtoBasisSet:

    def __init__(self, mol, basis_set='sto-3g', grid=None):

        elements = sorted(list(set([atom.atom for atom in mol.atoms])))
        r = requests.get('http://basissetexchange.org/api/basis/' + basis_set + '/format/JSON',
                         params={'elements': elements},
                         headers={'User-Agent': 'BSE Example Python Script', 'From': 'bse@molssi.org'})

        if r.status_code != 200:
            print(r.text)
            raise RuntimeError("Could not obtain data from the BSE. Check the error information above")

        response = json.loads(r.text)

        self.basis_functions = []

        for index, atom in enumerate(mol.atoms):
            for n, shell in enumerate(response['elements'][str(atom.atom)]['electron_shells']):
                for l in range(n + 1):
                    m_l = np.arange(-l, l + 1)

                    for ang_n in m_l:
                        self.basis_functions.append(cgto.ContractedGaussianOrbital(
                            origin=atom.coordinates,
                            shell=n,
                            ang_mom_numbers=ang_mom_dict[l][ang_n],
                            exponents=np.array(shell['exponents'], dtype=float),
                            coeffs=np.array(shell['coefficients'][l], dtype=float)
                        ))

        self.B = len(self.basis_functions)

        self.S = np.zeros((self.B, self.B), dtype=np.float64)
        self.T = np.zeros((self.B, self.B), dtype=np.float64)
        self.V = np.zeros((self.B, self.B), dtype=np.float64)
        self.met = np.zeros((self.B, self.B, self.B, self.B), dtype=np.float64)

        self.M = np.zeros((3, self.B, self.B), dtype=np.float64)
        self.mu = np.zeros(3)

        count = 0
        for i in range(self.B):
            for j in range(i + 1):
                self.S[i, j] = self.S[j, i] = cgto.S(self.basis_functions[i], self.basis_functions[j])
                self.T[i, j] = self.T[j, i] = cgto.T(self.basis_functions[i], self.basis_functions[j])

                self.M[0, i, j] = self.M[0, j, i] = cgto.Mu(self.basis_functions[i], self.basis_functions[j],
                                                            mol.center_of_charge, 'x')
                self.M[1, i, j] = self.M[1, j, i] = cgto.Mu(self.basis_functions[i], self.basis_functions[j],
                                                            mol.center_of_charge, 'y')
                self.M[2, i, j] = self.M[2, j, i] = cgto.Mu(self.basis_functions[i], self.basis_functions[j],
                                                            mol.center_of_charge, 'z')

                for atom in mol.atoms:
                    self.V[i, j] += -atom.charge * cgto.V(self.basis_functions[i], self.basis_functions[j], atom.coordinates)
                self.V[j, i] = self.V[i, j]

                for k in range(self.B):
                    for l in range(k + 1):
                        count += 1
                        self.met[i, j, k, l] = self.met[i, j, l, k] = self.met[j, i, k, l] = self.met[j, i, l, k] =\
                            self.met[k, l, i, j] = self.met[l, k, i, j] = self.met[k, l, j, i] = self.met[l, k, j, i]=\
                            cgto.ERI(self.basis_functions[i], self.basis_functions[j],
                                     self.basis_functions[k], self.basis_functions[l])
                        if count % 20 == 0: print('progress: ', (count / self.B ** 4 * 100), '%')

    @property
    def n_basis(self):
        return self.B

    @property
    def overlap(self):
        return self.S

    @property
    def kinetic(self):
        return self.T

    @property
    def extern(self):
        return self.V

    def coulomb(self, indices):
        return self.met[indices[0], indices[1], indices[2], indices[3]]

    def exchange(self, indices):
        return self.met[indices[0], indices[1], indices[2], indices[3]]


class PreComputedCgtoBasisSet:

    def __init__(self, s_path, t_path, v_path, eri_path):
        self.S = read_2d_matrix(s_path)
        self.T = read_2d_matrix(t_path)
        self.V = read_2d_matrix(v_path)
        self.met = read_4d_matrix(eri_path)

        self.B = self.S.shape[0]

        self.M = np.zeros((3, self.B, self.B), dtype=np.float64)
        self.mu = np.zeros(3)

    @property
    def n_basis(self):
        return self.B

    @property
    def overlap(self):
        return self.S

    @property
    def kinetic(self):
        return self.T

    @property
    def extern(self):
        return self.V

    def coulomb(self, indices):
        return self.met[get_4d_matrix_index(indices[0], indices[1], indices[2], indices[3])]

    def exchange(self, indices):
        return self.met[get_4d_matrix_index(indices[0], indices[1], indices[2], indices[3])]
