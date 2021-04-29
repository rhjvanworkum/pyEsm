import numpy as np

from py_esm.models.basis_set.CgtoBasisSet import CgtoBasisSet


class PESConstructor:

    def __init__(self, molecule, basis, axis, dofs, method, dimensions, num_dim):
        self.m = molecule
        self.basis = basis
        self.axis = axis
        self.dofs = dofs
        self.method = method
        self.output_array = np.zeros(dimensions).flatten()
        self.dimensions = dimensions
        self.num_dim = num_dim
        self.indices = np.zeros(num_dim, dtype=int)

    """ Sets element in the output array at the current index of the class """
    def set_value_in_output_array(self, value):
        index = np.ravel_multi_index(self.indices, self.dimensions)
        self.output_array[index] = value

    """ Get Element from an array ? """
    def get_element_from_array(self, array):
        element = 0

        for i in self.indices:
            element = array[i]

        return element

    """ loop through all the coordinates in all the DOF's """
    def loop_through_coordinates(self, n):
        if n >= 1:
            for x in range(self.dimensions[n - 1]):
                self.indices[n - 1] = x
                self.loop_through_coordinates(n - 1)
        else:
            for index, value in enumerate(self.indices):
                coordinate = self.axis[index][value]
                self.m.set_geometry(coordinate, self.dofs[index])

            # basis = CgtoBasisSet(self.m, self.basis)
            # basis
            # E, C,
            energy = self.method(self.m)

            self.set_value_in_output_array(energy)

    """ returns the PES data in the right dimensionality """
    @property
    def data(self):
        return self.output_array.reshape(self.dimensions)

    """ returns a certain slice of the PES data """
    def get_slice(self, indices, values):
        axis = self.output_array.reshape(self.dimensions)

        for idx, val in enumerate(values):
            if idx not in indices:
                axis = axis[val]
            else:
                axis = axis[:]

        return axis

    """ returns all gradient equals zero points """
    @property
    def extrema(self):
        return None

    """ returns all the curvature equals zero points """
    @property
    def transitions(self):
        return None