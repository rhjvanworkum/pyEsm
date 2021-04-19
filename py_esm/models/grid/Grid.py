import numpy as np

from py_esm.models.grid.AtomicGrid import AtomicGrid
from py_esm.helpers.grid import Bragg


class Grid:

    def __init__(self, atoms):
        atomic_grids = [AtomicGrid(atom) for atom in atoms]

        # todo: becke reweight?
        self.points = np.vstack([grid.points for grid in atomic_grids])
        self.n_points = self.points.shape[0]
        self.funcs = None

    def set_func(self, funcs):
        self.funcs = np.zeros((self.n_points, len(funcs)))
        for j, func in enumerate(funcs):
            for i, (x, y, z, w) in enumerate(self.points):
                self.funcs[i, j] = func(x, y, z)

    def get_rho(self, P):
        return 2 * np.einsum('ij,ik,jk->i', self.funcs, self.funcs, P)
