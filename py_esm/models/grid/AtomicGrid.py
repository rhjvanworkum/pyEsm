import numpy as np

from py_esm.helpers.grid import Bragg
from py_esm.utils.lebedev import lebedev
from py_esm.utils.legendre import leggauss
from py_esm.utils.linalg import get_magnitude


class AtomicGrid:

    def __init__(self, atom):

        grid = LegendreGrid(atom)

        self.points = []

        for r_rad, w_rad, n_ang_points in grid:
            for x_ang, y_ang, z_ang, w_ang in lebedev[n_ang_points]:
                w = w_rad * w_ang
                self.points.append((r_rad * x_ang + atom.coordinates[0],
                                    r_rad * y_ang + atom.coordinates[1],
                                    r_rad * z_ang + atom.coordinates[2],
                                    w))

        self.points = np.array(self.points, dtype=float)
        self.n_points = self.points.shape[0]


def LegendreGrid(atom):
    R_max = 0.5 * Bragg[atom.atom] * 1.8897259886  # Bragg is in angstrom, we want bohr * 1.8897259886

    n_radial = 32
    n_ang = 1

    radial_grid = leggauss(n_radial)

    grid = []

    # becke transformation
    for i in range(n_radial):
        x_rad, w_rad = radial_grid[i]
        r_rad = BeckeTransform(x_rad, get_magnitude(atom.coordinates), R_max)
        dr = 2 * R_max / pow(1 - x_rad, 2)
        vol = 4 * np.pi * r_rad * r_rad * dr
        n_ang_points = ang_mesh(float(i + 1) / n_radial, n_ang)
        grid.append((r_rad, w_rad * vol, n_ang_points))

    return grid


def BeckeTransform(x, r0, R_max):
    return r0 + R_max * (1.0 + x) / (1.0 - x)


def ang_mesh(fraction, fineness):
    ang_levels = [
        [6, 14, 26, 26, 14],
        [50, 50, 110, 50, 26],
        [50, 110, 194, 110, 50],
        [194, 194, 194, 194, 194]
    ]
    alevs = ang_levels[fineness]

    nang = alevs[0]
    if fraction > 0.4: nang = alevs[1]
    if fraction > 0.5: nang = alevs[2]
    if fraction > 0.7: nang = alevs[3]
    if fraction > 0.8: nang = alevs[4]
    return nang
