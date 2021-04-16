from py_esm.utils import elements, linalg


# pylint: disable=too-few-public-methods
class Atom:
    """
    Atom class, representing a single atom
    """

    def __init__(self, atom, coordinates):
        """
        Intializer
        :param atom: Atom Number, e.g. H -> 1
        :param coordinates: list of atom coordinates, xyz format
        """

        self.atom = atom
        self.atom_symbol = elements.get_element_symbol(atom)

        self.coordinates = linalg.list_to_array(coordinates)
        self.charge = atom
