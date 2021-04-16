element_dict = {
    'H':  1,
    'He': 2,
    'Li': 3,
    'Be': 4,
    'B':  5,
    'C':  6,
    'N':  7,
    'O':  8,
    'F':  9,
    'Ne': 10,
    'Na': 11,
    'Mg': 12,
    'Al': 13,
    'Si': 14,
    'P':  15,
    'S':  16,
    'Cl': 17,
    'Ar': 18
}

zeta_dict = {
    'H':  1.24,
    'He': 2.0925,
    'Li': 3.00,
    'Be': 4.00,
    'B':  5.00,
    'C':  6.00,
    'N':  7.00,
    'O':  8.00,
    'F':  9.00,
    'Ne': 10.00,
    'Na': 11.00,
    'Mg': 12.00,
    'Al': 13.00,
    'Si': 14.00,
    'P':  15.00,
    'S':  16.00,
    'Cl': 17.00,
    'Ar': 18.00
}

max_quantum_number = {
    'H':  1,
    'He': 1,
    'Li': 2,
    'Be': 2,
    'B':  3,
    'C':  3,
    'N':  3,
    'O':  3,
    'F':  3,
    'Ne': 3,
    'Na': 4,
    'Mg': 4,
    'Al': 5,
    'Si': 5,
    'P':  5,
    'S':  5,
    'Cl': 5,
    'Ar': 5
}

def get_element_symbol(atom : int):
    """
    Retrieves Element name , based on atom number
    :param atom: atom number
    :return: atom symbol string
    """
    values = list(element_dict.values())
    for idx, v in enumerate(values):
        if v == atom:
            return list(element_dict.keys())[idx]

def get_zeta_value(atom: int):
    """
    Retrieves zeta value, based on atom number
    :param atom: atom number
    :return: zeta value
    """
    return zeta_dict[get_element_symbol(atom)]
