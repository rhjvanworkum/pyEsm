"""
This module was only implemented for faster testing purposes
"""

# from openff.toolkit.typing.engines.smirnoff import ForceField
#
# forcefield = ForceField('openff-1.0.0.offxml')
#
#
# def get_parameter(mol, param, object):
#
#     if param == 'ProperTorsions':
#         isTorsion = True
#     else:
#         isTorsion = False
#
#     if param != "el": smarts = mol.get_smart_string(object, isTorsion)
#
#     match = []
#
#     if param == "el":
#         match = [forcefield.get_parameter_handler("Electrostatics").cutoff.value_in_unit(
#                              forcefield.get_parameter_handler("Electrostatics").cutoff.unit)]
#     else:
#         for parameter in forcefield.get_parameter_handler(param).parameters:
#             if matches(parameter.smirks, smarts):
#                 if param == 'Bonds':
#                     match = [parameter.k.value_in_unit(parameter.k.unit), parameter.length.value_in_unit(parameter.length.unit)]
#                 elif param == 'Angles':
#                     match = [parameter.k.value_in_unit(parameter.k.unit), parameter.angle.value_in_unit(parameter.angle.unit)]
#                 elif param == 'ProperTorsions':
#                     match.append(parameter.k)
#                     match.append(parameter.periodicity)
#                     match.append(parameter.phase)
#                 elif param == 'vdW':
#                     match = [parameter.epsilon.value_in_unit(parameter.epsilon.unit),
#                              parameter.rmin_half.value_in_unit(parameter.rmin_half.unit),
#                              forcefield.get_parameter_handler(param).cutoff.value_in_unit(
#                                  forcefield.get_parameter_handler(param).cutoff.unit)]
#
#     if len(match) != 0:
#         return match
#     else:
#         return "Did not find a match"
#
#
# def matches(smirks, smarts):
#     s_atoms, s_bonds = parse_smirk(smirks)
#     atoms, bonds = parse_smirk(smarts)
#
#     match_atoms = [False for i in range(len(s_atoms))]
#     match_bonds = [False for i in range(len(s_bonds))]
#
#     for i in range(len(s_atoms)):
#         s_a = s_atoms[i]
#         for a in atoms:
#             if s_a == a or s_a == '*':
#                 match_atoms[i] = True
#
#     for i in range(len(s_bonds)):
#         s_b = s_bonds[i]
#         for b in bonds:
#             if s_b == b or s_b == '~':
#                 match_bonds[i] = True
#
#     if False not in match_atoms and False not in match_bonds:
#         return True
#     else:
#         return False
#
#
# def parse_smirk(smirks):
#     atoms = ''
#     bonds = ''
#     type = 'none'
#     for i in range(len(smirks)):
#         if smirks[i] == '[':
#             type = 'atom'
#         if smirks[i] == ']':
#             type = 'bond'
#
#         if type == 'atom' and smirks[i] != '[':
#             if smirks[i] == ':':
#                 atoms += ' '
#                 type = 'none'
#             else:
#                 atoms += smirks[i]
#
#         if type == 'bond' and smirks[i] != ']':
#             bonds += smirks[i] + ' '
#
#     return atoms.split(' ')[:-1], bonds.split(' ')[:-1]