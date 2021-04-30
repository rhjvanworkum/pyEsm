"""
This File is used to serve the pyEsm module for the front-end on CodingChemistry.com
"""

from flask import Flask, jsonify, request
from flask_cors import CORS, cross_origin
from itertools import permutations
import math
import numpy as np
import json

from py_esm.models.Molecule import Molecule
from py_esm.models.PES import PESConstructor
from py_esm.models.basis_set.CgtoBasisSet import CgtoBasisSet
from py_esm.models.methods.hartree_fock import hf_mol, hf_mol_iterated
from py_esm.models.methods.mp2 import mp2_mol
from py_esm.models.methods.ccsd import ccsd_mol
from py_esm.models.methods.dft import dft_mol

app = Flask(__name__)
CORS(app)
app.config['CORS_HEADERS'] = 'Content-Type'

def unique(list, len):
    dict = []
    npr = math.factorial(len)

    for item in list:
        i = 0
        for c in permutations(item, len):
            if c in dict:
                break
            i += 1

        if (i == npr):
            dict.append(item)

    return dict

def read_url_arg(arg):
    args = []

    if len(arg) == 0:
        args = []
    elif len(arg) == 1:
        args = [int(arg)]
    else:
        args = [int(i) for i in arg.split(',')]

    return args


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


methods_dict = {
    'Hartree-Fock': hf_mol,
    'MP2': mp2_mol,
    'CCSD': ccsd_mol,
    'DFT-LSDA': dft_mol
}

global m
global pes


@app.route('/api/set_molecule', methods=['GET'])
@cross_origin()
def set_molecule():
    global m
    m = Molecule(request.args['smiles'])

    bonds = unique([(m.atoms[b[0]].atom, m.atoms[b[1]].atom) for b in m.bonds], 2)
    angles = unique([(m.atoms[a[0]].atom, m.atoms[a[1]].atom, m.atoms[a[2]].atom) for a in m.angles], 3)
    dof = bonds + angles

    message = {
        'atoms': [[a.atom, a.coordinates.tolist()] for a in m.atoms],
        'bonds': m.bonds,
        'bond_orders': m.bond_orders,
        'dofs': dof
    }

    response = jsonify(message)
    response.headers.add("Access-Control-Allow-Origin", "*")
    return response


@app.route('/api/set_geometry', methods=['POST'])
@cross_origin()
def set_geometry():
    global m

    data = request.json

    m.set_geometry(data["value"], data["dof"])

    print(' completed geometry request ')

    message = {
        'atoms': [[a.atom, a.coordinates.tolist()] for a in m.atoms],
    }

    response = jsonify(message)
    response.headers.add("Access-Control-Allow-Origin", "*")
    return response


@app.route('/api/view_orbitals', methods=['GET'])
@cross_origin()
def view_orbitals():
    global m

    basis_set = request.args['basis']

    basis = CgtoBasisSet(m, basis_set)
    MO, C = hf_mol_iterated(m, basis)

    bfss = []

    bf_pos_x = []
    bf_pos_y = []
    bf_pos_z = []

    for bfs in basis.basis_functions:
        bf_pos_x.append(bfs.origin[0])
        bf_pos_y.append(bfs.origin[1])
        bf_pos_z.append(bfs.origin[2])
        bfss.append([bfs.origin.tolist(), bfs.ang_mom, bfs.coeffs.tolist(), bfs.exponents.tolist()])

    message = {
        'E': [mo.tolist() for mo in MO],
        'C': [c.tolist() for c in C],
        'basis_functions': bfss,
        'box': [[np.min(bf_pos_x), np.min(bf_pos_y), np.min(bf_pos_z)],
                [np.max(bf_pos_x), np.max(bf_pos_y), np.max(bf_pos_z)]]
    }

    response = jsonify(message)
    response.headers.add("Access-Control-Allow-Origin", "*")
    return response


@app.route('/api/run_job', methods=['POST'])
@cross_origin()
def run_job():
    global m, pes, done

    done = False
    data = request.json

    num_axis = len(data['axis'])

    length = []
    axis = []
    dofs = []

    for i in range(num_axis):
        length.append(int((data['axis'][i]['max'] - data['axis'][i]['min']) / data['axis'][i]['interval'] + 1))
        axis.append(np.linspace(data['axis'][i]['min'], data['axis'][i]['max'], length[i]))
        dofs.append(data['axis'][i]['dof'])

    basis = data['basis']
    method = methods_dict[data['method']]

    pes = PESConstructor(m, basis, axis, dofs, method, length, num_axis)
    pes.loop_through_coordinates(num_axis)

    response = jsonify({'message': 'succes'})
    response.headers.add("Access-Control-Allow-Origin", "*")
    return response


@app.route('/api/get_slice', methods=['GET'])
@cross_origin()
def get_slice():
    global pes

    indices = read_url_arg(request.args['indices'])
    values = read_url_arg(request.args['values'])

    data = pes.get_slice(indices, values)

    return json.dumps({'output': data}, cls=NumpyEncoder)


if __name__ == '__main__':

    app.run(port=5000, debug=True, use_reloader=False)