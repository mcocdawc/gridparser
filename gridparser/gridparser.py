import pandas as pd
import chemcoord as cc
import numpy as np
import re
import io
# from . import _pandas_wrapper
from . import export


def _split(string, seperate_newline = True):
    if seperate_newline:
        return list(filter(None, re.split("[, =]+|(\\n)+", string)))
    else:
        return list(filter(None, re.split("[, =]+", string)))


@export
def parse_grid(file):
    """Parse an ASCII formatted MOLCAS grid file.

    Args:
        file (str): Filename, or path to file.

    Returns:
        dict: Dictionary with 4 keys:

            **grid_metadata**: Metadata of the grid.

            **orbitals**: Nested dictionary::

                orbitals[symmetry_charakter][n-th]

            returns the n-th orbital within a given symmetry_charakter.
            The symmetry_charakter is an integer as defined by MOLCAS
            depending on the used symmetry group in the grid file.
            Read more in the MOLCAS manual.
            ``orbitals[1][0]`` is a special case that returns the overall
            density. (Note that the density is totally symmetric so this
            assignment is even physical).

            **orbitals_metadata**: Nested dictionary::

                orbitals_metadata[symmetry_charakter][n-th]

            returns a dictionary for the n-th orbital within a given
            symmetry_charakter with the entries ``energy, occupation, status``.

            **molecule**: A chemcoord instance containing information about the
            coordinates of the molecule.
    """
    grid_metadata = {}
    orbitals_metadata = {}

    f = open(file, 'tr')
    for _ in range(3):
        line = _split(f.readline())

    grid_metadata['Natom'] = int(line[1])

    molecule_in = []
    for i in range(grid_metadata['Natom']):
        line = _split(f.readline())
        # The following removes numbers after element symbol
        element_symbol = re.search("[a-zA-Z]", line[0]).group()
        line[0] = element_symbol
        molecule_in.append(' '.join(line))

    molecule_in = ' '.join(molecule_in)
    molecule_in = str(grid_metadata['Natom']) + 2*'\n' + molecule_in
    molecule_in = io.StringIO(molecule_in)
    molecule = cc.xyz_functions.read_xyz(molecule_in)

    def get_value_and_correct_type(line):
        key = line[0]

        def get_string(line):
            return line[1]

        def get_integer(line):
            return int(line[1])

        def get_boolean(line):
            return bool(line[1])

        def get_floating(line):
            return float(line[1])

        def get_int_array(line):
            return np.array([int(number) for number in line[1:-1]])

        def get_float_array(line):
            return np.array([float(number) for number in line[1:-1]])

        def get_list(line):
            return line[1:-1]

        actions = {}
        actions['VERSION'] = get_string
        actions['N_of_MO'] = get_integer
        actions['N_of_Grids'] = get_integer
        actions['N_of_Points'] = get_integer
        actions['Block_Size'] = get_integer
        actions['N_Blocks'] = get_integer
        actions['Is_cutoff'] = get_boolean
        actions['CutOff'] = get_floating
        actions['N_P'] =  get_integer
        actions['N_INDEX'] = get_list
        actions['Net'] = get_int_array
        actions['Origin'] = get_float_array
        actions['Axis_1'] = get_float_array
        actions['Axis_2'] = get_float_array
        actions['Axis_3'] = get_float_array
        actions['GridName'] =  get_list
        actions['GridName'] = get_string
        # actions['Title'] = get_string  == is apparently not a part of metadata
        return actions[key](line)

    end_of_grid_metadata_reached = False
    while not end_of_grid_metadata_reached:
        line = _split(f.readline())
        if line[0] == 'GridName':
            end_of_grid_metadata_reached = True
        else:
            current_line = f.tell()
            key = line[0]
            grid_metadata[key] = get_value_and_correct_type(line)
    f.seek(current_line, 0)

    orbitals = {}
    order_of_orbitals = []

    for _ in range(grid_metadata['N_of_Grids']):
        line = _split(f.readline())
        try:
            symmetry_charakter = int(line[1])
            number_of_order = int(line[2])
        except ValueError:
            # because int('Density') just does not work
            symmetry_charakter = 1
            number_of_order = 0

        order_of_orbitals.append((symmetry_charakter, number_of_order))

        value = np.full((grid_metadata['N_of_Points'], 4), '0', dtype='a30')
        try:
            orbitals[symmetry_charakter][number_of_order] = value
        except KeyError:
            orbitals[symmetry_charakter] = {number_of_order : value}

        try:
            orbitals_metadata[symmetry_charakter][number_of_order] = {}
        except KeyError:
            orbitals_metadata[symmetry_charakter] = {number_of_order : {}}

        finally:
            current = orbitals_metadata[symmetry_charakter][number_of_order]
            if line[1] == 'Density':
                current['energy'] = np.nan
                current['occupation'] = np.nan
                current['status'] = 'not defined'

            else:
                current['energy'] = float(line[3])
                # # The following REGEX removes leading and trailing
                # # parentheses
                current['occupation'] = float(re.sub('[\(\)]', '', line[4]))
                current['status'] = line[5]


    for current_block in range(grid_metadata['N_Blocks'] - 1):
        offset = current_block * grid_metadata['Block_Size']
        for orbital in range(grid_metadata['N_of_Grids']):
            symmetry_charakter, number_of_order = order_of_orbitals[orbital]
            f.readline() # omit Title = ...
            for i in range(grid_metadata['Block_Size']):
                current_array = orbitals[symmetry_charakter][number_of_order]
                line = f.readline()
                current_array[offset + i, -1] = line

    last_block_size = (
        grid_metadata['N_P']
        - grid_metadata['Block_Size'] * (grid_metadata['N_Blocks'] - 1))
    offset = (grid_metadata['N_Blocks'] - 1) * grid_metadata['Block_Size']
    for orbital in range(grid_metadata['N_of_Grids']):
        symmetry_charakter, number_of_order = order_of_orbitals[orbital]
        f.readline() # omit Title = ...
        for i in range(last_block_size):
            current_array = orbitals[symmetry_charakter][number_of_order]
            line = f.readline()
            current_array[offset + i, -1] = line

    for symmetry_charakter in orbitals.keys():
        for number_of_order in orbitals[symmetry_charakter].keys():
            orbitals[symmetry_charakter][number_of_order] = \
            pd.DataFrame(
                orbitals[symmetry_charakter][number_of_order].astype('f8'),
                columns=['x', 'y', 'z', 'value'])

    value_to_return = {
        'grid_metadata' : grid_metadata,
        'orbitals_metadata' : orbitals_metadata,
        'orbitals' : orbitals,
        'molecule' : molecule}
    return value_to_return
