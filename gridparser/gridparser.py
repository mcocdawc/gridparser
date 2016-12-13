import pandas as pd
import chemcoord as cc
import numpy as np
import re
import io
import line_profiler
from line_profiler import LineProfiler
# from . import _pandas_wrapper
from . import export


def split(string, seperate_newline = True):
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

            **metadata**: Metadata of the grid.

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
    metadata = {}
    orbitals_metadata = {}
    orbitals = {}

    f = open(file, 'tr')
    for _ in range(2):
        f.readline()
    line = split(f.readline())
    metadata['Natom'] = int(line[1])

    molecule_in = []
    for _ in range(metadata['Natom']):
        line = split(f.readline())
        # The following removes numbers after element symbol
        element_symbol = re.search("[a-zA-Z]", line[0]).group()
        line[0] = element_symbol
        molecule_in.append(' '.join(line))

    molecule_in = ' '.join(molecule_in)
    molecule_in = str(metadata['Natom']) + (2 * '\n') + molecule_in
    molecule_in = io.StringIO(molecule_in)
    molecule = cc.read(molecule_in, filetype='xyz')

    def get_value_and_correct_type(line):
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
        return actions[line[0]](line)

    end_of_metadata_reached = False
    while not end_of_metadata_reached:
        line = split(f.readline())
        if line[0] == 'GridName':
            end_of_metadata_reached = True
        else:
            current_line = f.tell()
            key = line[0]
            metadata[key] = get_value_and_correct_type(line)
    f.seek(current_line, 0)


    order_of_orbitals = []
    for _ in range(metadata['N_of_Grids']):
        line = split(f.readline())
        try:
            symmetry_charakter = int(line[1])
            number_of_order = int(line[2])
        except ValueError:
            # because int('Density') just does not work
            # Density is totally symmetric => symmetry_charakter=1
            # grid files are 1 indexed => 0 is magic number
            symmetry_charakter = 1
            number_of_order = 0

        order_of_orbitals.append((symmetry_charakter, number_of_order))

        value = np.full((metadata['N_of_Points'], 4), '0', dtype='a30')
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



    dimension Value(0:Net1, 0:Net2, 0:Net3,1:N_of_Grids)

i1=i2=i3=0;

for (ib=0; ib<N_Blocks; ib++)
{
 for (ig=0; ig<N_of_Grids; ig++)
 {
   read_string_title;
   ix=Block_Size; if (ib==N_Blocks) { ix=N_of_Points-Block_Size*(N_Blocks-1);}
   for (ip=0; ip<ix; ip++)
   {
     Q=read_number;
     Value[i1,i2,i3,ig]=Q;
     /* increment indexes in 3d array */
     i3++;
     if(i3 > Net1) {
         i3=0; i2++
     }
     if(i2 > Net2) {
         i2=0; i1++
    }

   }
 }
}

 (i1,i2,i3) has coords in a.u.:
 X=Orig[1]+i1*Axis(1,1)/Net1
 Y=Orig[2]+i2*Axis(2,2)/Net2
 Z=Orig[3]+i3*Axis(3,3)/Net3



    # TODO from here on performance
    last_block_size = (metadata['N_P']
                       - metadata['Block_Size'] * (metadata['N_Blocks'] - 1))

    for current_block in range(metadata['N_Blocks'] - 1):
        offset = current_block * metadata['Block_Size']
        for orbital in range(metadata['N_of_Grids']):
            symmetry_charakter, number_of_order = order_of_orbitals[orbital]
            current_array = orbitals[symmetry_charakter][number_of_order]
            f.readline() # omit Title = ...
            for i in range(metadata['Block_Size']):
                line = f.readline()
                current_array[offset + i, -1] = line

    offset = (metadata['N_Blocks'] - 1) * metadata['Block_Size']
    for orbital in range(metadata['N_of_Grids']):
        symmetry_charakter, number_of_order = order_of_orbitals[orbital]
        current_array = orbitals[symmetry_charakter][number_of_order]
        f.readline() # omit Title = ...
        for i in range(last_block_size):
            line = f.readline()
            current_array[offset + i, -1] = line

    for symmetry_charakter in orbitals.keys():
        for number_of_order in orbitals[symmetry_charakter].keys():
            orbitals[symmetry_charakter][number_of_order] = \
            pd.DataFrame(
                orbitals[symmetry_charakter][number_of_order].astype('f8'),
                columns=['x', 'y', 'z', 'value'])

    value_to_return = {
        'metadata' : metadata,
        'orbitals_metadata' : orbitals_metadata,
        'orbitals' : orbitals,
        'molecule' : molecule}
    return value_to_return
