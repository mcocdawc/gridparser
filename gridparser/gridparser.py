import pandas as pd
import chemcoord as cc
import numpy as np
import re
import io
# import line_profiler
# from line_profiler import LineProfiler
# from . import _pandas_wrapper
from . import export

from scipy.constants import physical_constants

BOHR_TO_A = physical_constants['Bohr radius'][0] * 1e10

def split(string, seperate_newline = True):
    if seperate_newline:
        return list(filter(None, re.split("[, =]+|(\\n)+", string)))
    else:
        return list(filter(None, re.split("[, =]+", string)))

@export
class Grid():
    def __init__(self, structure, metadata, orbitals):
        self.structure  = structure
        self.metadata = metadata
        self._orbital_template = self._give_orbital_template()
        self._orbitals = orbitals

    def _give_orbital_template(self):
        fractional_coord = np.empty([self.metadata['N_of_Points'], 3],
                                    dtype='f8')
        i1 = i2 = i3 = 0
        for row in range(fractional_coord.shape[0]):
            fractional_coord[row, :3] = [i3, i2, i1]
            i1 = i1 + 1
            if i1 > self.metadata['Net'][0]:
                i1 = 0
                i2 = i2 + 1
            if i2 > self.metadata['Net'][1]:
                i2 = 0
                i3 = i3 + 1

        basis = self.metadata['Axis'] / (self.metadata['Net'] + 1)
        increment = (basis[None, :] * fractional_coord[:, None]).sum(axis=1)
        location = (increment + self.metadata['Origin']) * BOHR_TO_A

        orbital = pd.DataFrame(np.empty([self.metadata['N_of_Points'], 4]),
                               columns=['x', 'y', 'z', 'value'])
        orbital.loc[:, 'x' : 'z'] = location
        return orbital

    def __repr__(self):
        treat_density = 0 in self._orbitals[1]
        string_list = ['1 Electronic density\n'] if treat_density else []
        for symm_char in self._orbitals.keys():
            N_orb = len(self._orbitals[symm_char])
            N_orb = N_orb - 1 if symm_char == 1 and treat_density else N_orb
            text = '{0} orbitals in symmetry character {1}\n'.format(
                N_orb, symm_char)
            string_list.append(text)
        return ''.join(string_list)

    def __add__(self, other):
        pass

    def __radd__(self, other):
        pass

    def __rmatmul__(self, other):
        pass

    def give_orbital(self, symmetry_char, iorb):
        orbital = self._orbital_template.copy()
        orbital['value'] = self._orbitals[symmetry_char][iorb]
        return orbital
        # return orbital, energy

    @classmethod
    def parse_grid(cls, file):
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
#          molecule = cc.read(molecule_in, filetype='xyz')
        molecule = cc.Cartesian.read_xyz(molecule_in)

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
        metadata['Axis'] = np.array([metadata[axis] for axis in ('Axis_1', 'Axis_2', 'Axis_3')]).T
        for axis in ('Axis_1', 'Axis_2', 'Axis_3'):
            del metadata[axis]

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

            value = np.full(metadata['N_of_Points'], '0', dtype='a30')
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



        last_block_size = (metadata['N_P']
                           - metadata['Block_Size'] * (metadata['N_Blocks'] - 1))
        for ib in range(metadata['N_Blocks']):
            offset = ib * metadata['Block_Size']
            for ig in range(metadata['N_of_Grids']):
                symmetry_charakter, number_of_order = order_of_orbitals[ig]
                current_array = orbitals[symmetry_charakter][number_of_order]
                f.readline() # omit Title = ...
                if ib == (metadata['N_Blocks'] - 1):
                    ix = last_block_size
                else:
                    ix = metadata['Block_Size']
                for ip in range(ix):
                    current_array[offset + ip] = f.readline()

        for sym_char in orbitals.keys():
            for iorb in orbitals[sym_char].keys():
                orbitals[sym_char][iorb] = orbitals[sym_char][iorb].astype('f8')
        # return orbitals, metadata
        # return metadata
        return cls(molecule, metadata, orbitals)
        # return orbitals

        # for symmetry_charakter in orbitals.keys():
        #     for number_of_order in orbitals[symmetry_charakter].keys():
        #         orbitals[symmetry_charakter][number_of_order] = \
        #         pd.DataFrame(
        #             orbitals[symmetry_charakter][number_of_order].astype('f8'),
        #             columns=['x', 'y', 'z', 'value'])
        #
        # value_to_return = {
        #     'metadata' : metadata,
        #     'orbitals_metadata' : orbitals_metadata,
        #     'orbitals' : orbitals,
        #     'molecule' : molecule}
        # return value_to_return




@export
class Plane:
    __array_priority__ = 10.
    def __init__(self, r, normal_v):
        self.r = r
        self.normal_v = normal_v / np.linalg.norm(normal_v)

    def __repr__(self):
        return "Hesse normal form\nr:\n{0}\nnormal_v:\n{1}".format(
            self.r.__repr__(), self.normal_v.__repr__())

    def __add__(self, other):
        return self.__class__(self.r + other, self.normal_v)

    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        return self.__class__(self.r * other, self.normal_v)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __rmatmul__(self, other):
        det = np.linalg.det(other)
        if np.isclose(det, [-1, 1]).any():
            return self.__class__(other @ self.r, other @ self.normal_v)
        else:
            raise TypeError('Only orthogonal transformations are defined.')

    @classmethod
    def span(cls, r, v1, v2):
        return cls(r, np.cross(v1, v2))
