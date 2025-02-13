"""
Function to read the addmittance files as obtained with the tool
"""
"""
Copyright (C) 2024  Francisco Javier Cifuentes Garcia

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

__all__ = ['read_admittance']

from os import listdir
import numpy as np  # Numerical python functions

class Admittance:
    def __init__(self, variables, admittance, f, node=False):
        self.variables = variables  # Variable names including the block, side and d,q,dc
        self.vars = []  # Names of the variables without the block side for variable pairing
        self.y = admittance  # Admittance data
        self.f = f  # Frequency data
        # Determine if it is an AC, DC or ACDC matrix and extract block names
        self.blocks_info = {}  # Keys are the block names, contents: side (1 or 2) and type (AC, or DC)
        self.blocks = []
        y_type = []
        for name_loop in variables:
            if "dc" == name_loop.split("_")[-1]:
                y_type.append("DC")
            else:
                y_type.append("AC")
            # Remove the ending: _d, _dc, _q and add it as a dict key
            b = "_".join(name_loop.split("_")[:-1])
            self.vars.append("_".join([b[:-2],name_loop.split("_")[-1]]))
            if b not in self.blocks: self.blocks.append(b)  # "-".join(b.split("-")[:-1])
            if self.blocks[-1] not in list(self.blocks_info.keys()):
                self.blocks_info[self.blocks[-1]] = {"type": y_type[-1], "side": b[-1]}  # b..split("-")[-1]
        y_type = list(set(y_type))
        if len(y_type) != 1: y_type = ["ACDC"]
        # print(self.blocks[-1],self.blocks_info[self.blocks[-1]])  # Only the name of the blocks as they appear in the matrix (.txt file)
        # print(self.vars)
        self.y_type = y_type[0]  # AC, DC or ACDC
        if node or self.y_type == "ACDC":
            #  or (admittance.shape[1] == 2 and self.y_type == "AC") or (admittance.shape[1] == 1 and self.y_type == "DC")
            self.node = True
        else:
            self.node = False

def read_admittance(path=None, involved_blocks=None, file_name=None, file_root=None):
    if file_name is None:
        if (involved_blocks or path) is None:
            print('\nError: One or more required arguments are missing. \n')
            return
        else:
            # Look for the text file involving the indicated blocks
            file_name = [file for file in listdir(path) if (file.endswith("#.txt") and all(x in file for x in involved_blocks))]
            if file_root is not None: file_name = [file for file in file_name if file.startswith(file_root+"#")]
            file_name = file_name[0]
            # and (file.count("#") == len(involved_blocks)+2) Just as many blocks as involved
    else:
        if path is None:
            print('\nError: File path is missing. \n')
            return

    # Read the variable names
    with open(path + '\\' + file_name, 'r') as f:
        variables = f.readline().strip('\n').split()

    # Load the data
    data = np.loadtxt(path + '\\' + file_name, dtype='cdouble', skiprows=1)
    freq = np.real(data[:, 0])  # Extract frequency column
    data = data[:, 1:]  # Remove frequency column

    return Admittance(variables[1:], data.reshape(data.shape[0], int(np.sqrt(data.shape[1])), int(np.sqrt(data.shape[1]))), freq)

read_admittance.__doc__ = """
Function to read the addmittance text files as obtained with the tool.

The result is an Admittance object which contains the admittance matrix, frequencies, variable names, type of matrix, etc. See the Admittance class for more information.
Either file_name or involved_blocks need to be specified.

Required arguments
        path                (string) Directory where the file is located.
        file_name           (string) Full name of the file, including extension. Either this or the involved_blocks argument is required.

Optional arguments
        involved_blocks     (list of strings) It contains the names of the PSCAD blocks with the side at the end separated with a hyphen, e.g. ["BlockA-2","BlockB-1"] denotes the admittance between side 2 of BlockA and side 1 of BlockB.
        file_root           (string) Starting part of the file name. It can be used when several study cases are performed with different names as the involved_blocks variable is the same.

Returns
        Admittance object.

Examples
# Using path and file_name
pth = r"C:/Users/fcifuent/Desktop/Z-tool/Examples/Energy_hub/Results stable"
file = "ISGT_stable#Y_AC#g2-2#g3-2#MMC2-1#MMC3-1#.txt"
admittance = read_admittance(path=pth, file_name=file)

# Using path and involved_blocks
blocks = ["g2-2","MMC2-1","MMC3-1,"g3-2"]
admittance = read_admittance(path=pth, involved_blocks=blocks)

# Access the data
print(admittance.y.shape)
print(admittance.variables)
print(admittance.blocks)
print(admittance.blocks_info)
print(admittance.node)

"""