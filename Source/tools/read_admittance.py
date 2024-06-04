"""
Function to read the addmittance as obtained with the tool

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
    def __init__(self, variables, admittance, f):
        self.variables = variables  # Variable names including the block, side and d,q,dc
        self.vars = []  # Names of the variables without the block side for variable pairing
        self.y = admittance  # Admittance data
        self.f = f  # Frequency data
        # Determine if it is an AC, DC or ACDC matrix and extract block names
        self.blocks_info = {}  # Keys are the block names, contents: side (1 or 2) and type (AC, or DC)
        self.blocks = []
        y_type = []
        for name_loop in variables:
            if "dc" in name_loop.split("_")[-1]:
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
        if self.y_type == "ACDC" or (admittance.shape[1] == 2 and self.y_type == "AC") or (admittance.shape[1] == 1 and self.y_type == "DC"):
            self.node = True
        else:
            self.node = False

def read_admittance(path=None, involved_blocks=None, file_name=None, file_root=None):
    # involved_blocks is a list of strings containing the names of the blocks with sides included
    # Either admittance_type, involved_blocks and path (and name root) can be specified or full file_name
    if file_name is None:
        if (involved_blocks or path) is None:
            print('\nError: One or more required arguments are missing. \n')
            return
        else:
            # Look for the text file involving the indicated blocks
            file_name = [file for file in listdir(path) if (file.endswith("#.txt") and all(x in file for x in involved_blocks))]
            if file_root is not None: file_name = [file for file in file_name if file.startswith(file_root+"#")]
            # and (file.count("#Y_"+admittance_type+"#") > 0) indicated admittance type, admittance_type=None
            # and (file.count("#") == len(involved_blocks)+2) Just as many blocks as involved
    else:
        if path is None:
            print('\nError: File path is missing. \n')
            return

    # Read the variable names
    with open(path + '\\' + file_name[0], 'r') as f:
        variables = f.readline().strip('\n').split()

    # Load the data
    # print("Loading",path + '\\' + file_name[0])
    data = np.loadtxt(path + '\\' + file_name[0], dtype='cdouble', skiprows=1)
    freq = np.real(data[:, 0])  # Extract frequency column
    data = data[:, 1:]  # Remove frequency column

    return Admittance(variables[1:], data.reshape(data.shape[0], int(np.sqrt(data.shape[1])), int(np.sqrt(data.shape[1]))), freq)


# # Testing code
# pth = r"C:\Users\fcifuent\Desktop\KU Leuven\Projects\Z based analysis\Freq measurement\P2P test system\Tool_test\Results"
# # name = "Ytest#Y_AC#3-AC-2#4-AC-1#.txt"
# # y = read_admittance(path=pth, file_name=name)
# block = ["4-AC-2"]
# y = read_admittance(path=pth, involved_blocks=block)
# print(y.y.shape[1])
# print(y.variables)
# print(y.blocks)
# print(y.blocks_info)
# print(y.node)