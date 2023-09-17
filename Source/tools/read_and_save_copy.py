__all__ = ['single_s', 'multiple_s']

import numpy as np  # Numerical python functions
from os import listdir

AC_scan_variables = ['V_DUTac:1', 'V_DUTac:2', 'I_DUTacA1:1', 'I_DUTacA1:2', 'I_DUTacA2:1', 'I_DUTacA2:2']
DC_scan_variables = ['V_DUTdc', 'I_DUTdcA1', 'I_DUTdcA2']

def single_s(original_folder=None, target_filename=None, new_folder=None, output=False, save=False,
             ScanBlock=None, scantype="AC"):
    if scantype is "AC":
        variables = AC_scan_variables
    else:
        variables = DC_scan_variables
    with open(original_folder + '\\' + target_filename + ".inf", 'r') as info_file:  # Open the info file
        # Continue from here:
        # Identify the variables to be retrieved based on ScanBlock.parameters()['block_id']
        cols = [0]
        names = ["time"]
        counter = 1
        for line in info_file.readlines():
            if line.split()[2].split('"')[1] in variables:
                cols.append(counter)  # Get variable absolute position and name
                names.append(line.split()[2].split('"')[1])
            counter = counter + 1
    # In addition to the voltage/current signals we also need the time at the column 0
    values = np.loadtxt(original_folder + '\\' + target_filename + "_01.out", skiprows=1, usecols=cols)
    permutation = [0]  # The first column is always time
    for name in variables: permutation.append(names.index(name))  # Get permutation vector
    idx = np.empty_like(permutation)
    idx[permutation] = np.arange(len(permutation))
    values[:] = values[:, idx]  # Permute the columns: sorted signals
    if save: np.savetxt(new_folder + '\\' + target_filename + '.txt', values, delimiter='\t', comments='')
    if output: return values


def multiple_s(n_sim=None, original_folder=None, target_filename=None, new_folder=None, output=False, save=False,
               scantype="AC"):
    # Filter the file names to identify the target multiple output files
    files = [file for file in listdir(original_folder) if
             (file.endswith(".out") and (file.count(target_filename) > 0))]  # end in .out and contain name
    # More file filtering: file name followed by an _ and another two _ near the end
    files_filtered = [file for file in files if (file[len(target_filename)] == '_' and file.count("_", len(target_filename)) == 2 and len(file) > len(target_filename))]
    # If multiple output vars, it only reads the 1st file
    if len(files_filtered) != n_sim: files_filtered = [file for file in files_filtered if file.endswith("_01.out")]
    # Sort the files from low to high simulation: split by '_' and take the rank number (position 2 from the end)
    files_filtered.sort(key=lambda file_name: int(file_name.split('_')[-2]))
    if scantype is "AC":
        n = 6  # Number of simulation variables to be retreived (excluding time)
        variables = AC_scan_variables
    else:
        n = 3  # Number of simulation variables to be retreived (excluding time)
        variables = DC_scan_variables
    with open(original_folder + '\\' + target_filename + "_01.inf", 'r') as info_file:  # Open the first info file
        cols = [0]
        names = ["time"]
        counter = 1
        for line in info_file.readlines():
            if line.split()[2].split('"')[1] in variables:
                cols.append(counter)  # Get variable position (column) and names
                names.append(line.split()[2].split('"')[1])
            counter = counter + 1
    # Only for the first file, it loads also the first column (time)
    first_values = np.loadtxt(original_folder + '\\' + files_filtered[0], skiprows=1, usecols=cols)
    permutation = [0]
    for name in variables: permutation.append(names.index(name))  # Get permutation vector
    idx = np.empty_like(permutation)
    idx[permutation] = np.arange(len(permutation))
    first_values[:] = first_values[:, idx]  # Permute the columns: sorted signals acordingly to the names in "variables"
    values = np.empty((first_values.shape[0], n * len(files_filtered) + 1), dtype='d')  # Preallocation of memory
    values[:, 0:n + 1] = first_values  # Save the first file
    del first_values, files  # Clear some memory

    # Compute the permutation indexes for the rest of the files
    permutation = [k-1 for k in permutation[1:]]  # Get permutation vector without the time column (first)
    idx = np.empty_like(permutation)
    idx[permutation] = np.arange(len(permutation))

    j = n + 1  # File counter
    for file in files_filtered[1:n_sim]:
        # Save the rest without the time / first column
        data = np.loadtxt(original_folder + '\\' + file, skiprows=1, usecols=cols[1:])  # Not time column
        values[:, j:j + n] = data[:, idx]  # Permute the columns: sorted signals acordingly to the names in "variables"
        j += n
    if save: np.savetxt(new_folder + '\\' + target_filename + '.txt', values, delimiter='\t', comments='')
    if output: return values


multiple_s.__doc__ = """
Function that reads Multiple Simulations results and saves them into a dedicated file and/or loads them into memory for further processing.
The function accepts several input arguments to customize the reading:
Required
        n_sim		 Number of simulation results to be read: total number of simulations.
        original_folder	 Absolute path of the folder where the results can be found.
        target_filename	 Common name of the files to be read.       

Optional
        output		 Boolean variable to control if the read results should be returned as an array. output = True returns an array with the results.
        save		 Boolean variable to control the saving of the read data. False means the data is not saved.
        new_folder	 Absolute path of the destination folder in case save = True.
        output_filename	 Output file name: no extension needed. By default it saves the data into .txt files where the first column is the time and the rest
                         are the concatenated extracted variables from each simulation in the order that they where performed.

Furthermore, there is an additional internal parameter that indicates how many variables are read from each file. For impedance computation
only two currents and voltages are needed, thus n = 4 variables.

"""
single_s.__doc__ = """
Function that reads Single Simulation results and saves them into a dedicated file and/or loads them into memory for further processing.
The function accepts several input arguments to customize the reading:
Required
        original_folder	 Absolute path of the folder where the results can be found.
        target_filename	 Common name of the file to be read.       

Optional
        output		 Boolean variable to control if the read results should be returned as an array. output = True returns an array with the results.
        save		 Boolean variable to control the saving of the read data. False means the data is not saved.
        new_folder	 Absolute path of the destination folder in case save = True.
        output_filename	 Output file name: no extension needed. By default it saves the data into .txt files where the first column is the time and the rest
                         are the extracted variables.

Furthermore, there is an additional internal parameter that indicates how many variables are read from each file. For impedance computation
only two currents and voltages are needed, thus n = 4 variables.

"""
