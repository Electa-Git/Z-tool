__all__ = ['single_s', 'multiple_s']

"""Documentation for this function is not updated"""

import numpy as np  # Numerical python functions
from os import listdir

AC_scan_variables = ['V_DUTac:1', 'V_DUTac:2', 'I_DUTacA1:1', 'I_DUTacA1:2', 'I_DUTacA2:1', 'I_DUTacA2:2']
DC_scan_variables = ['V_DUTdc', 'I_DUTdcA1', 'I_DUTdcA2']

def single_s(out_files=None, save_folder=None, save=False, files=None, zblocks=None, new_file_name=None):
    save_time = True  # Save the time vector only once for space-saving
    for file_num in files:
        # Load each target file and for each block related to the file asign the corresponding data to the block
        if file_num < 10:  # If the file number is less than 10, then it adds 0 before the file number
            values = np.loadtxt(out_files + "_0" + str(file_num) + ".out", skiprows=1)
        else:
            values = np.loadtxt(out_files + "_" + str(file_num) + ".out", skiprows=1)
        for block in zblocks:
            if save_time:
                block.snapshot_data["time"] = values[:, 0]  # Retreive the time vector
                save_time = False  # Only for the first scanning block (to save mem space)
            for block_file in block.files_to_open:
                if file_num == block_file:
                    # If the block needs data from the file, then use the block's target columns for this file
                    for col in block.relative_cols[file_num]:
                        ch = col - 1 + 10 * (file_num - 1)  # Absolute output channel number
                        # print(block.out_vars_names[ch], ch_var_names[ch])
                        block.snapshot_data[block.out_vars_names[ch]] = values[:, col - 1]  # Retreived data

    if save:
        new_file_name = save_folder + '\\' + new_file_name + '.txt'
        var_names = ["time"]
        data = zblocks[0].snapshot_data["time"].reshape(-1, 1)  # Retreive the time vector data
        for block in zblocks:
            for name in list(block.out_vars_names.values()):
                if name != "time":  # Do not save the time vector
                    data = np.append(data, block.snapshot_data[name].reshape(-1, 1), axis=1)
                    var_names.append(name)
        np.savetxt(new_file_name, data, delimiter='\t', header="\t".join(var_names))

def multiple_s(n_sim=None, out_folder=None, file_name=None, save_folder=None, save=False, tar_files=None, zblocks=None):
    # # Filter the file names to identify the target multiple output files
    # files = [file for file in listdir(out_folder) if
    #          (file.endswith(".out") and (file.count(file_name) > 0))]  # end in .out and contain name
    # # More file filtering: file name followed by an _ and another two _ near the end
    # files_filtered = [file for file in files if (file[len(file_name)] == '_' and
    #                                              file.count("_", len(file_name)) == 2 and
    #                                              len(file) > len(file_name))]
    # # Sort the files from low to high simulation: split by '_' and take the rank number (position 2 from the end)
    # files_filtered.sort(key=lambda file_name: int(file_name.split('_')[-2]))
    # for i in files_filtered: print(i)  # Just for debugging

    # New function
    save_time = True  # Save the time vector only once improving the memory usage
    root_name = out_folder + '\\' + file_name  # Root of the filename
    sim_type = "_"+file_name.split("_")[-1]  # Use the end of the file name (sim type) to re-name the output variables
    # print("Perturbation type: ",sim_type)
    for sim in range(1,n_sim+1):
        # For each simulation
        for file_num in tar_files:
            # Load each target file and for each block related to the file and asign the corresponding data to the block
            if file_num < 10:  # If the file number is less than 10, then it adds 0 before the file number
                if sim < 10:  # If the simulation number is less than 10, then it adds a 0 before the simultation number
                    values = np.loadtxt(root_name + "_0" + str(sim) + "_0" + str(file_num) + ".out", skiprows=1)
                else:
                    values = np.loadtxt(root_name + "_" + str(sim) + "_0" + str(file_num) + ".out", skiprows=1)
            else:
                if sim < 10:
                    values = np.loadtxt(root_name + "_0" + str(sim) + "_" + str(file_num) + ".out", skiprows=1)
                else:
                    values = np.loadtxt(root_name + "_" + str(sim) + "_" + str(file_num) + ".out", skiprows=1)
            if save_time:
                zblocks[0].perturbation_data["time"] = values[:, 0]  # Retreive the time vector
                print(" Time vector start",zblocks[0].perturbation_data["time"][0],"and end",zblocks[0].perturbation_data["time"][-1])
                save_time = False  # Only for the first scanning block (to save memory)
            for block in zblocks:
                # For each z-tool block related to the file, asign the corresponding data to the block
                for block_file in block.files_to_open:
                    if file_num == block_file:
                        # print("Simulation",sim, "block:", block.name, 'File:', block_file)
                        # If the block needs data from the file, then use the block's target columns for this file
                        for col in block.relative_cols[file_num]:
                            ch = col - 1 + 10 * (file_num - 1)  # Absolute output channel number for the column
                            # print(" Channel",ch,"variable",block.out_vars_names[ch],"initial value",values[0, col - 1])
                            block.perturbation_data[sim-1][block.out_vars_names[ch]+sim_type] = values[:, col - 1]  # Data

    if save:
        file_name = save_folder+'\\'+file_name
        for sim in range(n_sim):
            var_names = ["time"]
            data = zblocks[0].perturbation_data["time"].reshape(-1, 1)  # Retreive the time vector data
            for block in zblocks:
                for name in list(block.out_vars_names.values()):
                    if name != "time":  # Do not include the time vector because it is already added
                        data = np.append(data, block.perturbation_data[sim][name+sim_type].reshape(-1, 1), axis=1)
                        var_names.append(name)
            np.savetxt(file_name + "_" + str(sim) + '.txt', data, delimiter='\t', header="\t".join(var_names))

multiple_s.__doc__ = """
Function that reads Multiple Simulations results and saves them into a dedicated file and/or loads them into memory for further processing.
The function accepts several input arguments to customize the reading:
Required
        n_sim		 Number of simulation results to be read: total number of simulations.
        original_folder	 Absolute path of the folder where the results can be found.
        file_name	 Common name of the files to be read.       

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
        file_name	 Common name of the file to be read.       

Optional
        output		 Boolean variable to control if the read results should be returned as an array. output = True returns an array with the results.
        save		 Boolean variable to control the saving of the read data. False means the data is not saved.
        new_folder	 Absolute path of the destination folder in case save = True.
        output_filename	 Output file name: no extension needed. By default it saves the data into .txt files where the first column is the time and the rest
                         are the extracted variables.

Furthermore, there is an additional internal parameter that indicates how many variables are read from each file. For impedance computation
only two currents and voltages are needed, thus n = 4 variables.

"""
