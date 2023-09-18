""" Documentation is not updated """

import time as t  # Relative time
from os import listdir, chdir, getcwd, path, makedirs
import numpy as np  # Numerical python functions
from mhi.pscad import launch  # PSCAD automation functions
from Source.tools import *

# Output channels for each scanning block: blockid is assummed to be the first output when retreiving the data /!\
AC_scan_variables = ['blockid', 'IDUTacA1:1', 'IDUTacA1:2','VDUTac:1', 'VDUTac:2', 'IDUTacA2:1', 'IDUTacA2:2','theta']
DC_scan_variables = ['blockid', 'IDUTdcA1', 'IDUTdcA2','VDUTdc']

class Network:
    def __init__(self, name_blocks_involved, scan_type, adj_matrix):
        self.names = name_blocks_involved  # List with full names including the sides
        self.scan_type = scan_type  # AC or DC scan
        self.adj_matrix = adj_matrix  # The adjacent matrix = zeros correspond to y = 0 (disconnection)
        if scan_type == "AC":
            self.runs = 2*len(name_blocks_involved)  # Number of needed runs for the network scan
            perturbations = ["_d","_q"]
        else:
            self.runs = len(name_blocks_involved)  # Number of needed runs for the network scan
            perturbations = ["_dc"]
        self.blocks_idx = None  # Dict key: self.names, pointing at the key of the associated blocks in ScanBlocksTool
        self.remaining_scans = [sources+scan_type for sources in self.names for scan_type in perturbations]
        self.all_scans = [sources + scan_type for sources in self.names for scan_type in perturbations]
        self.enforce = True  # Enforce network connectivity when computing the admittance

    def updateScans(self,done):
        self.remaining_scans.remove(done)

class Graph:
    def __init__(self, V):
        self.V = V
        self.adj = [[] for i in range(V)]

    # Depth-first search algorithm method
    def DFSUtil(self, temp, v, visited):
        # Mark the current vertex as visited
        visited[v] = True
        # Add the vertex to list
        temp.append(v)
        # Repeat for all vertices adjacent to this vertex v
        for i in self.adj[v]:
            if not visited[i]:
                # Update the list
                temp = self.DFSUtil(temp, i, visited)
        return temp

    # Add an undirected edge
    def addEdge(self, v, w):
        self.adj[v].append(w)
        self.adj[w].append(v)

    # Method to retrieve connected components in an undirected graph
    def connectedComponents(self):
        visited = []
        cc = []
        for i in range(self.V):
            visited.append(False)
        for v in range(self.V):
            if not visited[v]:
                temp = []
                cc.append(self.DFSUtil(temp, v, visited))
        return cc

class Scanblock:
    def __init__(self, pscad_block, name, block_id):
        self.pscad_block = pscad_block
        self.name = name
        self.block_id = block_id     # Unique block identifier number
        self.out_vars_ch = None      # Output variables channel number
        self.out_vars_names = {}     # Output variables channel name keyed by absolute channel number
        self.files_to_open = None    # Files' number containing the data of each scan block
        self.relative_cols = {}      # Dictionary of relative columns lists for each file; keys = files_to_open
        self.snapshot_data = dict()  # Snapshot recordings, the keys are the signal names out_vars_names[ch]
        self.perturbation_data = None  # Dict of dicts: 1-key = sim#, 2-keys = names out_vars_names[ch] +"_d","_q",...
        if "AC" in self.pscad_block.defn_name[1]:
            self.type = "AC"
            self.var_names = ['VDUTac', 'IDUTacA1', 'IDUTacA2', 'theta']  # Root of the variable names
            self.group = "ACscan"
        else:
            self.var_names = ['VDUTdc', 'IDUTdcA1', 'IDUTdcA2']  # Root part of the variable names
            self.type = "DC"
            self.group = "DCscanPM"

def frequency_sweep_ACDC_full(t_snap=None, t_sim=None, t_step=None, sample_step=None, v_perturb_mag=None,
                              freq=None, f_points=None, f_base=None, f_min=None, f_max=None, working_dir=None,
                              freq_text_file='frequencies.txt', snapshot_file='Snapshot',
                              take_snapshot=True, dt_injections=None, topology=None,
                              project_name='DUT', workspace_name='DUTscan', scanid=None,
                              fortran_ext=r'.gf46', num_parallel_sim=8, component_parameters=None,
                              results_folder=None, output_files='Perturbation', compute_yz=False, save_td=False,
                              fft_periods=1, start_fft=None):
    # Debugging control
    run_sim = False
    verbose = True
    """ --- Input data handling --- """
    # CHECK the following if: it does not work... "or" operator gives a bool as output not None
    if (t_snap or t_step or start_fft or v_perturb_mag or ((f_points or f_base or f_max or f_min) and freq)) is None:
        print(
            'One or more required arguments are missing!! \n Check the function documentation by typing help('
            'frequency_sweep) \n')
        return
    # If the sample time is not provided, it is set to half of the minimum required value (multiple of step_time)
    if sample_step is None:
        if f_max is None: f_max = max(freq)
        sample_step = round(t_step * np.floor((1e6 * 0.5 * 0.5 / f_max + t_step / 2) / t_step), 3)  # [us]
    if t_sim is None: t_sim = start_fft + fft_periods / f_base
    if dt_injections is None: dt_injections = start_fft  # Time to reach steady-state after the system decoupling
    if working_dir is None:
        working_dir = getcwd() + '\\'  # Location of the PSCAD workspace
    else:
        chdir(working_dir.encode('ascii', 'backslashreplace'))  # ('unicode_escape'))  Location of the PSCAD workspace
        working_dir = getcwd() + '\\'
    print('\nRunning from ' + working_dir + '\n')
    out_dir = working_dir + project_name + fortran_ext  # Output files directory
    # The snapshot and simulation times must be a multiple of the sampling time
    t_snap_internal = round(np.ceil((t_snap + dt_injections) / (sample_step * 1e-6)) * sample_step * 1e-6, 6)
    t_sim_internal = round(np.ceil(t_sim / (sample_step * 1e-6)) * sample_step * 1e-6, 6)

    # Create frequency list if it is not provided
    if freq is None: freq = create_freq.loglist(f_min=f_min, f_max=f_max, f_points=f_points, f_base=f_base)
    # Frequency vector to XY table
    f_points = len(freq)
    with open(freq_text_file, 'w') as f:  # Create the .txt file
        f.write('! This file stores the frequency swept values \n')  # Write header
        for j in range(f_points): f.write(str(j + 1) + '\t' + str(freq[j]) + '\n')  # Write values
        f.write('ENDFILE:')  # Write end of the file
    f.close()
    # Create the folder if it does not exist
    if (results_folder is not None) and (not path.exists(results_folder)): makedirs(results_folder)
    if (results_folder is None) and (save_td or compute_yz): results_folder = working_dir

    """ --- Main program --- """
    t0 = t.time()  # Initial time
    print('Launching PSCAD')
    pscad = launch(minimize=True)
    wait4pscad(time=1, pscad=pscad)
    t.sleep(5)  # Wait a bit more just in case PSCAD is still loading some stuff
    if not pscad.licensed():
        certificates = pscad.get_available_certificates()
        keys = list(certificates.keys())
        pscad.get_certificate(certificates[keys[0]])  # Get the first available certificate
    pscad.load(working_dir + workspace_name + ".pswx")  # Load workspace
    project = pscad.project(project_name)  # Where the model is defined
    main = project.canvas('Main')  # Main of the project
    project.focus()  # Switch PSCAD’s focus to the project

    print(' Setting parameters and configuration')
    # Get components and set global parameters for all simulations
    if scanid is None:
        # If a single block is in the canvas, only that one needs to be scanned
        blocks = main.find_first("Z_tool:ACscan")  # This assumes a single scan block in the main canvas
        if blocks is None: blocks = main.find_first("Z_tool:DCscanPM")  # If it did not find a ACscan it is DCscan
        scanid = blocks.parameters()['Name']  # Retrieve the scan block name for identification
        # The rest of the code works well based on scanid

    ScanBlocksAC = []
    ScanBlocksDC = []
    for identification in scanid:
        blocks = main.find_all(Name=identification)
        blocks_tool = [block for block in blocks if "Z_tool" in block.defn_name[0]]  # Filter only Z-tool components
        if "DC" in blocks_tool[0].defn_name[1]:
            ScanBlocksDC.append(blocks_tool[0])
        else:
            ScanBlocksAC.append(blocks_tool[0])
        # if len(blocks_tool) > 1:  # If more than one scan block shares name with
        #     blocksname = blocks_tool[0].defn_name[1] + blocks_tool[1].defn_name[1]  # Second part of the block name
        #     if "AC" in blocksname and "DC" in blocksname:
        #         ScanBlocksACDC.append(blocks_tool[0])
        #         ScanBlocksACDC.append(blocks_tool[1])
        #         blocks_tool[0].parameters(Tdecoupling=t_snap, T_inj=t_snap_internal, selector=0)
        #         blocks_tool[1].parameters(Tdecoupling=t_snap, T_inj=t_snap_internal, selector=0)
        #     else:
        #         print('Error: The list of scaning blocks is inconsistent. \n')
        #         return
        # elif "DC" in blocks_tool[0].defn_name[1]:
        #     ScanBlocksDC.append(blocks_tool[0])
        #     blocks_tool[0].parameters(Tdecoupling=t_snap, T_inj=t_snap_internal, selector=0)
        # else:
        #     ScanBlocksAC.append(blocks_tool[0])
        #     blocks_tool[0].parameters(Tdecoupling=t_snap, T_inj=t_snap_internal, selector=0)
    # ScanBlocksAC = list(set(ScanBlocksAC))
    # ScanBlocksDC = list(set(ScanBlocksDC))

    ScanBlocksAC_names = [block.parameters()['Name'] for block in ScanBlocksAC]

    if ScanBlocksAC and ScanBlocksDC:
        scantype = "ACDC"
        ScanBlocks = ScanBlocksAC + ScanBlocksDC
        group = "ACscan" + "DCscanPM"
    elif ScanBlocksAC:
        scantype = "AC"
        ScanBlocks = ScanBlocksAC
        group = "ACscan"
    else:
        scantype = "DC"
        ScanBlocks = ScanBlocksDC
        group = "DCscanPM"
    print(' Type of scan:', scantype)
    ScanBlocks.sort(key=lambda x: x.parameters()['Name'][:-3], reverse=False)  # Sort the blocks by their "bus" number
    ScanBlocks_id = [i for i in range(1, len(ScanBlocksAC) + len(ScanBlocksDC) + 1)]  # Unique scan block_id signals
    # Create a list with the active scan block objects containing rich information about each block
    ScanBlocksTool = []
    ScanBlocksTool_names = []  # Name of the blocks as they appear in ScanBlocksTool to avoid excesive iterations
    for idx, block in enumerate(ScanBlocks):
        # Set snapshot parameters and block ID in the scan blocks
        block.parameters(Tdecoupling=t_snap, T_inj=t_snap_internal, selector=0, block_id=ScanBlocks_id[idx])
        ScanBlocksTool.append(Scanblock(block, block.parameters()['Name'], int(block.parameters()['block_id'])))
        if verbose: print(" Scan block type ",block.defn_name[1])
        ScanBlocksTool[idx].perturbation_data = {i: {} for i in range(f_points)}  # Dict of dicts
        ScanBlocksTool_names.append(ScanBlocksTool[idx].name)
        if verbose: print(" Scan block names",ScanBlocksTool[idx].name)
        # The following lines identify the ACDC scan points
        # ScanBlocks_type.append(ScanBlocksTool[idx].type)  # List with block's scan type
        # Alternative: if "DC" in block.defn_name[1]: ScanBlocks_type.append("DC") # No need for AC or DC in block name
        # if "AC" in block.defn_name[1]: ScanBlocks_type.append("AC")
        # if idx > 0:
        #     if block.parameters()['Name'][:-3] == ScanBlocks[idx - 1].parameters()['Name'][:-3]:
        #         # If two buses have the same number, then it is an ACDC bus
        #         ScanBlocks_type[idx] = "ACDC"
        #         ScanBlocks_type[idx - 1] = "ACDC"

    # Read the topology matrix
    Ytopology = np.loadtxt(topology, skiprows=1, comments=["#","%","!"])
    # len(ScanBlocksTool)*2 by len(ScanBlocksTool)*2 # nameA-1 nameA-2 nameB-1 nameB-2 ... x nameA-1 nameA-2 nameB-1 ...
    # 0 means no interconnection, 1 means connection between the edges: diagonals are single-sided / shunt
    with open(topology, 'r') as f:
        block_names_Y = f.readline().strip('\n').split()
    if verbose: print(block_names_Y)

    # Create the undirected graph - adjacent matrix but diagonals can be 1
    g = Graph(len(Ytopology))
    for row, name in enumerate(block_names_Y):
        if verbose: print("Block names based on topology file", name)
        for col, edge in enumerate(Ytopology[row]):
            if int(edge) == 1: g.addEdge(row, col)

    # Obtain the connected components of the graph
    cc = g.connectedComponents()  # List of lists with blocks # involved in each scan

    # Extract the interconnected blocks for network scan (remove shunts a.k.a. unconnected vertices)
    scans_network = [c for c in cc if len(c) != 1]  # List of lists with blocks # in each scan with more than 1 block
    # After this for loop, the list is filtered & enhanced by a network scans object: num of runs, names, adj matrix...
    passive_networks_scans = []  # This variable stores said lists
    for idx, net in enumerate(scans_network):
        network_names = [block_names_Y[element] for element in net]
        if len(net) > 2:
            # Multiterminal network
            if network_names[0][:-2] in ScanBlocksAC_names:
                passive_networks_scans.append(Network(network_names, "AC", Ytopology[net, :][:, net]))
                if verbose: print("   AC network scan involving",network_names)
            else:
                passive_networks_scans.append(Network(network_names, "DC", Ytopology[net, :][:, net]))
                if verbose: print("   DC network scan involving", network_names)
        else:
            # Point to point: it needs to check that it is not a AC/DC converter
            types = []
            for name in network_names:
                if name[:-2] in ScanBlocksAC_names:
                    types.append("AC")
                else:
                    types.append("DC")
            if types[0] == types[1]:
                # If the block type are the same, then they are interconnecting a AC or DC network
                passive_networks_scans.append(Network(network_names, types[0], Ytopology[net, :][:, net]))
                if verbose: print("   "+types[0]+" network scan involving", network_names)

    # Retrieve the indexes of ScanBlocksTool for each network based on ScanBlocksTool_names
    for net in passive_networks_scans:
        net.blocks_idx = {name: ScanBlocksTool_names.index(name[:-2]) for name in net.names}  # Assumes no name duplicate

    runs = [net.runs for net in passive_networks_scans]
    max_runs = max(runs)
    bottleneck_scan = passive_networks_scans[runs.index(max_runs)]

    del scans_network  # Clean some useless variables here!

    if component_parameters is not None:
        Parameters = [main.find('master:const', 'Param1'), main.find('master:const', 'Param2'),
                      main.find('master:const', 'Param3'), main.find('master:const', 'Param4'),
                      main.find('master:const', 'Param5'), main.find('master:const', 'Param6')]
        for i in range(min(len(component_parameters), len(Parameters))): Parameters[i].parameters(
            Value=component_parameters[i])

    # Iterate over output channel components to disable all but Z-tool's scaning blocks
    all_pgb = project.find_all("master:pgb")  # Find all output channels in the project
    scan_vars = ['blockid', 'VDUTac', 'IDUTacA1', 'IDUTacA2', 'VDUTdc', 'IDUTdcA1', 'IDUTdcA2', 'theta']  # Target outputs
    for pgb in all_pgb:
        if not (pgb.parameters()['Name'] in scan_vars):  pgb.disable()  # Disable the non-selected outputs

    # Set simulation-specific parameters
    if 'Perturbation' not in pscad.simulation_sets():
        simset = pscad.create_simulation_set('Perturbation')
        simset.add_tasks(project_name)
        print(' A simulation set has been created')
    else:
        simset = pscad.simulation_set('Perturbation')
    simset_task = simset.tasks()[0]  # The task is extracted

    if take_snapshot:  # It runs the snapshots
        print(' Running snapshot')
        t1 = t.time()
        simset_task.parameters(volley=1, affinity_type='DISABLE_TRACING',
                               ammunition=1)  # affinity_type = 'DISABLE_TRACING' disables the plotting
        simset_task.overrides(duration=t_snap_internal + t_sim_internal, time_step=t_step, plot_step=sample_step,
                              start_method=0, timed_snapshots=1, snapshot_file=snapshot_file + '.snp',
                              snap_time=t_snap_internal, save_channels_file=snapshot_file + '.out', save_channels=1)
        if run_sim: simset.run()
        print(' Snapshot completed in', round((t.time() - t1), 2), 'seconds')
    else:  # It performs the unperturbed simulation starting from the given snapshot (not fully tested yet)
        print(' Running steady-state simulation')
        t1 = t.time()
        ScanBlock.parameters(selector=0)  # No injection for the steady state
        simset_task.parameters(volley=1, affinity_type='DISABLE_TRACING',
                               ammunition=1)  # affinity_type = 'DISABLE_TRACING' disables the plotting
        simset_task.overrides(duration=t_sim_internal, time_step=t_step, plot_step=sample_step, start_method=1,
                              timed_snapshots=0, startup_inputfile=snapshot_file + '.snp',
                              save_channels_file=snapshot_file + '.out', save_channels=1)
        if run_sim: simset.run()
        print(' Steady-state simulation completed in', round((t.time() - t1), 2), 'seconds')

    if take_snapshot and (save_td or compute_yz):
        wait4pscad(time=1, pscad=pscad)
        t1 = t.time()
        ch_var_names = dict()
        # Identify the variables to be retrieved, their output channels and associated Z-scan blocks
        out_filename = out_dir + "\\" + simset_task.overrides()['save_channels_file'][:-4]  # Snapshot output filename
        with open(out_filename + ".inf", 'r') as info_file:
            out_num = []  # Ztool's variables output channel number
            names = []  # Ztool's variables output names
            counter = 1  # Total number of PSCAD output signals
            for line in info_file.readlines():
                if line.split()[3].split('"')[1] in group:  # If the output channel corresponds to a Ztool variable
                    out_num.append(counter)  # Get variable's output channel number
                    # Same as out_num.append(int(line.split()[0].split('(')[1].split(')')[0]))
                    names.append(line.split()[2].split('"')[1])  # Get output variable name
                    # ch_var_names.__setitem__(key=counter, value=names[-1])  # Var name entry with the channel num as key
                    ch_var_names[counter] = names[-1]  # Var name entry with the channel num as key
                counter = counter + 1
        block_id_out_num = [out_num[i] for i in range(len(names)) if "blockid" in names[i]]  # block_id outputs numbers
        # out_files = int(np.ceil(counter / 10))  # Only 10 output channels per .out file
        files_to_open = [int(np.ceil(block_id_out / 10)) for block_id_out in block_id_out_num]  # With block_id outputs
        files_to_open = list(set(files_to_open))  # File's number to be opened
        block_id_out_signal = []
        for file_num in files_to_open:
            # Select the columns to be read relative to each file and only for the block_id signal
            cols = [num + 1 - 10 * (file_num - 1) for num in block_id_out_num if
                    int(np.ceil(num / 10)) == file_num]
            if file_num < 10:  # If the file number is less than 10, then it adds 0 before the file number
                # values = np.loadtxt(out_filename + "0" + str(file_num) + ".out", skiprows=1, max_rows=2, usecols=cols)
                values = read_one_line(out_filename + "_0" + str(file_num) + ".out", nline=1)  # Read the first value
            else:
                # values = np.loadtxt(out_filename + str(file_num) + ".out", skiprows=1, max_rows=2, usecols=cols)
                values = read_one_line(out_filename + "_" + str(file_num) + ".out", nline=1)  # Read the first value
            for idx, signal in enumerate(values):
                if idx + 1 in cols:
                    block_id_out_signal.append(int(float(signal)))

        # Map the block_id_out_signal to the blocks in the list ScanBlocksTool
        all_files_to_open = []  # A list containing all the output files number that need to be read
        for idx, id_signal in enumerate(block_id_out_signal):  # Loop over the identification signals
            for block in ScanBlocksTool:  # And check for each active scaning block
                # If the ids match, then define the first and last output channel numbers for the measurement block
                if block.block_id == id_signal:
                    ch0 = block_id_out_num[idx]  # Start channel
                    if verbose: print("Block name:",block.name,"type",block.type)
                    if "AC" == block.type:
                        ch1 = ch0 + len(AC_scan_variables)  # End channel number containing the scan block signals
                    else:
                        ch1 = ch0 + len(DC_scan_variables)  # Idem but for DC scan blocks
                    block.out_vars_ch = [i for i in range(ch0,ch1)]  # Output channel numbers for this block

                    # Dict with output channel number as the key and output channel name as the content
                    for ch in block.out_vars_ch:
                        name_ch = ch_var_names[ch].split('_')
                        if len(name_ch) > 1:
                            # There is an underscore
                            if ":" in name_ch[1]:
                                # Remove the additional numbering and add back the end of the name
                                block.out_vars_names[ch] = name_ch[0] + ":" + name_ch[1].split(":")[1]
                            else:
                                block.out_vars_names[ch] = name_ch[0]
                        else:
                            block.out_vars_names[ch] = name_ch[0]
                    if verbose:
                        print(" Block ", block.name, ": target output channels ", block.out_vars_ch)
                        print(" Block ", block.name, ": output names ", [block.out_vars_names[ch] for ch in block.out_vars_ch])

                    # block.out_vars_names.__setitem__(key=ch, value=ch_var_names.get(ch))
                    files_to_open = [int(np.ceil(block_out / 10)) for block_out in block.out_vars_ch]
                    files_to_open = list(set(files_to_open))  # File's number to be opened (no repetitions)
                    block.files_to_open = files_to_open  # Number of the files containing the block's outputs
                    for f2o in files_to_open: all_files_to_open.append(f2o)
                    for file_num in files_to_open:
                        # Select the columns to be read relative to each file and only for the signals of the scan block
                        cols = [num + 1 - 10 * (file_num - 1) for num in block.out_vars_ch if
                                int(np.ceil(num / 10)) == file_num]
                        block.relative_cols[file_num] = cols
                        if verbose: print(" Block ", block.name, " output file:", file_num,", columns:",cols)

        all_files_to_open = list(set(all_files_to_open))  # Get rid of repetitions

        # Save snapshot run results
        read_and_save.single_s(out_files=out_filename, save_folder=results_folder,
                               save=save_td, files=all_files_to_open, zblocks=ScanBlocksTool,
                               new_file_name=simset_task.overrides()['save_channels_file'][:-4])

        # Save snapshot run results (original code - no function)
        # save_time = True  # Save the time vector only once for space-saving
        # for file_num in all_files_to_open:
        #     # Load each target file and for each block related to the file asign the corresponding data to the block
        #     if file_num < 10:  # If the file number is less than 10, then it adds 0 before the file number
        #         values = np.loadtxt(out_filename + "_0" + str(file_num) + ".out", skiprows=1)
        #     else:
        #         values = np.loadtxt(out_filename + "_" + str(file_num) + ".out", skiprows=1)
        #     for block in ScanBlocksTool:
        #         if save_time:
        #             block.snapshot_data["time"] = values[:,0]  # Retreive the time vector
        #             save_time = False  # Only for the first scanning block (to save mem space)
        #         for block_file in block.files_to_open:
        #             if file_num == block_file:
        #                 # If the block needs data from the file, then use the block's target columns for this file
        #                 for col in block.relative_cols[file_num]:
        #                     ch = col - 1 + 10 * (file_num - 1)  # Absolute output channel number
        #                     # print(block.out_vars_names[ch], ch_var_names[ch])
        #                     block.snapshot_data[block.out_vars_names[ch]] = values[:, col-1]  # Retreived data
        #
        # if save_td:
        #     file_name = results_folder + '\\' + simset_task.overrides()['save_channels_file'] + '.txt'
        #     var_names = ["time"]
        #     data = ScanBlocksTool[0].snapshot_data["time"].reshape(-1,1)  # Retreive the time vector data
        #     for block in ScanBlocksTool:
        #         for name in list(block.out_vars_names.values()):
        #             if name != "time":  # Do not save the time vector
        #                 data = np.append(data, block.snapshot_data[name].reshape(-1,1), axis=1)
        #                 var_names.append(name)
        #     np.savetxt(file_name, data, delimiter='\t', header="\t".join(var_names))

        if take_snapshot:
            # Remove the snapshot time offset
            initial_row = find_nearest(ScanBlocksTool[0].snapshot_data["time"], t_snap_internal)
            if verbose: print("Snapshot w/o time offset",ScanBlocksTool[0].snapshot_data["time"][initial_row],t_snap_internal)
            ScanBlocksTool[0].snapshot_data["time"] = ScanBlocksTool[0].snapshot_data["time"][initial_row:] - ScanBlocksTool[0].snapshot_data["time"][initial_row]

            for block in ScanBlocksTool:
                for name in list(block.out_vars_names.values()):
                    block.snapshot_data[name] = block.snapshot_data[name][initial_row:]
        print(' Unperturbed simulation results collected in', round((t.time() - t1), 2), 'seconds')

    if scantype is "AC":
        # AC-type bus scan
        # d-axis injection
        print('\n Running single frequency d-axis injection simulations')
        t1 = t.time()
        idx_selected_blocks = []
        for idx, block in enumerate(ScanBlocksTool):
            if block.type == "AC":
                block.pscad_block.parameters(V_perturb_mag=v_perturb_mag,selector=1)  # d-axis injection
                idx_selected_blocks.append(idx)
            else:
                block.pscad_block.parameters(V_perturb_mag=v_perturb_mag,selector=0)  # No injection
        simset_task.parameters(volley=num_parallel_sim, affinity_type='DISABLE_TRACING', ammunition=f_points)
        simset_task.overrides(duration=t_sim_internal, time_step=t_step, plot_step=sample_step, start_method=1,
                              timed_snapshots=0, startup_inputfile=snapshot_file + '.snp',
                              save_channels_file=output_files + '_d.out', save_channels=1)
        if run_sim: simset.run()
        print(' d-axis injection finished in', round((t.time() - t1), 2), 'seconds')
        if save_td or compute_yz:
            wait4pscad(time=1, pscad=pscad)
            t2 = t.time()
            read_and_save.multiple_s(n_sim=f_points, out_folder=out_dir, save_folder=results_folder, save=save_td,
                                     tar_files=all_files_to_open, zblocks=[ScanBlocksTool[ind] for ind in idx_selected_blocks],
                                     file_name=simset_task.overrides()['save_channels_file'][:-4])
            print(' d-axis injection results collected in', round((t.time() - t2), 2), 'seconds\n')

        # q-axis injection
        print(' Running single frequency q-axis injection simulations')
        t1 = t.time()
        for block in ScanBlocksTool:
            if block.type == "AC":
                block.pscad_block.parameters(V_perturb_mag=v_perturb_mag,selector=2)  # q-axis injection
            else:
                block.pscad_block.parameters(V_perturb_mag=v_perturb_mag,selector=0)  # No injection
        simset_task.parameters(volley=num_parallel_sim, affinity_type='DISABLE_TRACING', ammunition=f_points)
        simset_task.overrides(duration=t_sim_internal, time_step=t_step, plot_step=sample_step, start_method=1,
                              timed_snapshots=0, startup_inputfile=snapshot_file + '.snp',
                              save_channels_file=output_files + '_q.out', save_channels=1)
        if run_sim: simset.run()
        print(' q-axis injection finished in', round((t.time() - t1), 2), 'seconds')
        if save_td or compute_yz:
            wait4pscad(time=1, pscad=pscad)
            t2 = t.time()
            read_and_save.multiple_s(n_sim=f_points, out_folder=out_dir, save_folder=results_folder, save=save_td,
                                     tar_files=all_files_to_open, zblocks=[ScanBlocksTool[ind] for ind in idx_selected_blocks],
                                     file_name=simset_task.overrides()['save_channels_file'][:-4])
            print(' q-axis injection results collected in', round((t.time() - t2), 2), 'seconds\n')

    elif scantype is "DC":
        # DC-side injection
        print(' Running single frequency DC-side injection simulations')
        t1 = t.time()
        idx_selected_blocks = []
        for idx, block in enumerate(ScanBlocksTool):
            if block.type == "AC":
                block.pscad_block.parameters(V_perturb_mag=v_perturb_mag, selector=0)  # No dq-axis injection
            else:
                idx_selected_blocks.append(idx)
                block.pscad_block.parameters(V_perturb_mag=v_perturb_mag, selector=1)  # DC-side injection
        simset_task.parameters(volley=num_parallel_sim, affinity_type='DISABLE_TRACING', ammunition=f_points)
        simset_task.overrides(duration=t_sim_internal, time_step=t_step, plot_step=sample_step, start_method=1,
                              timed_snapshots=0, startup_inputfile=snapshot_file + '.snp',
                              save_channels_file=output_files + '_dc.out', save_channels=1)
        if run_sim: simset.run()
        print(' DC-side injection finished in', round((t.time() - t1), 2), 'seconds')
        if save_td or compute_yz:
            wait4pscad(time=1, pscad=pscad)
            t2 = t.time()
            read_and_save.multiple_s(n_sim=f_points, out_folder=out_dir, save_folder=results_folder, save=save_td,
                                     tar_files=all_files_to_open, zblocks=[ScanBlocksTool[ind] for ind in idx_selected_blocks],
                                     file_name=simset_task.overrides()['save_channels_file'][:-4])
            print(' DC-side injection results collected in', round((t.time() - t2), 2), 'seconds\n')

    else:
        # ACDC-type scan
        # d-axis injection
        print('\n Running single frequency d-axis injection simulations')
        t1 = t.time()
        for block in ScanBlocksTool:
            if block.type == "AC":
                block.pscad_block.parameters(V_perturb_mag=v_perturb_mag,selector=1)  # d-axis injection
            else:
                block.pscad_block.parameters(V_perturb_mag=v_perturb_mag,selector=0)  # No injection
        simset_task.parameters(volley=num_parallel_sim, affinity_type='DISABLE_TRACING', ammunition=f_points)
        simset_task.overrides(duration=t_sim_internal, time_step=t_step, plot_step=sample_step, start_method=1,
                              timed_snapshots=0, startup_inputfile=snapshot_file + '.snp',
                              save_channels_file=output_files + '_d.out', save_channels=1)
        if run_sim: simset.run()
        print(' d-axis injection finished in', round((t.time() - t1), 2), 'seconds')
        if save_td or compute_yz:
            wait4pscad(time=1, pscad=pscad)
            t2 = t.time()
            read_and_save.multiple_s(n_sim=f_points, out_folder=out_dir, save_folder=results_folder, save=save_td,
                                     tar_files=all_files_to_open, zblocks=ScanBlocksTool,
                                     file_name=simset_task.overrides()['save_channels_file'][:-4])
            print(' d-axis injection results collected in', round((t.time() - t2), 2), 'seconds\n')

        # q-axis injection
        print(' Running single frequency q-axis injection simulations')
        t1 = t.time()
        for idx, block in enumerate(ScanBlocksTool):
            if block.type == "AC":
                block.pscad_block.parameters(selector=2)  # q-axis injection
            else:
                block.pscad_block.parameters(selector=0)  # No injection
        simset_task.parameters(volley=num_parallel_sim, affinity_type='DISABLE_TRACING', ammunition=f_points)
        simset_task.overrides(duration=t_sim_internal, time_step=t_step, plot_step=sample_step, start_method=1,
                              timed_snapshots=0, startup_inputfile=snapshot_file + '.snp',
                              save_channels_file=output_files + '_q.out', save_channels=1)
        if run_sim: simset.run()
        print(' q-axis injection finished in', round((t.time() - t1), 2), 'seconds')
        if save_td or compute_yz:
            wait4pscad(time=1, pscad=pscad)
            t2 = t.time()
            read_and_save.multiple_s(n_sim=f_points, out_folder=out_dir, save_folder=results_folder, save=save_td,
                                     tar_files=all_files_to_open, zblocks=ScanBlocksTool,
                                     file_name=simset_task.overrides()['save_channels_file'][:-4])
            print(' q-axis injection results collected in', round((t.time() - t2), 2), 'seconds\n')

        # DC-side injection
        print(' Running single frequency DC-side injection simulations')
        t1 = t.time()
        for block in ScanBlocksTool:
            if block.type == "AC":
                block.pscad_block.parameters(selector=0)  # No dq-axis injection
            else:
                block.pscad_block.parameters(selector=1)  # DC-side injection
        simset_task.parameters(volley=num_parallel_sim, affinity_type='DISABLE_TRACING', ammunition=f_points)
        simset_task.overrides(duration=t_sim_internal, time_step=t_step, plot_step=sample_step, start_method=1,
                              timed_snapshots=0, startup_inputfile=snapshot_file + '.snp',
                              save_channels_file=output_files + '_dc.out', save_channels=1)
        if run_sim: simset.run()
        print(' DC-side injection finished in', round((t.time() - t1), 2), 'seconds')
        if save_td or compute_yz:
            wait4pscad(time=1, pscad=pscad)
            t2 = t.time()
            read_and_save.multiple_s(n_sim=f_points, out_folder=out_dir, save_folder=results_folder, save=save_td,
                                     tar_files=all_files_to_open, zblocks=ScanBlocksTool,
                                     file_name=simset_task.overrides()['save_channels_file'][:-4])
            print(' DC-side injection results collected in', round((t.time() - t2), 2), 'seconds\n')

    if compute_yz:
        t2 = t.time()
        print(' Computing admittances')
        # Snapshot data time-aligment
        t_0 = ScanBlocksTool[0].perturbation_data["time"][0]
        if verbose:
            print('  Simulation and snapshot initial time:', t_0, ScanBlocksTool[0].snapshot_data["time"][0])
            print("  Shape snapshot:",ScanBlocksTool[0].snapshot_data["time"].shape,"perturb:",ScanBlocksTool[0].perturbation_data["time"].shape)

        # Shift the snapshot data by one time-step if needed
        if t_0 != 0.0 and round(ScanBlocksTool[0].snapshot_data["time"][0],10) != round(t_0,10):
            print("  Shifting snapshot")
            ScanBlocksTool[0].snapshot_data["time"] = ScanBlocksTool[0].snapshot_data["time"][1:]
            for block in ScanBlocksTool:
                for name in list(block.out_vars_names.values()):
                    block.snapshot_data[name] = block.snapshot_data[name][1:]  # Time has been modified already

        # Time vector sampling time and start FFT index
        dt = ScanBlocksTool[0].perturbation_data["time"]  # Dummy variable time vector
        dt = np.mean([dt[i + 1] - dt[i] for i in range(min(len(dt), 100))])  # Sampling time [s]
        if verbose: print(" Sampling time: ",dt)
        start_idx = find_nearest(ScanBlocksTool[0].perturbation_data["time"], start_fft)
        if verbose: print(" FFT time: ", start_fft, "Time vector value: ", ScanBlocksTool[0].perturbation_data["time"][start_idx])

        Ytopology_scan = np.copy(Ytopology)  # Make a copy to modify as the different scans are being done
        for idx, name in enumerate(block_names_Y):
            if sum(Ytopology_scan[idx,:]) == 1:
                # Only for point-to-point connections
                nz = np.nonzero(Ytopology_scan[idx,:])[0]  # Find the index for the other scan block
                # The loop looks for the point to point blocks, but it can be improved (see passives filt)
                for block in ScanBlocksTool:
                    if block.name == name[:-2]:
                        block_type0 = block.type  # The current block
                        block0 = block
                        if verbose: print("P2P block name", block.name)
                    if block.name == block_names_Y[int(nz)][:-2]:
                        block_type1 = block.type  # The other scan block
                        block1 = block
                        if verbose: print("P2P block name", block.name)
                if block_type0 != block_type1:
                    t1 = t.time()
                    # If one block is AC and the other DC: AC/DC converter scan
                    # sides = [str((idx % 2)+1),str((int(nz) % 2)+1)]  # Sides to retreive the data from each block
                    sides = [name[-1:], block_names_Y[int(nz)][-1:]]
                    if verbose:
                        print("  Name of the blocks: ",block0.name, block1.name)
                        print("  Side of the blocks: ",sides)
                    yz_computation.admittance(f_base=f_base, frequencies=freq, fft_periods=fft_periods,dt=dt,
                                              start_idx=start_idx,
                                              zblocks=[block0,block1], sides=sides, scantype="ACDC",
                                              results_folder=results_folder, results_name=output_files)
                    # Update the scan matrix to indicate that no scan from and to these two ports is pending
                    Ytopology_scan[nz,idx] = 0
                    Ytopology_scan[idx,nz] = 0
                    print('  Admittance between',block0.name+sides[0],'and',block1.name+sides[1],'computed in',round((t.time() - t1), 2),'seconds')
                else:
                    # Both are AC or DC scan blocks
                    if block0.name == block1.name:
                        # It is the same name = one-side scan
                        t1 = t.time()
                        yz_computation.admittance(f_base=f_base, frequencies=freq, fft_periods=fft_periods, dt=dt,
                                                  start_idx=start_idx,
                                                  zblocks=block0, sides=name[-1:], scantype=block0.type,
                                                  results_folder=results_folder, results_name=output_files)
                        # Update the scan matrix to indicate that no scan from and to these two ports is pending
                        Ytopology_scan[nz, idx] = 0
                        Ytopology_scan[idx, nz] = 0
                        print('  Admittance at',name,'computed in', round((t.time() - t1), 2), 'seconds \n')

        print(' Admittance computation finished in ', round((t.time() - t2), 2), 'seconds')

    """ Perform the simulations for the scan of the passive networks based on the topology information """
    sim_select = {"_d": 1, "_q": 2}  # Dict containing the AC simulation type based on the name ending
    # Disable all perturbations
    for block in ScanBlocksTool:
        block.pscad_block.parameters(V_perturb_mag=v_perturb_mag, selector=0)
    # Shorter simulation times can be set (EMT transients are faster) but for simplicity the former settings are used
    # The PSCAD project folder (i.e. project_name.gf46) can be cleared here to decrease memory usage in large projects

    if (len(passive_networks_scans) != 0) and compute_yz:
        # Passive network scan
        print("\nScan of AC and/or DC networks")
        t1 = t.time()
        for run in range(1, max_runs + 1):
            t2 = t.time()
            idx_selected_blocks = []  # This list changes every run so only the necessary blocks store the PSCAD results
            # Configure every PSCAD block involved in the network scan according to the remaining scans
            for network_scan in passive_networks_scans:
                if len(network_scan.remaining_scans) > 0:
                    # There are still scans to perform: update the selected Z-blocks and iterate over the blocks & names
                    for block_idx in list(network_scan.blocks_idx.values()): idx_selected_blocks.append(block_idx)
                    if network_scan.scan_type == "AC":
                        # If the scan is AC then set up the d/q injection or no injection
                        for name in network_scan.names:
                            if network_scan.remaining_scans[0][:-2] == name:
                                # Enable the perturbation at this block
                                ScanBlocksTool[network_scan.blocks_idx[name]].pscad_block.parameters(
                                    selector=sim_select[network_scan.remaining_scans[0][-2:]])
                                if verbose: print("  AC perturbation at", name)
                            else:
                                # Disable the rest of the blocks from this network
                                ScanBlocksTool[network_scan.blocks_idx[name]].pscad_block.parameters(selector=0)
                    else:
                        # If the scan is DC then set up the dc injection or no injection
                        for name in network_scan.names:
                            if network_scan.remaining_scans[0][:-3] == name:
                                if verbose: print("  DC perturbation at",name)
                                # Enable the perturbation at this block
                                ScanBlocksTool[network_scan.blocks_idx[name]].pscad_block.parameters(selector=1)
                            else:
                                # Disable the rest of the blocks from this network
                                ScanBlocksTool[network_scan.blocks_idx[name]].pscad_block.parameters(selector=0)

                    # Update the pending scans
                    network_scan.updateScans(network_scan.remaining_scans[0])

            # Perform the simulations and label the output data by using the "run" number
            simset_task.parameters(volley=num_parallel_sim, affinity_type='DISABLE_TRACING',ammunition=f_points)
            simset_task.overrides(duration=t_sim_internal, time_step=t_step, plot_step=sample_step,start_method=1,
                                  timed_snapshots=0, startup_inputfile=snapshot_file + '.snp',
                                  save_channels_file=output_files + "_" + str(run) + '.out', save_channels=1)
            if run_sim: simset.run()
            print(' Run',str(run)+'/'+str(max_runs),'finished in', round((t.time() - t2), 2), 'seconds')
            if save_td or compute_yz:
                wait4pscad(time=1, pscad=pscad)
                t2 = t.time()
                read_and_save.multiple_s(n_sim=f_points, out_folder=out_dir, save_folder=results_folder, save=save_td,
                                         tar_files=all_files_to_open,
                                         zblocks=[ScanBlocksTool[ind] for ind in idx_selected_blocks],
                                         file_name=simset_task.overrides()['save_channels_file'][:-4])
                print(' Results collected in', round((t.time() - t2), 2), 'seconds\n')

        # Compute the admittance for each network individually
        for network_scan in passive_networks_scans:
            t2 = t.time()
            idx_selected_blocks = []
            sides_selected_blocks = []
            for name in network_scan.names:
                idx_selected_blocks.append(network_scan.blocks_idx[name])
                sides_selected_blocks.append(name[-1])
            print(" Computing admittance for the network:",' '.join(network_scan.names))
            if verbose:
                aux = [ScanBlocksTool[ind].name+"-"+sides_selected_blocks[idx] for idx, ind in enumerate(idx_selected_blocks)]
                print("  PSCAD blocks, name with side:",' '.join(aux))
            yz_computation.admittance(f_base=f_base, frequencies=freq, fft_periods=fft_periods, dt=dt,
                                      start_idx=start_idx, scantype="Network", results_folder=results_folder,
                                      zblocks=[ScanBlocksTool[ind] for ind in idx_selected_blocks],
                                      sides=[side for side in sides_selected_blocks], network=network_scan,
                                      results_name=output_files)

            print('  Admittance matrix involving'," ".join(network_scan.names),'computed in',round((t.time() - t2), 2), 'seconds')

    # project.save()  # Save the project changes
    print(' Quitting PSCAD')
    wait4pscad(time=1, pscad=pscad)
    pscad.quit()  # Quit PSCAD
    print('\nTotal execution time', round((t.time() - t0) / 60, 2), 'minutes\n')

def wait4pscad(time=1, pscad=None):
    busy = pscad.is_busy()
    while busy:
        t.sleep(time)  # Wait a bit
        busy = pscad.is_busy()

def find_nearest(array, value):  # Efficient function to find the nearest value to a given one and its position
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or np.abs(value - array[idx - 1]) < np.abs(value - array[idx])):
        return idx - 1
    else:
        return idx


def read_one_line(file_path, nline):
    with open(file_path, 'r') as file:
        for line_number, line in enumerate(file):
            if line_number > nline + 1:  # Offset by empty header line
                break
            if line_number == nline + 1: selected_data = line.split()
    return selected_data


frequency_sweep_ACDC_full.__doc__ = """
Author: Francisco Javier Cifuentes García
V0.1 [03/08/2022]
PSCAD automation
The values of the frequency list are multiples of the base frequency.
The simulation configuration is set to perform one snapshot and then the frequency sweeps.
The function accepts several input arguments to customize the frequency sweep:

Required
        t_snap          Time when the snapshot is taken [s].
        take_snapshot		Bool: Does the user want to take a snapshot? Default = False. A previous snapshot can still be used for the Y computation.
                                The snapshot simulation runs for t_snap_internal + t_sim_internal so as to save the steady-state unperturbed waveforms.
        t_sim			Duration of each frequency injection simulation [s].
        t_step			Simulation time step [us].
        sample_step		Sample time of the output channels [us].
        v_perturb_mag	Voltage perturbation in per unit [pu]
        freq			Frequencies to perform the injections [Hz]. Alternatively, the user can provide info to compute the list.
        dt_injections   Time to reach steady-state after a perturbation [s]. Alternatively, start_fft can be provided
      
Optional
        f_points			Number of frequency perturbation points.
        f_base		 	Base frequency [Hz].
        f_min			Start frequency [Hz].
        f_max			End frequency [Hz].
        fft_periods 		Number of periods used to compute the FFT. Default = 1.
        start_fft		Time for the DUT to reach the new steady-state (injections) [s] . 
        freq_text_file	 	File name with the perturbation frequencies. Default = 'frequencies.txt'.
        project_name	 	Name of the project. Default 'Impedance_testing_single_frequ'.
        workspace_name	 	Name of the workspace. Default 'STATCOM_Worksace'.
        fortran_ext	 	Fortran extension. Default r'.gf46'.
        num_parallel_sim 	Number of parallel simulations. Default 8.
	component_parameters	List of component parameter values to be modified in PSCAD. E.g. [ParamVal_1, ParamVal_2, 3.0148, 0,0,0]
	working_dir	 	Working directory (in case the python file is not in the same folder as the PSCAD project).
        snapshot_file		Name of the snapshot.
        output_files	 	Name of the output files
        results_folder	 	Absolute path where the formated results will be stored. (If not specified, they are not saved)
                         	If specified, it can also perform the analysis of the results: admitance / impedance measurement.
        compute_yz	 	Compute the admittance and save the results. If no results folder is specified then it saves the data in working_dir.
        save_td  		Bool: If set to True, the several files of time domain data are saved into more compact .txt files for each independent
                                perturbation. The format is [time Vd(f1) Vq(f1) Id(f1) Iq(f1) ... Vd(f_max) Vq(f_max) Id(f_max) Iq(f_max)].
"""
