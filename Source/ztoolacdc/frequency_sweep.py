""" This script contains the functions to perform different frequency scans """
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
__all__ = ['frequency_sweep', 'frequency_sweep_TF']

import time as t  # Relative time
from os import listdir, chdir, getcwd, path, makedirs
import numpy as np  # Numerical python functions
import matplotlib.pyplot as plt  # Plot library
import pickle  # Handle for interactive plots
from mhi.pscad import launch  # PSCAD automation functions
from . import create_freq, read_and_save, yz_computation

# Output channels for each scanning block: blockid is assummed to be the first output when retreiving the data /!\
AC_scan_variables = ['blockid', 'IDUTacA1:1', 'IDUTacA1:2','VDUTac:1', 'VDUTac:2', 'IDUTacA2:1', 'IDUTacA2:2','theta']
DC_scan_variables = ['blockid', 'IDUTdcA1', 'IDUTdcA2','VDUTdc']

class Network:
    def __init__(self, name_blocks_involved, scan_type, adj_matrix):
        self.names = name_blocks_involved  # List of names including the sides
        self.names_wo_sides = [names[:-2] for names in self.names] # List of names without the sides
        self.adj_matrix = adj_matrix  # The adjacent matrix = zeros correspond to y = 0 (disconnection)
        if scan_type == "AC":
            self.scan_type = scan_type  # AC scan
            perturbations = ["_d","_q"]
            self.all_scans = [name[:-2] + perturbation_type for name in self.names for perturbation_type in perturbations]  # Block name (without sides) + perturbation type
        elif scan_type == "DC":
            self.scan_type = scan_type  # DC scan
            perturbations = ["_dc"]
            self.all_scans = [name[:-2] + perturbation_type for name in self.names for perturbation_type in perturbations]  # Block name (without sides) + perturbation type
        else:
            self.all_scans = []
            for idx, name in enumerate(self.names):
                if scan_type[idx] == "AC":
                    self.all_scans.append(name[:-2] + "_d")
                    self.all_scans.append(name[:-2] + "_q")
                else:
                    self.all_scans.append(name[:-2] + "_dc")        
            self.scan_type = "ACDC"
        
        self.runs = len(self.all_scans)  # Number of needed runs for the network scan
        self.blocks_idx = None  # Dict key: self.names, pointing at the key of the associated blocks in ScanBlocksTool
        self.scan_per_run = {}  # Dictionary with the block+side+perturbation for the runs (keys) in runs_list
        self.runs_list = []  # Keys to the scan_per_run dict indicating in which runs the scans of this subsystem are performed
        self.enforce = True  # Enforce network connectivity based on the topology file when computing and plotting the admittance

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
        
        if "AC" in "".join(self.pscad_block.defn_name):
            self.type = "AC"
            self.var_names = ['VDUTac', 'IDUTacA1', 'IDUTacA2', 'theta']  # Root of the variable names
            self.group = "ACscan"
            self.area = None  # Corresponding AC area identificator
            self.theta = 0.0  # Steady-state voltage angle
            self.is_ref = 0  # Is the angle of this component the reference for its area?
        elif "DC" in "".join(self.pscad_block.defn_name):
            self.var_names = ['VDUTdc', 'IDUTdcA1', 'IDUTdcA2']  # Root part of the variable names
            self.type = "DC"
            self.group = "DCscanPM"
            self.area = None  # Corresponding DC area identificator (not used yet)
        else:
            self.var_names = ['output_TF', 'input_TF']  # Root part of the variable names
            self.type = "TF"
            self.group = "TFscan"

def visualize_graph(G, node_names, save_dir, file_name):
    # Very basic function that plots the graph given by the adjancent matrix (also with the node names)
    # Better results and flexibility can be obtained using NetworkX
    num_nodes = G.shape[0]
    # Get node positions using circular layout
    pos = {i: (np.cos(2 * np.pi * i / num_nodes), np.sin(2 * np.pi * i / num_nodes)) for i in range(num_nodes)}

    # Draw nodes
    plt.figure()
    plt.scatter(*zip(*pos.values()), color='skyblue', s=700)

    # Draw edges
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            if G[i][j] == 1:
                plt.plot([pos[i][0], pos[j][0]], [pos[i][1], pos[j][1]], color='gray')

    # Draw additional nodes for self-loops
    for i in range(num_nodes):
        if G[i][i] == 1:
            plt.scatter(pos[i][0], pos[i][1], color='red', s=200, zorder=10)

    # Draw node labels above the nodes
    for i, name in enumerate(node_names):
        plt.text(pos[i][0], pos[i][1] + 0.1, name, fontsize=12, ha='center', va='bottom')

    plt.axis('off')
    # plt.show()
    plt.savefig(save_dir + r'\\' + file_name + '_network_visualization.pdf',format="pdf", bbox_inches="tight")

    with open(save_dir + r'\\' + file_name + "_network_visualization.pickle", 'wb') as f:
        pickle.dump(plt.gcf(), f)  # Save the file for allowing plot interaction

def create_scan_schedule(passive_networks_scans):
    # Avoid linearly-dependent perturbations of subsystems sharing scan blocks by defining a scan run schedule
    all_scans = list(set(scan for network_aux in passive_networks_scans for scan in network_aux.all_scans))  # Name of the block with perturbation (no sides) 

    # Greedy algorithm: check all scans one by one and add it to the first run without an existing scan of the same subsystem
    scheduled_scans = [[] for i in range(len(all_scans))]  # Worst case: no scans can be done simultaneously
    involved_networks = [[] for i in range(len(all_scans))]  # Involved networks in each scan run
    
    for scan in all_scans:
        scan_scheduled = False

        # Find the networks with the current scan block
        networks_with_current_block = []
        for network in passive_networks_scans:
            # If the name of the block (no sides, no perturbation) is in the network, then add it to the networks_with_current_block list
            if "_".join(scan.split("_")[:-1]) in network.names_wo_sides:
                networks_with_current_block.append(network)

        for run in range(len(scheduled_scans)):
            if not scan_scheduled and not any(element in involved_networks[run] for element in networks_with_current_block):
                # Schedule the scan/block to this run if it does not belong to a system with blocks in the current run already
                scan_scheduled = True
                scheduled_scans[run].append(scan)
                # Update the involved networks' runs_list and scan_per_run
                for network in networks_with_current_block:
                    involved_networks[run].append(network)  # List to avoid overlapping (networks sharing perturbation blocks)
                    if scan in network.all_scans:
                        # Use the name without side and perturbation to access the corresponding list of names with sides
                        index_of_name = network.names_wo_sides.index("_".join(scan.split("_")[:-1]))
                        # Add the name with sides and perturbation to the scan per run list
                        network.runs_list.append(run+1)
                        network.scan_per_run[run+1] = network.names[index_of_name] + "_" + scan.split("_")[-1]

def frequency_sweep(t_snap=None, t_sim=None, t_step=None, sample_step=None, v_perturb_mag=None,freq=None, f_points=None,
                    f_base=None, f_min=None, f_max=None, working_dir=None, multi_freq_scan=False,
                    snapshot_file=None, dedicated_SS_sim=False, take_snapshot=True, dt_injections=None, scan_actives = True,
                    topology=None, project_name=None, workspace_name=None,fortran_ext=r'.gf46', num_parallel_sim=8,
                    scan_passives=True, edge_dq_sym=False, edge_sym=False, component_parameters=None,
                    results_folder=None, output_files='Perturbation', compute_yz=True, save_td=False,
                    fft_periods=1, start_fft=None, pscad_plot=0, show_powerflow=False, visualize_network=False,
                    run_sim=True,verbose=False):

    # Debugging control: run_sim enables or dissables all PSCAD simulations; verbose enables or disables script run info
    """ --- Input data handling --- """
    # CHECK the following if: it does not work... "or" operator gives a bool as output not None
    if (t_snap or t_step or start_fft or v_perturb_mag or project_name or workspace_name or ((f_points or f_base or f_max or f_min) and freq)) is None:
        print('One or more required arguments are missing!! \n Check the function documentation by typing help(frequency_sweep) \n')
        return

    # Create frequency list if it is not provided
    if freq is None: freq = create_freq.loglist(f_min=f_min, f_max=f_max, f_points=f_points, f_base=f_base)

    # If the sample time is not provided, it is set to half of the minimum required value (multiple of step_time)
    if sample_step is None:
        if f_max is None: f_max = max(freq)
        sample_step = round(t_step * np.floor((1e6 * 0.5 * 0.5 / f_max + t_step / 2) / t_step), 3)  # [us]
    if t_sim is None: t_sim = start_fft + fft_periods / f_base
    if dt_injections is None: dt_injections = start_fft  # dt_injections
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

    # Create the folder to store the results if it does not exist or use the current one
    if (results_folder is not None) and (not path.exists(results_folder)): makedirs(results_folder)
    if (results_folder is None) and (save_td or compute_yz): results_folder = working_dir

    # Folder for the scan options
    if not path.exists(working_dir+"Scan_options"): makedirs(working_dir+"Scan_options")

    # Frequency vector to PSCAD XY table: single frequency perturbations
    freq_text_file='frequencies.txt'
    f_points = len(freq)
    with open(working_dir+"Scan_options\\"+freq_text_file, 'w') as f:  # Create the .txt file
        f.write('! This file stores the frequency sweep values in Hz \n')  # Write header
        for j in range(f_points): f.write(str(j + 1) + '\t' + str(freq[j]) + '\n')  # Write values
        f.write('ENDFILE:')  # Write end of the file
    f.close()

    # Frequency vector to several PSCAD XY tables: multi-frequency perturbations
    f_points_per_file = int(np.ceil(f_points / num_parallel_sim))  
    N_files = 8  # 8 frequencies at the same time by default
    freq_multi = np.pad(freq, (0, f_points_per_file*N_files - f_points), 'constant') # Add zeros if not enough frequencies
    freq_multi = freq_multi.reshape(N_files, f_points_per_file)
    for freq_file_num in range(N_files):
        with open(working_dir+"Scan_options\\"+freq_text_file[:-4]+str(freq_file_num+1)+".txt", 'w') as f:  # Create the .txt file
            f.write('! This file stores the multi-frequency sweep values in Hz \n')  # Write header
            for j in range(f_points_per_file):
                f.write(str(j + 1) + '\t' + str(freq_multi[freq_file_num][j]) + '\n')  # Write values
            f.write('ENDFILE:')  # Write end of the file is apparently not needed for this component
        f.close()

    """ --- Main program --- """
    t0 = t.time()  # Initial time
    print('Launching PSCAD')
    pscad = launch(minimize=True)
    wait4pscad(time=2, pscad=pscad)
    t.sleep(5)  # Wait a bit more just in case PSCAD is still loading
    if not pscad.licensed():
        certificates = pscad.get_available_certificates()  # Retreive PSCAD license certificate
        keys = list(certificates.keys())
        pscad.get_certificate(certificates[keys[0]])  # Get the first available certificate
    pscad.load(working_dir + workspace_name + ".pswx")  # Load workspace
    project = pscad.project(project_name)  # Where the model is defined
    main = project.canvas('Main')  # Main of the project
    project.focus()  # Switch PSCAD’s focus to the project

    print(' Setting parameters and configuration')
    # Get components and set global parameters for all simulations
    if topology is None:
        # If a single block is in the canvas, only that one needs to be scanned
        blocks = main.find_first("Z_tool:ACscan")  # This assumes a single scan block in the main canvas
        if blocks is None: blocks = main.find_first("Z_tool:DCscanPM")  # If it did not find a ACscan look for PM DCscan
        if blocks is None: blocks = main.find_first("Z_tool:DCscan")
        scanid = [blocks.parameters()['Name']]  # Retrieve the scan block name for identification
        Ytopology = np.identity(2)
        block_names_Y = [scanid[0]+"-1",scanid[0]+"-2"]
        print("Warning: no topology has been specified, this assumes a single scan block is available |",scanid[0])
    else:
        # Read the topology matrix
        Ytopology = np.loadtxt(topology, skiprows=1, comments=["#", "%", "!"])
        # len(ScanBlocksTool)*2 by len(ScanBlocksTool)*2 # nameA-1 nameA-2 nameB-1 nameB-2 ... x nameA-1 nameA-2 nameB-1 ...
        # 0 means no interconnection, 1 means connection between the edges: diagonals are single-sided / shunt
        with open(topology, 'r') as f:
            block_names_Y = f.readline().strip('\n').split()
        scanid = [name[:-2] for name in block_names_Y[::2]]  # List of PSCAD block names / identifiers
    if verbose:
        print("Detailed block names:",block_names_Y)
        print("PSCAD block names:", scanid)

    ScanBlocksAC = []
    ScanBlocksDC = []
    for identification in scanid:
        blocks = main.find_all(Name=identification)
        # blocks_tool = [block for block in blocks if "Z_tool" in block.defn_name[0]]  # Filter only Z-tool components
        blocks_tool = [block for block in blocks if "Z_tool" in "".join(block.defn_name)]  # Filter only Z-tool components
        if "DC" in "".join(blocks_tool[0].defn_name):
            ScanBlocksDC.append(blocks_tool[0])
        else:
            ScanBlocksAC.append(blocks_tool[0])

    ScanBlocksAC_names = [block.parameters()['Name'] for block in ScanBlocksAC]
    ScanBlocksDC_names = [block.parameters()['Name'] for block in ScanBlocksDC]

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
    if verbose:
        print("AC scan blocks", ScanBlocksAC)
        print("DC scan blocks", ScanBlocksDC)
    print(' Type of scan:', scantype)

    # Create the undirected graph - adjacent matrix but diagonals can be 1
    g = Graph(len(Ytopology))
    for row, name in enumerate(block_names_Y):
        for col, edge in enumerate(Ytopology[row]):
            if int(edge) == 1: g.addEdge(row, col)
    if verbose: print("Block names:"," | ".join(block_names_Y))
    # Obtain the connected components of the graph
    cc = g.connectedComponents()  # List of lists with blocks # involved in each scan

    # Visualization (in progress)
    if visualize_network:
        CC_for_plot = np.zeros((Ytopology.shape[0] // 2, Ytopology.shape[0] // 2))  # Initialize new matrix
        for i in range(0, Ytopology.shape[0], 2):
            for j in range(0, Ytopology.shape[0], 2):
                CC_for_plot[i//2, j//2] = Ytopology[i, j] + Ytopology[i, j+1] + Ytopology[i+1, j] + Ytopology[i+1, j+1]
        if verbose: print(" Connected network plot\n",CC_for_plot,"\n",[blocks_name[:-2] for k, blocks_name in enumerate(block_names_Y) if k % 2])
        visualize_graph(CC_for_plot,[blocks_name[:-2] for k, blocks_name in enumerate(block_names_Y) if k % 2],results_folder,output_files)

    # Extract the interconnected blocks for network scan (remove shunts a.k.a. unconnected vertices)
    scans_network = [c for c in cc if len(c) != 1]  # List of lists with blocks # in each scan with more than 1 block
    # After this for loop, the list is filtered & enhanced by a network scans object: num of runs, names, adj matrix...
    passive_networks_scans = []  # This variable stores said lists
    acdc_converters_blocks = []  # Stores the block names of the AC-side of AC/DC converters
    for idx, net in enumerate(scans_network):
        network_names = [block_names_Y[element] for element in net]
        if len(net) > 2:
            # len(net) > 2 changed to len(net) >= 2 to allow for AC/DC converters' scan
            # Multiterminal network
            # if network_names[0][:-2] in ScanBlocksAC_names:
            if all(elem in ScanBlocksAC_names for elem in [block_name[:-2] for block_name in network_names]):
                passive_networks_scans.append(Network(network_names, "AC", Ytopology[net, :][:, net]))
                if verbose: print("   AC network scan involving", network_names)
            elif all(elem in ScanBlocksDC_names for elem in [block_name[:-2] for block_name in network_names]):
                passive_networks_scans.append(Network(network_names, "DC", Ytopology[net, :][:, net]))
                if verbose: print("   DC network scan involving", network_names)
            else:
                scan_type = []
                for block_name in network_names:
                    if block_name[:-2] in ScanBlocksAC_names:
                        scan_type.append("AC")
                    else:
                        scan_type.append("DC")
                passive_networks_scans.append(Network(network_names, scan_type, Ytopology[net, :][:, net]))
                if verbose: print("   AC/DC network scan involving", network_names)

        else:
            # Point to point: it needs to check that it is not an AC/DC converter
            types = []
            names = []
            for name in network_names:
                if name[:-2] in ScanBlocksAC_names:
                    types.append("AC")
                    names.append(name)
                else:
                    types.append("DC")
                    names.append(name)
            if types[0] == types[1]:
                # If the block type are the same, then they are interconnecting an AC or DC network
                passive_networks_scans.append(Network(network_names, types[0], Ytopology[net, :][:, net]))
                if verbose: print("   " + types[0] + " network scan involving", network_names)
            else:
                # AC / DC converter: decouples AC grids in terms of reference frames / angles
                acdc_converters_blocks.append(names[0])
                acdc_converters_blocks.append(names[1])

    acdc_converters_blocks = list(set(acdc_converters_blocks))  # Get rid of repetitions
    # Create a new undirected graph for AC areas identification
    if verbose: print("AC/DC converters at buses:", acdc_converters_blocks)
    g_new = Graph(len(Ytopology))
    for row, name in enumerate(block_names_Y):
        for col, edge in enumerate(Ytopology[row]):
            if col == row or block_names_Y[col][:-2] == name[:-2]:
                g_new.addEdge(row, col)  # Same scan block
                # if verbose: print("Redundant OR test check",col == row, block_names_Y[col][:-2] == name[:-2])
            if int(edge) == 1 and name not in acdc_converters_blocks: g_new.addEdge(row, col)  # Not AC/DC
    # Obtain the connected components of the new graph
    cc_new = g_new.connectedComponents()  # List of lists with blocks in the same area
    if verbose: print("CC for area id: ",cc_new)

    # TODO Asign the unique block id based on the names in the topology file, i.e. first block is 1, second is 2, etc
    ScanBlocks.sort(key=lambda x: x.parameters()['Name'][:-3], reverse=False)  # Sort the blocks by their "bus" number
    ScanBlocks_id = [i for i in range(1, len(ScanBlocksAC) + len(ScanBlocksDC) + 1)]  # Unique scan block_id signals
    # Create a list with the active scan block objects containing rich information about each block
    ScanBlocksTool = []
    ScanBlocksTool_names = []  # Name of the blocks as they appear in ScanBlocksTool to avoid excesive iterations
    for idx, block in enumerate(ScanBlocks):
        # Set snapshot parameters and block ID in the scan blocks
        block.parameters(Tdecoupling=t_snap, T_inj=t_snap_internal, selector=0, block_id=ScanBlocks_id[idx])
        ScanBlocksTool.append(Scanblock(block, block.parameters()['Name'], int(block.parameters()['block_id'])))
        # if verbose: print(" Scan block type ",block.defn_name[1])
        for area_id, blocks in enumerate(cc_new):
            if ScanBlocksTool[-1].name in [block_names_Y[num][:-2] for num in blocks] and ScanBlocksTool[-1].type == "AC":
                ScanBlocksTool[-1].area = area_id  # Only for AC areas (so far)
        ScanBlocksTool[idx].perturbation_data = {i: {} for i in range(f_points)}  # Dict of dicts
        ScanBlocksTool_names.append(ScanBlocksTool[idx].name)
        if verbose:
            if ScanBlocksTool[idx].type == "AC":
                print("AC scan block",ScanBlocksTool[idx].name,'with block_id',int(block.parameters()['block_id']),"at area",ScanBlocksTool[idx].area)
            else:
                print("DC scan block",ScanBlocksTool[idx].name,'with block_id',int(block.parameters()['block_id']))
        # The following lines identify the ACDC scan points
        # ScanBlocks_type.append(ScanBlocksTool[idx].type)  # List with block's scan type
        # Alternative: if "DC" in block.defn_name[1]: ScanBlocks_type.append("DC") # No need for AC or DC in block name
        # if "AC" in block.defn_name[1]: ScanBlocks_type.append("AC")
        # if idx > 0:
        #     if block.parameters()['Name'][:-3] == ScanBlocks[idx - 1].parameters()['Name'][:-3]:
        #         # If two buses have the same number, then it is an ACDC bus
        #         ScanBlocks_type[idx] = "ACDC"
        #         ScanBlocks_type[idx - 1] = "ACDC"
        # Retrieve the indexes of ScanBlocksTool for each network based on ScanBlocksTool_names

    for net in passive_networks_scans:
        net.blocks_idx = {name: ScanBlocksTool_names.index(name[:-2]) for name in net.names}  # Assumes no equal names
        if edge_dq_sym and net.scan_type == "AC":
            # If the matrices are assummed dq-symmetric, i.e. Ydq = -Yqd and Ydd = Yqq, then the number of runs is halved
            net.all_scans = net.all_scans[::2]  # Only d-axis perturbations: int(net.runs / 2)
            net.runs = len(net.all_scans)
        if edge_sym:
            # If the matrices are assummed symetric, i.e. Y = transpose(Y), then the number of runs is np.ceil((N+1)/2)
            net.runs = np.ceil((net.runs + 1)/2)  # TODO Perturbations and matrix computations for Y symmetric

    # Set variable component parameters
    if component_parameters is not None:
        for i in range(len(component_parameters)):
            parameter_component = main.find('master:const', component_parameters[i][0])
            parameter_component.parameters(Value=component_parameters[i][1])

    # Disable output channel components (all but Z-tool's scaning blocks) and animated displays
    if not pscad_plot:
        scan_vars = ['blockid','VDUTac','IDUTacA1','IDUTacA2','VDUTdc','IDUTdcA1','IDUTdcA2','theta']  # Target outputs
        all_pgb = project.find_all("master:pgb")  # Find all output channels in the project
        for pgb in all_pgb:
            if not (pgb.parameters()['Name'] in scan_vars):  pgb.disable()  # Disable the non-selected outputs
        all_multimeters = project.find_all("master:multimeter")  # Find all multimeters in the project
        for multimeter in all_multimeters: multimeter.parameters(Dis=0)  # Animated display is disabled (0)

    # Set simulation-specific parameters
    if 'Perturbation' not in pscad.simulation_sets():
        simset = pscad.create_simulation_set('Perturbation')
        simset.add_tasks(project_name)
        if verbose: print(' A simulation set has been created')
    else:
        simset = pscad.simulation_set('Perturbation')
    simset_task = simset.tasks()[0]  # The task is extracted

    if take_snapshot:  # It runs the snapshots
        print(' Running snapshot')
        t1 = t.time()
        if snapshot_file is None: snapshot_file = 'Snapshot'
        simset_task.parameters(volley=1, affinity=1, ammunition=1)
        if dedicated_SS_sim:
            # Run first just to take a snapshot without saving the data
            # This because with dll files the states can be more efficiently saved only at the end of the (snapshot) run
            simset_task.overrides(duration=t_snap_internal, time_step=t_step, plot_step=sample_step,
                                  start_method=0, timed_snapshots=1, snapshot_file=snapshot_file + '.snp',
                                  snap_time=t_snap_internal, save_channels_file=snapshot_file + '.out', save_channels=0)
            if run_sim: simset.run()
            # Run again to record the steady-state waveforms (no perturbations)
            simset_task.overrides(duration=t_sim_internal, time_step=t_step, plot_step=sample_step, start_method=1,
                                  timed_snapshots=0, startup_inputfile=snapshot_file + '.snp',
                                  save_channels_file=snapshot_file + '.out', save_channels=1)
            if run_sim: simset.run()
        else:
            # Run the snapshot and a bit longer to record the steady-state waveforms
            simset_task.overrides(duration=t_snap_internal + t_sim_internal, time_step=t_step, plot_step=sample_step,
                                  start_method=0, timed_snapshots=1, snapshot_file=snapshot_file + '.snp',
                                  snap_time=t_snap_internal, save_channels_file=snapshot_file + '.out', save_channels=1)
            if run_sim: simset.run()
        print(' Snapshot completed in', round((t.time() - t1), 2), 'seconds')
    elif snapshot_file is not None:
        # It performs the unperturbed simulation starting from the given snapshot (not fully tested yet)
        t1 = t.time()
        simset_task.parameters(volley=1, affinity=1, ammunition=1)  # affinity_type = "DISABLE_TRACING" disables the plotting
        simset_task.overrides(duration=t_sim_internal, time_step=t_step, plot_step=sample_step, start_method=1,
                              timed_snapshots=0, startup_inputfile=snapshot_file + '.snp',
                              save_channels_file=snapshot_file+'.out', save_channels=1)
        if dedicated_SS_sim:
            print(' Running steady-state simulation from a given snapshot')
            if run_sim: simset.run()
            print(' Steady-state simulation completed in', round((t.time() - t1), 2), 'seconds')
        # Otherwise, just read the outputs from previous scan

    # Identify the variables to be retrieved, their output channels and associated Z-tool scan blocks
    wait4pscad(time=1, pscad=pscad)
    t1 = t.time()
    ch_var_names = dict()
    out_filename = out_dir + "\\" + snapshot_file  # Snapshot output filename
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
    if verbose: print("Output channel variable names \n","\n".join(names))
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
                    print(" Block ", block.name, "\n \t Target output channels ", block.out_vars_ch)
                    print(" \t Output names: ", [block.out_vars_names[ch] for ch in block.out_vars_ch])

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
                    if verbose: print("\t Block ", block.name, " output file:", file_num,", columns:",cols)

    all_files_to_open = list(set(all_files_to_open))  # Get rid of repetitions

    # Save snapshot and steady-state run results
    read_and_save.single_s(out_files=out_filename, save_folder=results_folder,
                           save=save_td, files=all_files_to_open, zblocks=ScanBlocksTool,
                           new_file_name=simset_task.overrides()['save_channels_file'][:-4])

    # Remove the snapshot time offset and show power flow
    if dedicated_SS_sim:
        # snapshot_data contains already the steady-state simulation data starting from a dedicated snapshot file
        initial_row = 0
    else:
        # Snapshot data contains both the cold-start data and the steady-state waveforms
        initial_row = find_nearest(ScanBlocksTool[0].snapshot_data["time"], t_snap_internal)
        if verbose: print(" Snapshot w/o time offset", ScanBlocksTool[0].snapshot_data["time"][initial_row], t_snap_internal)
        # ScanBlocksTool[0].snapshot_data["time"] = ScanBlocksTool[0].snapshot_data["time"][initial_row:] - ScanBlocksTool[0].snapshot_data["time"][initial_row]

    """ ---------- Save power flow and remove time offset---------- """
    prec = 3  # Saving precision
    powerflow = ["Block \t [kV] \t [MW] \t [MVAr] \t [rad] \t Area"]
    for block in ScanBlocksTool:
        if block.type == "AC":
            block.theta = block.snapshot_data["theta"][initial_row]  # Save the voltage angle
            Vd = block.snapshot_data["VDUTac:1"][initial_row]  # Peak L-N values
            Vq = block.snapshot_data["VDUTac:2"][initial_row]
            Id = block.snapshot_data["IDUTacA2:1"][initial_row]  # Peak
            Iq = block.snapshot_data["IDUTacA2:2"][initial_row]
            P = (Vd*Id + Vq*Iq)*3/2
            Q = (Vd*Iq - Vq*Id)*3/2
            V = np.sqrt(Vd**2 + Vq**2) * np.sqrt(3/2)  # L-L magnitude
            # print(" ",block.name,"-\t V =",round(V,3),"kV,\t angle =",round(block.theta,4),"rad,\t P =",round(P,3),"MW,\t Q =",round(Q,3),"MVAr")
            powerflow.append(block.name+"\t"+format(V, f".{prec}f")+"\t"+format(P, f".{prec}f")+"\t"+format(Q, f".{prec}f")+"\t"+format(block.theta, f".{prec}f")+"\t"+str(block.area))
        else:
            V = block.snapshot_data["VDUTdc"][initial_row]
            P = V*block.snapshot_data["IDUTdcA2"][initial_row]
            powerflow.append(block.name + "\t" + format(V, f".{prec}f") + "\t" + format(P, f".{prec}f"))

        for name in block.snapshot_data.keys():
            # Remove snapshot time offset
            if name == "time": 
                block.snapshot_data[name] = block.snapshot_data[name][initial_row:] - block.snapshot_data[name][initial_row]
            else:
                block.snapshot_data[name] = block.snapshot_data[name][initial_row:]

    if show_powerflow: print(" ------ Power flow after the snapshot ------")
    with open(results_folder+r'\\'+output_files+"_power_flow.txt", 'w') as powerflowfile:
        for item in powerflow:
            powerflowfile.write(item + '\n')
            if show_powerflow: print(" "," \t ".join(item.split()))

    print('\n Unperturbed simulation results collected in', round((t.time() - t1), 2), 'seconds')

    # Store the terminal angle and the area for each AC block
    with open(results_folder+r'\\'+output_files+"_angles.txt", 'w') as file:
        file.write('\t'.join(map(str, ["Bus", "Area", "Angle [rad]"])) + '\n')
        for block in ScanBlocksTool:
            if block.type == "AC":
                file.write('\t'.join(map(str,[block.name, block.area, block.theta]))+'\n')

    if scantype == "AC" and scan_actives:
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

        simset_task.parameters(volley=num_parallel_sim, affinity=1, ammunition=f_points)
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

        if not edge_dq_sym:
            # q-axis injection
            print(' Running single frequency q-axis injection simulations')
            t1 = t.time()
            for block in ScanBlocksTool:
                if block.type == "AC":
                    block.pscad_block.parameters(V_perturb_mag=v_perturb_mag,selector=2)  # q-axis injection
                else:
                    block.pscad_block.parameters(V_perturb_mag=v_perturb_mag,selector=0)  # No injection
            
            simset_task.parameters(volley=num_parallel_sim, affinity=1, ammunition=f_points)
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

    elif scantype == "DC" and scan_actives:
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
        simset_task.parameters(volley=num_parallel_sim, affinity=1, ammunition=f_points)
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

    elif scan_actives:
        # ACDC-type scan
        # d-axis injection
        print('\n Running single frequency d-axis injection simulations')
        t1 = t.time()
        for block in ScanBlocksTool:
            if block.type == "AC":
                block.pscad_block.parameters(V_perturb_mag=v_perturb_mag,selector=1)  # d-axis injection
            else:
                block.pscad_block.parameters(V_perturb_mag=v_perturb_mag,selector=0)  # No injection
        simset_task.parameters(volley=num_parallel_sim, affinity=1, ammunition=f_points)
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
        for block in ScanBlocksTool:
            if block.type == "AC":
                block.pscad_block.parameters(selector=2)  # q-axis injection
            else:
                block.pscad_block.parameters(selector=0)  # No injection
        simset_task.parameters(volley=num_parallel_sim, affinity=1, ammunition=f_points)
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
        simset_task.parameters(volley=num_parallel_sim, affinity=1, ammunition=f_points)
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

    if compute_yz and scan_actives:
        t2 = t.time()
        print(' Computing admittances')
        # Snapshot data time-aligment
        t_0 = ScanBlocksTool[0].perturbation_data["time"][0]
        if verbose:
            print('  Simulation and snapshot initial time:', t_0, ScanBlocksTool[0].snapshot_data["time"][0])
            print("  Shape snapshot:",ScanBlocksTool[0].snapshot_data["time"].shape,"perturb:",ScanBlocksTool[0].perturbation_data["time"].shape)

        # Shift the snapshot data by one time-step if needed
        if t_0 != 0.0 and round(ScanBlocksTool[0].snapshot_data["time"][0],10) != round(t_0,10) and \
                ScanBlocksTool[0].snapshot_data["time"].shape != ScanBlocksTool[0].perturbation_data["time"].shape:
            print("  Shifting snapshot")
            # ScanBlocksTool[0].snapshot_data["time"] = ScanBlocksTool[0].snapshot_data["time"][1:]
            for block in ScanBlocksTool:
                for name in block.snapshot_data.keys():
                    block.snapshot_data[name] = block.snapshot_data[name][1:]  # Time has been modified already

        # Start FFT index and sampling time used in the FFT
        dt = round(sample_step * 1e-6, 12)  # Sampling time [s]
        if verbose: print(" Sampling time: ", sample_step, " [us]")
        start_idx = find_nearest(ScanBlocksTool[0].perturbation_data["time"], start_fft)
        if verbose: print(" FFT time: ", start_fft, "Time vector value: ", ScanBlocksTool[0].perturbation_data["time"][start_idx])

        Ytopology_scan = np.copy(Ytopology)  # Make a copy to modify as the different scans are being done
        for idx, name in enumerate(block_names_Y):
            if sum(Ytopology_scan[idx,:]) == 1:
                # Only for AC/DC converters and shunt AC or DC components
                nz = np.nonzero(Ytopology_scan[idx,:])[0]  # Find the index for the other scan block
                # The loop looks for the point to point blocks, but it can be improved (see passives filt)
                for block in ScanBlocksTool:
                    if block.name == name[:-2]:
                        block_type0 = block.type  # The current block
                        block0 = block
                    if block.name == block_names_Y[int(nz)][:-2]:
                        block_type1 = block.type  # The other scan block
                        block1 = block
                if block_type0 != block_type1:
                    # If one block is AC and the other DC: AC/DC converter
                    t1 = t.time()
                    if block0.type == "DC":
                        sides = [name[-1], block_names_Y[int(nz)][-1]]
                        zblocks_pair = [block0,block1]  # DC block and AC block
                    else:
                        sides = [block_names_Y[int(nz)][-1], name[-1]]
                        zblocks_pair = [block1,block0]  # DC block and AC block

                    yz_computation.admittance(f_base=f_base, frequencies=freq, fft_periods=fft_periods,dt=dt,
                                              start_idx=start_idx,
                                              zblocks=zblocks_pair, sides=sides, scantype="ACDC",
                                              results_folder=results_folder, results_name=output_files)
                    # Update the scan matrix to indicate that no scan from and to these two ports is pending
                    Ytopology_scan[nz,idx] = 0
                    Ytopology_scan[idx,nz] = 0
                    print('  Admittance between',zblocks_pair[0].name+"-"+sides[0],'and',zblocks_pair[1].name+"-"+sides[1],'computed in',round((t.time() - t1), 2),'seconds')
                else:
                    # Both are AC or DC scan blocks: shunt components (one-side scan)
                    if block0.name == block1.name:
                        t1 = t.time()
                        yz_computation.admittance(f_base=f_base, frequencies=freq, fft_periods=fft_periods, dt=dt,
                                                  start_idx=start_idx, exploit_dq_sym=edge_dq_sym,
                                                  zblocks=block0, sides=name[-1], scantype=block0.type,
                                                  results_folder=results_folder, results_name=output_files)
                        # Update the scan matrix to indicate that no scan from and to these two ports is pending
                        Ytopology_scan[nz, idx] = 0
                        Ytopology_scan[idx, nz] = 0
                        print('  Admittance at',name,'computed in',round((t.time() - t1), 2),'seconds')

        print(' Admittance computation finished in', round((t.time() - t2), 2), 'seconds')

    """ Perform the simulations for the scan of the passive networks based on the topology information """
    # Shorter simulations can be set as EMT transients are usually faster, but for simplicity former settings are used
    # The PSCAD project folder (i.e. project_name.gf46) can be cleared here to decrease memory usage in large projects
    sim_select = {"d": 1, "q": 2, "dc": 1}  # Dict containing the perturbation type based on the name ending

    if verbose: print("Passive network scans:",[network_aux.names for network_aux in passive_networks_scans])
    if (len(passive_networks_scans) != 0) and compute_yz and scan_passives:
        # Passive network scan
        print("\nScan of AC and/or DC networks")
        t1 = t.time()

        # Create a scan schedule based on blocks connected to the same vs different subsystems
        create_scan_schedule(passive_networks_scans)
        max_runs = max([max(net.runs_list) for net in passive_networks_scans])
        if verbose:
            if verbose: print("Scheduled perturbations per run")
            for run in range(1, max_runs + 1):
                scans_this_run = []
                for network_scan in passive_networks_scans:
                    if run in network_scan.runs_list: scans_this_run.append(network_scan.scan_per_run[run])
                if verbose: print(" Run",str(run)+'/'+str(max_runs)+":",", ".join(scans_this_run))
            print("Starting multi-terminal scan with a total of",max_runs,"runs")

        for run in range(1, max_runs + 1):
            if verbose: print(' Run',str(run)+'/'+str(max_runs))
            # Disable all perturbations
            for block in ScanBlocksTool: block.pscad_block.parameters(V_perturb_mag=v_perturb_mag, selector=0)
            t2 = t.time()
            idx_selected_blocks = []  # This list changes every run so only the necessary blocks store the PSCAD results
            # Configure every PSCAD block involved in the network scan for this run
            for network_scan in passive_networks_scans:
                if run in network_scan.runs_list:                 
                    # If this network is scanned during this run: update the selected Z-blocks and iterate over the blocks & names
                    for block_idx in list(network_scan.blocks_idx.values()): idx_selected_blocks.append(block_idx)

                    # Enable the perturbation at the block of this run
                    block_side_perturbation = network_scan.scan_per_run[run]
                    block_name_w_side = "_".join(block_side_perturbation.split("_")[:-1])
                    perturbation = block_side_perturbation.split("_")[-1]  # d, q, or dc perturbation
                    ScanBlocksTool[network_scan.blocks_idx[block_name_w_side]].pscad_block.parameters(selector=sim_select[perturbation], single_frequency = 0 if multi_freq_scan else 1)
                    if verbose: print("\t",perturbation,"perturbation at",block_name_w_side)

            # Perform the simulations and label the output data by using the "run" number
            simset_task.parameters(volley=num_parallel_sim, affinity=1, ammunition=f_points_per_file if multi_freq_scan else f_points)
            simset_task.overrides(duration=t_sim_internal, time_step=t_step, plot_step=sample_step, start_method=1,
                                  timed_snapshots=0, startup_inputfile=snapshot_file + '.snp',
                                  save_channels_file=output_files + "_" + str(run) + '.out', save_channels=1)
            
            if run_sim: simset.run()

            print(' Run',str(run)+'/'+str(max_runs),'finished in', round((t.time() - t2), 2), 'seconds')
            if save_td or compute_yz:
                wait4pscad(time=1, pscad=pscad)
                t2 = t.time()
                read_and_save.multiple_s(n_sim= f_points_per_file if multi_freq_scan else f_points, out_folder=out_dir, save_folder=results_folder, save=save_td,
                                         tar_files=all_files_to_open,
                                         zblocks=[ScanBlocksTool[ind] for ind in idx_selected_blocks],
                                         file_name=simset_task.overrides()['save_channels_file'][:-4])
                if verbose: print(' Results collected in', round((t.time() - t2), 2), 'seconds\n')

        print(' Network scans completed in', round((t.time() - t1), 2), 'seconds\n')

        # Compute the admittance for each subnetwork individually
        for network_scan in passive_networks_scans:
            t2 = t.time()
            print(' Computing admittances')
            # Snapshot data time-aligment
            idx0 = idx_selected_blocks[0]
            t_0 = ScanBlocksTool[idx0].perturbation_data["time"][0]
            if verbose:
                print('  Simulation and snapshot initial time:', t_0, ScanBlocksTool[idx0].snapshot_data["time"][0])
                print("  Shape snapshot:",ScanBlocksTool[idx0].snapshot_data["time"].shape,"perturb:",ScanBlocksTool[idx0].perturbation_data["time"].shape)

            # Shift the snapshot data by one time-step if needed
            if t_0 != 0.0 and round(ScanBlocksTool[idx0].snapshot_data["time"][0],10) != round(t_0,10) and \
                    ScanBlocksTool[idx0].snapshot_data["time"].shape != ScanBlocksTool[idx0].perturbation_data["time"].shape:
                if verbose: print("  Shifting snapshot")
                for block in ScanBlocksTool:
                    for name in block.snapshot_data.keys():
                        block.snapshot_data[name] = block.snapshot_data[name][1:] 

            # Start FFT index and sampling time used in the FFT
            dt = round(sample_step * 1e-6, 12)  # Sampling time [s]
            if verbose: print(" Sampling time: ", sample_step, " [us]")
            start_idx = find_nearest(ScanBlocksTool[idx0].perturbation_data["time"], start_fft)
            if verbose: print(" FFT time: ", start_fft, "Time vector value: ", ScanBlocksTool[idx0].perturbation_data["time"][start_idx])

            idx_selected_blocks = []
            sides_selected_blocks = []
            for name in network_scan.names:
                idx_selected_blocks.append(network_scan.blocks_idx[name])
                sides_selected_blocks.append(name[-1])
            print(" Computing admittance for the network",', '.join(network_scan.names))

            if multi_freq_scan:
                yz_computation.admittance_multi_freq(f_base=f_base, freq_multi=freq_multi, fft_periods=fft_periods, dt=dt,
                                                     start_idx=start_idx, results_folder=results_folder,
                                                     zblocks=[ScanBlocksTool[ind] for ind in idx_selected_blocks],
                                                     sides=[side for side in sides_selected_blocks], network=network_scan,
                                                     results_name=output_files, exploit_dq_sym=edge_dq_sym)
            else:
                # yz_computation.admittance(f_base=f_base, frequencies=freq, fft_periods=fft_periods, dt=dt,
                #                         start_idx=start_idx, scantype="Network", results_folder=results_folder,
                #                         zblocks=[ScanBlocksTool[ind] for ind in idx_selected_blocks],
                #                         sides=[side for side in sides_selected_blocks], network=network_scan,
                #                         results_name=output_files, exploit_dq_sym=edge_dq_sym)
                yz_computation.admittance_generic(f_base=f_base, frequencies=freq, fft_periods=fft_periods, dt=dt,
                                                  start_idx=start_idx, results_folder=results_folder,
                                                  zblocks=[ScanBlocksTool[ind] for ind in idx_selected_blocks],
                                                  sides=[side for side in sides_selected_blocks], network=network_scan,
                                                  results_name=output_files, exploit_dq_sym=edge_dq_sym)
                
            print(' Admittance matrix involving',", ".join(network_scan.names),'computed in',round((t.time() - t2), 2), 'seconds')

    # Quit PSCAD
    wait4pscad(time=1, pscad=pscad)
    pscad.quit()

    print('\nTotal execution time', round((t.time() - t0) / 60, 2), 'minutes\n')

def frequency_sweep_TF(t_snap=None, t_sim=None, t_step=None, sample_step=None, v_perturb_mag=None, freq=None, f_points=None,
                       f_base=None, f_min=None, f_max=None, working_dir=None, multi_freq_scan=False,
                       snapshot_file=None, dedicated_SS_sim=False, take_snapshot=True, dt_injections=None,
                       project_name=None, workspace_name=None, fortran_ext=r'.gf46', num_parallel_sim=8,
                       component_parameters=None, results_folder=None, output_files='Perturbation', save_td=False,
                       fft_periods=1, start_fft=None, run_sim=True, verbose=False, pscad_plot=False):

    # Debugging control: run_sim enables or dissables all PSCAD simulations; verbose enables or disables script run info
    """ --- Input data handling --- """
    # CHECK the following if: it does not work... "or" operator gives a bool as output not None
    if (t_snap or t_step or start_fft or v_perturb_mag or project_name or workspace_name or ((f_points or f_base or f_max or f_min) and freq)) is None:
        print('One or more required arguments are missing!! \n Check the function documentation by typing help(frequency_sweep) \n')
        return

    # Create frequency list if it is not provided
    if freq is None: freq = create_freq.loglist(f_min=f_min, f_max=f_max, f_points=f_points, f_base=f_base)

    # If the sample time is not provided, it is set to half of the minimum required value (multiple of step_time)
    if sample_step is None:
        if f_max is None: f_max = max(freq)
        sample_step = round(t_step * np.floor((1e6 * 0.5 * 0.5 / f_max + t_step / 2) / t_step), 3)  # [us]
    if t_sim is None: t_sim = start_fft + fft_periods / f_base
    if dt_injections is None: dt_injections = start_fft  # dt_injections
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

    # Create the folder to store the results if it does not exist or use the current one
    if (results_folder is not None) and (not path.exists(results_folder)): makedirs(results_folder)
    if results_folder is None: results_folder = working_dir

    # Folder for the scan options
    if not path.exists(working_dir+"Scan_options"): makedirs(working_dir+"Scan_options")

    # Frequency vector to PSCAD XY table: single frequency perturbations
    freq_text_file='frequencies.txt'
    f_points = len(freq)
    with open(working_dir+"Scan_options\\"+freq_text_file, 'w') as f:  # Create the .txt file
        f.write('! This file stores the frequency sweep values in Hz \n')  # Write header
        for j in range(f_points): f.write(str(j + 1) + '\t' + str(freq[j]) + '\n')  # Write values
        f.write('ENDFILE:')  # Write end of the file
    f.close()

    # Frequency vector to several PSCAD XY tables: multi-frequency perturbations
    f_points_per_file = int(np.ceil(f_points / num_parallel_sim))  
    N_files = 8  # 8 frequencies at the same time by default
    freq_multi = np.pad(freq, (0, f_points_per_file*N_files - f_points), 'constant') # Add zeros if not enough frequencies
    freq_multi = freq_multi.reshape(N_files, f_points_per_file)
    for freq_file_num in range(N_files):
        with open(working_dir+"Scan_options\\"+freq_text_file[:-4]+str(freq_file_num+1)+".txt", 'w') as f:  # Create the .txt file
            f.write('! This file stores the multi-frequency sweep values in Hz \n')  # Write header
            for j in range(f_points_per_file):
                f.write(str(j + 1) + '\t' + str(freq_multi[freq_file_num][j]) + '\n')  # Write values
            f.write('ENDFILE:')  # Write end of the file is apparently not needed for this component
        f.close()

    """ --- Main program --- """
    t0 = t.time()  # Initial time
    print('Launching PSCAD')
    pscad = launch(minimize=True)
    wait4pscad(time=1, pscad=pscad)
    t.sleep(5)  # Wait a bit more just in case PSCAD is still loading
    if not pscad.licensed():
        certificates = pscad.get_available_certificates()  # Retreive PSCAD license certificate
        keys = list(certificates.keys())
        pscad.get_certificate(certificates[keys[0]])  # Get the first available certificate
    pscad.load(working_dir + workspace_name + ".pswx")  # Load workspace
    project = pscad.project(project_name)  # Where the model is defined
    main = project.canvas('Main')  # Main of the project
    project.focus()  # Switch PSCAD’s focus to the project

    print(' Setting parameters and configuration')
    # Look for the TF scan blocks in the main canvas
    blocks = main.find_all("Z_tool:TFscan")

    ScanBlocks_id = [i for i in range(1, len(blocks) + 1)]  # Unique scan block_id signals
    # Create a list with the active scan block objects containing rich information about each block
    ScanBlocksTool = []
    ScanBlocksTool_names = []  # Name of the blocks as they appear in ScanBlocksTool to avoid excesive iterations
    
    for idx, block in enumerate(blocks):
        # Set snapshot parameters and block ID in the scan blocks
        block.parameters(Tdecoupling=t_snap, T_inj=t_snap_internal, selector=0, block_id=ScanBlocks_id[idx])
        ScanBlocksTool.append(Scanblock(block, block.parameters()['Name'], int(block.parameters()['block_id'])))
        ScanBlocksTool[idx].perturbation_data = {i: {} for i in range(f_points)}  # Dict of dicts
        ScanBlocksTool_names.append(ScanBlocksTool[idx].name)
        print("TF scan block",ScanBlocksTool[idx].name,'with block_id',int(block.parameters()['block_id']))

    # Set variable component parameters
    if component_parameters is not None:
        for i in range(len(component_parameters)):
            parameter_component = main.find('master:const', component_parameters[i][0])
            parameter_component.parameters(Value=component_parameters[i][1])

    # Disable output channel components (all but Z-tool's scaning blocks) and animated displays
    if not pscad_plot:
        all_pgb = project.find_all("master:pgb")  # Find all output channels in the project
        for pgb in all_pgb:
            if not (pgb.parameters()['Name'] in ['blockid','inputTF', 'outputTF']):  pgb.disable()  # Disable the non-selected outputs
        all_multimeters = project.find_all("master:multimeter")  # Find all multimeters in the project
        for multimeter in all_multimeters: multimeter.parameters(Dis=0)  # Animated display is disabled (0)

    # Set simulation-specific parameters
    if 'Perturbation' not in pscad.simulation_sets():
        simset = pscad.create_simulation_set('Perturbation')
        simset.add_tasks(project_name)
        if verbose: print(' A simulation set has been created')
    else:
        simset = pscad.simulation_set('Perturbation')
    simset_task = simset.tasks()[0]  # The task is extracted

    if take_snapshot:  # It runs the snapshots
        print(' Running snapshot')
        t1 = t.time()
        if snapshot_file is None: snapshot_file = 'Snapshot'
        simset_task.parameters(volley=1, affinity=1, ammunition=1)
        if dedicated_SS_sim:
            # Run first just to take a snapshot without saving the data
            # This because with dll files the states can be more efficiently saved only at the end of the (snapshot) run
            simset_task.overrides(duration=t_snap_internal, time_step=t_step, plot_step=sample_step,
                                  start_method=0, timed_snapshots=1, snapshot_file=snapshot_file + '.snp',
                                  snap_time=t_snap_internal, save_channels_file=snapshot_file + '.out', save_channels=0)
            if run_sim: simset.run()
            # Run again to record the steady-state waveforms (no perturbations)
            simset_task.overrides(duration=t_sim_internal, time_step=t_step, plot_step=sample_step, start_method=1,
                                  timed_snapshots=0, startup_inputfile=snapshot_file + '.snp',
                                  save_channels_file=snapshot_file + '.out', save_channels=1)
            if run_sim: simset.run()
        else:
            # Run the snapshot and a bit longer to record the steady-state waveforms
            simset_task.overrides(duration=t_snap_internal + t_sim_internal, time_step=t_step, plot_step=sample_step,
                                  start_method=0, timed_snapshots=1, snapshot_file=snapshot_file + '.snp',
                                  snap_time=t_snap_internal, save_channels_file=snapshot_file + '.out', save_channels=1)
            if run_sim: simset.run()
        print(' Snapshot completed in', round((t.time() - t1), 2), 'seconds')
        
    elif snapshot_file is not None:
        # It performs the unperturbed simulation starting from the given snapshot (not fully tested yet)
        t1 = t.time()
        simset_task.parameters(volley=1, affinity=1, ammunition=1)  # affinity_type = "DISABLE_TRACING" disables the plotting
        simset_task.overrides(duration=t_sim_internal, time_step=t_step, plot_step=sample_step, start_method=1,
                              timed_snapshots=0, startup_inputfile=snapshot_file + '.snp',
                              save_channels_file=snapshot_file+'.out', save_channels=1)
        if dedicated_SS_sim:
            print(' Running steady-state simulation from a given snapshot')
            if run_sim: simset.run()
            print(' Steady-state simulation completed in', round((t.time() - t1), 2), 'seconds')
        # Otherwise, just read the outputs from previous scan

    # Identify the variables to be retrieved, their output channels and associated Z-tool scan blocks
    wait4pscad(time=1, pscad=pscad)
    t1 = t.time()
    ch_var_names = dict()
    out_filename = out_dir + "\\" + snapshot_file  # Snapshot output filename
    with open(out_filename + ".inf", 'r') as info_file:
        out_num = []  # Ztool's variables output channel number
        names = []  # Ztool's variables output names
        counter = 1  # Total number of PSCAD output signals
        for line in info_file.readlines():
            if line.split()[3].split('"')[1] in "TFscan":  # If the output channel corresponds to a Ztool variable
                out_num.append(counter)  # Get variable's output channel number
                names.append(line.split()[2].split('"')[1])  # Get output variable name
                ch_var_names[counter] = names[-1]  # Var name entry with the channel num as key
            counter = counter + 1
    if verbose: print("Output channel variable names \n","\n".join(names))
    block_id_out_num = [out_num[i] for i in range(len(names)) if "blockid" in names[i]]  # block_id outputs numbers
    files_to_open = [int(np.ceil(block_id_out / 10)) for block_id_out in block_id_out_num]  # With block_id outputs
    files_to_open = list(set(files_to_open))  # File's number to be opened
    block_id_out_signal = []
    for file_num in files_to_open:
        # Select the columns to be read relative to each file and only for the block_id signal
        cols = [num + 1 - 10 * (file_num - 1) for num in block_id_out_num if
                int(np.ceil(num / 10)) == file_num]
        if file_num < 10:  # If the file number is less than 10, then it adds 0 before the file number
            values = read_one_line(out_filename + "_0" + str(file_num) + ".out", nline=1)  # Read the first value
        else:
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
                ch1 = ch0 + len(['blockid','input_TF', 'output_TF'])  # End channel number containing the TFscan block signals
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
                    print(" Block ", block.name, "\n \t Target output channels ", block.out_vars_ch)
                    print(" \t Output names: ", [block.out_vars_names[ch] for ch in block.out_vars_ch])

                files_to_open = [int(np.ceil(block_out / 10)) for block_out in block.out_vars_ch]
                files_to_open = list(set(files_to_open))  # File's number to be opened (no repetitions)
                block.files_to_open = files_to_open  # Number of the files containing the block's outputs
                for f2o in files_to_open: all_files_to_open.append(f2o)
                for file_num in files_to_open:
                    # Select the columns to be read relative to each file and only for the signals of the scan block
                    cols = [num + 1 - 10 * (file_num - 1) for num in block.out_vars_ch if int(np.ceil(num / 10)) == file_num]
                    block.relative_cols[file_num] = cols
                    if verbose: print("\t Block ", block.name, " output file:", file_num,", columns:",cols)

    all_files_to_open = list(set(all_files_to_open))  # Get rid of repetitions

    # Save snapshot and steady-state run results
    read_and_save.single_s(out_files=out_filename, save_folder=results_folder,
                           save=save_td, files=all_files_to_open, zblocks=ScanBlocksTool,
                           new_file_name=simset_task.overrides()['save_channels_file'][:-4],scan_vars=['inputTF', 'outputTF'])

    # Remove the snapshot time offset and show power flow
    if dedicated_SS_sim:
        # snapshot_data contains already the steady-state simulation data starting from a dedicated snapshot file
        initial_row = 0
    else:
        # Snapshot data contains both the cold-start data and the steady-state waveforms
        initial_row = find_nearest(ScanBlocksTool[0].snapshot_data["time"], t_snap_internal)
        if verbose: print(" Snapshot w/o time offset", ScanBlocksTool[0].snapshot_data["time"][initial_row], t_snap_internal)

    for name in block.snapshot_data.keys():
        # Remove snapshot time offset
        if name == "time": 
            block.snapshot_data[name] = block.snapshot_data[name][initial_row:] - block.snapshot_data[name][initial_row]
        else:
            block.snapshot_data[name] = block.snapshot_data[name][initial_row:]
            
    print('\n Unperturbed simulation results collected in', round((t.time() - t1), 2), 'seconds')

    if multi_freq_scan:
        print(' Running multi-frequency perturbation simulations')
    else:
        print(' Running single-frequency perturbation simulations')
    t1 = t.time()
    idx_selected_blocks = []
    for idx, block in enumerate(ScanBlocksTool):
        idx_selected_blocks.append(idx)
        block.pscad_block.parameters(V_perturb_mag=v_perturb_mag, selector=1, single_frequency = 0 if multi_freq_scan else 1)

    simset_task.parameters(volley=num_parallel_sim, affinity=1, ammunition=f_points)
    simset_task.overrides(duration=t_sim_internal, time_step=t_step, plot_step=sample_step, start_method=1,
                            timed_snapshots=0, startup_inputfile=snapshot_file+'.snp',
                            save_channels_file=output_files+'_TF.out', save_channels=1)
    if run_sim: simset.run()
    print(' Simulation finished in', round((t.time() - t1), 2), 'seconds')

    wait4pscad(time=1, pscad=pscad)
    t2 = t.time()
    read_and_save.multiple_s(n_sim=f_points, out_folder=out_dir, save_folder=results_folder, save=save_td,
                             tar_files=all_files_to_open, zblocks=[ScanBlocksTool[ind] for ind in idx_selected_blocks],
                             file_name=simset_task.overrides()['save_channels_file'][:-4], scan_vars=['inputTF', 'outputTF'])
    print(' Simulation results collected in', round((t.time() - t2), 2), 'seconds\n')

    t2 = t.time()
    print(' Computing transfer functions')
    # Snapshot data time-aligment
    t_0 = ScanBlocksTool[0].perturbation_data["time"][0]
    if verbose:
        print('  Simulation and snapshot initial time:', t_0, ScanBlocksTool[0].snapshot_data["time"][0])
        print("  Shape snapshot:",ScanBlocksTool[0].snapshot_data["time"].shape,"perturb:",ScanBlocksTool[0].perturbation_data["time"].shape)

    # Shift the snapshot data by one time-step if needed
    if t_0 != 0.0 and round(ScanBlocksTool[0].snapshot_data["time"][0],10) != round(t_0,10) and \
            ScanBlocksTool[0].snapshot_data["time"].shape != ScanBlocksTool[0].perturbation_data["time"].shape:
        print("  Shifting snapshot")
        # ScanBlocksTool[0].snapshot_data["time"] = ScanBlocksTool[0].snapshot_data["time"][1:]
        for block in ScanBlocksTool:
            for name in block.snapshot_data.keys():
                block.snapshot_data[name] = block.snapshot_data[name][1:]  # Time has been modified already

    # Start FFT index and sampling time used in the FFT
    dt = round(sample_step * 1e-6, 12)  # Sampling time [s]
    if verbose: print(" Sampling time: ", sample_step, " [us]")
    start_idx = find_nearest(ScanBlocksTool[0].perturbation_data["time"], start_fft)
    if verbose: print(" FFT time: ", start_fft, "Time vector value: ", ScanBlocksTool[0].perturbation_data["time"][start_idx])

    for idx, name in enumerate(ScanBlocksTool_names):
        t1 = t.time()
        yz_computation.SISO_TF(f_base=f_base, frequencies=freq, fft_periods=fft_periods, dt=dt,
                               start_idx=start_idx, zblocks=[ScanBlocksTool[idx]],
                               results_folder=results_folder, results_name=output_files)
        print('  TF at', name,'computed in',round((t.time() - t1), 2),'seconds \n')

    print(' SISO TF computation finished in', round((t.time() - t2), 2), 'seconds')

    # Quit PSCAD
    wait4pscad(time=1, pscad=pscad)
    pscad.quit()

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


frequency_sweep.__doc__ = """
This function performs sinusoidal perturbations of electrical quantities directly inserted in the network effectively decoupling the system into smaller subsystems.

The default simulation configuration is to run a cold-start, decouple the subsystems, take a snapshot, and then run the frequency sweeps in different cores.
The function accepts several input arguments to customize the frequency scan:

Required
        workspace_name  Name of the PSCAD workspace: only the name, no path and without the extension (which is .pswx)
        project_name    Name of the PSCAD project under study (only name, no path and no extension which is .pscx)
        t_snap          Time when the snapshot is taken [s].
        t_sim			Duration of each frequency injection simulation [s].
        t_step			Simulation time step [us].
        start_fft		Time for the system to reach the sinusoidal steady-state after the injections [s]. 
        v_perturb_mag	Voltage perturbation in per unit w.r.t. the steady-state fundamental value, e.g. 0.01
        f_points		Number of frequency perturbation points. Ideally a multiple of the possible PSCAD multi-core simulations number.
        f_base		 	Base frequency (determines frequency resolution) [Hz]. The values of the frequency list are multiples of this base frequency.
        f_min			Start frequency [Hz].
        f_max			End frequency [Hz].
        topology        .txt absolute path and file name of the network topology between scan blocks (see examples for instructions on how to build it).
        output_files    Name root of the scan results output files.
      
Optional
        working_dir	 	        Absolute path of the PSCAD workspace in case the python file is not in the same folder as the PSCAD project.
        output_files	 	    Root name of the output files
        results_folder	 	    Absolute path where the tool results want to be stored. If not specified, they are saved in the working directory.
        fortran_ext	 	        Fortran extension. Default r'.gf46'.
        num_parallel_sim 	    Number of parallel simulations. Default = 8.
        
        sample_step		        Sample time of the output channels [us]. The default value is computed based on the Nyquist frequency.
        snapshot_file		    Name of the snapshot file so it can be re-used or in case a previous snapshot is used.
        take_snapshot	        Bool: Does the user want to take a snapshot? Default = True. A previous snapshot can still be used if snapshot_file is specified.
                                The snapshot simulation runs for t_snap_internal + t_sim_internal so as to also save the steady-state unperturbed waveforms.
        dedicated_SS_sim        Bool flag to perform a dedicated simulation just to record the steady-state waveforms

        scan_actives            Bool flag to scan the active components, i.e. AC and/or DC components identified by 1 in the diagonal entries of the topology file and/or AC/DC converters. Default = True.
        scan_passives           Bool flag to scan the passive networks. Default = True. It can be set to False in case the edge matrix does not need to be scanned.
        edge_dq_sym             Bool to consider the symmetry of passive AC networks so as to halve the necessary perturbations and thus their computation time.
        
        freq			        Frequency list to perform the injections [Hz]. Alternatively, the user can provide info to compute the list.
        fft_periods 		    Number of periods used to compute the FFT. Default = 1.
        multi_freq_scan         Bool: If set to True, it performs multi-sine injections to reduce the total computation time. Default = False.
        dt_injections           Additional simulation time in seconds after the system decoupling (steady-state). It can be set to zero.
        
        show_powerflow          Bool: do you want to print the power flow after the snapshot (steady-state)? Default = False.
        visualize_network       Bool: create a primitive graph representing the provided network topology. Default = False.
        
	    component_parameters	List of two-value lists containing the component parameter name and value to be modified in PSCAD. E.g. [["BRK_time", 2.50], ["P_load", 50]]
	    
        pscad_plot              Binary to disable PSCAD plots so as to speed-up the simulations. Default = 0 i.e. no plots
        compute_yz	 	        Compute the admittance and save the results. If no results folder is specified then it saves the data in working_dir. Default = True
        save_td  		        Bool: If set to True, many files of time domain data are saved into more compact .txt files for each independent
                                perturbation. The format is [time Vd(f1) Vq(f1) Id(f1) Iq(f1) ... Vd(f_max) Vq(f_max) Id(f_max) Iq(f_max)].
        run_sim                 Bool flag to run or not run all PSCAD simulations. True = runs PSCAD, False = does NOT run PSCAD
        verbose                 Bool flag to show detailed debugging and processing information.
"""

frequency_sweep_TF.__doc__ = """
This function performs sinusoidal perturbations of a single signal (TF scan block) to identify transfer functions.

The default simulation configuration is to run a cold-start, sample & hold the original input (opening the loop), take a snapshot, and then run the frequency sweeps in different cores.
The function accepts several input arguments to customize the frequency scan:

Required
        workspace_name  Name of the PSCAD workspace (only the name, no path and without the extension)
        project_name    Name of the PSCAD project under study 
        t_snap          Time when the snapshot is taken [s].
        t_sim			Duration of each frequency injection simulation [s].
        t_step			Simulation time step [us].
        start_fft		Time for the system to reach the sinusoidal steady-state after the injections [s]. 
                        Currently, the same time is considered regardless of the perturbation frequency but this could be improved in the future.
        v_perturb_mag	Amplitude perturbation in per unit w.r.t. the steady-state average value, e.g. 0.01
        f_points		Number of frequency perturbation points. Ideally a multiple of the possible PSCAD multi-core simulations number.
        f_base		 	Base frequency (determines frequency resolution) [Hz]. The values of the frequency list are multiples of this base frequency.
        f_min			Start frequency [Hz].
        f_max			End frequency [Hz].
      
Optional
        working_dir	 	        Absolute path of the PSCAD workspace in case the python file is not in the same folder as the PSCAD project.
        output_files	 	    Root name of the results output files
        results_folder	 	    Absolute path where the tool's results are to be stored. If not specified, they are saved in the working directory.
        fortran_ext	 	        Fortran extension. Default r'.gf46'.
        num_parallel_sim 	    Number of parallel simulations. Default = 8.
        
        sample_step		        Sample time of the output channels [us]. The default value is computed based on the Nyquist frequency.
        snapshot_file		    Name of the snapshot file so it can be re-used or in case a previous snapshot is used.
        take_snapshot	        Bool: Does the user want to take a snapshot? Default = True. A previous snapshot can still be used if snapshot_file is specified.
                                The snapshot simulation runs for t_snap_internal + t_sim_internal so as to also save the steady-state unperturbed waveforms.
        dedicated_SS_sim        Bool flag to perform a dedicated simulation just to record the steady-state waveforms

        freq			        Frequency list to perform the injections [Hz]. Alternatively, the user can provide info to compute the list.
        fft_periods 		    Number of periods used to compute the FFT. Default = 1.
        multi_freq_scan         Bool: If set to True, it performs multi-sine injections to reduce the total computation time. Default = False.
        dt_injections           Additional simulation time in seconds after the system decoupling (steady-state). It can be set to zero.
                
	    component_parameters	List of two-value lists containing the component parameter name and value to be modified in PSCAD. E.g. [["BRK_time", 2.50], ["P_load", 50]]
	    
        pscad_plot              Binary to disable PSCAD plots so as to speed-up the simulations. Default = 0 i.e. no plots
        save_td  		        Bool: If set to True, many files of time domain data are saved into more compact .txt files for each independent
                                perturbation. The format is [time Vd(f1) Vq(f1) Id(f1) Iq(f1) ... Vd(f_max) Vq(f_max) Id(f_max) Iq(f_max)].
        run_sim                 Bool flag to run or not run all PSCAD simulations. True = runs PSCAD, False = does NOT run PSCAD
        verbose                 Bool flag to show detailed debugging and processing information.
"""