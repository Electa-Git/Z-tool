"""
This program contains several functions for frequency-domain stability analysis including:
    0) Loading and building of (multi-infeed) subsystem matrices for subsequent small-signal analysis
    1) Generalized Nyquist Criterion (GNC) application to determine system stability: via eigenvalue decomposition and via the determinant
    2) Eigenvalue Decomposition (EVD) of the closed-loop matrix to determine oscillatory modes and bus participation factors
    3) Passivity index (for the application of the passivity theorem) and singular value decomposition (small-gain theorem) of target matrices
    4) A main stability_analysis function to apply all the previously described to a specific system

Copyright (C) 2026  Francisco Javier Cifuentes Garcia

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

__all__ = ['stability_analysis','passivity','nyquist','small_gain','EVD','nyquist_det','loci_sensitivity','mode_estimation']

import numpy as np
from scipy.optimize import linear_sum_assignment, least_squares
from scipy.signal import argrelmax
import matplotlib.pyplot as plt
from .read_admittance import read_admittance
from os import path, makedirs
from warnings import warn
import pickle
from .plot_utils import bode_plot, sparsity_plot
from matplotlib import rcParams  # Text's parameters for plots
rcParams['mathtext.fontset'] = 'cm'  # Font selection
rcParams['font.family'] = 'STIXGeneral'  # 'cmu serif'

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

def stability_analysis(topology=None, results_folder=None, file_root=None, indentations=[], node_blocks=None, rotate_edge=False, rotate_node=False, reference_buses=None, relative_angles=True,
                       check_conditioning=False, condition_number_th=10e6, make_plot=True, save_pickle=False, save_results=True, save_Y=False, save_loop_gain=False,
                       verbose=True, run_nyquist=True, run_nyquist_det=False, run_EVD=True, run_EVD_PFs=True, run_EVD_PFs_extended=False, run_passivity=True, run_small_gain=True,
                       run_GNC_sensitivity=False, normalize_GNC_sensitivity=False, run_PMD=False, run_sigma=False, modal_estimation_nyquist=False, modal_estimation_EVD=False,
                       order_maxima=4, extra_poles=0, samples_fitting=12, Ibase={}, Vbase={}, PMD_zeta_threshold=0.25):
    # This function loads and builds the edge and node admittance matrices and applies the most common stability analysis functions
    # 0) Firstly, read the terminal angle information for the AC blocks if rotations are required
    if rotate_edge or rotate_node:
        block_area_angle = []  # List for each block containing a list as [bus/block name, area_id, terminal angle in rad]
        with open(results_folder+r'\\'+file_root+'_angles.txt', 'r') as file:
            next(file)  # First line contains the header
            for line in file:
                content = [str(line.strip().split("\t")[0]),int(line.strip().split("\t")[1]),float(line.strip().split("\t")[2])]
                block_area_angle.append(content) # For each element in block_area_angle, the first element is the block name, the second is the area number, and the third is the angle
        
        # Define a reference bus for each area when using relative angles. Otherwise, use the absolute angles as in the '_angles.txt' file.
        if relative_angles:
            reference_angle = {i: -100 for i in set([bus[1] for bus in block_area_angle])}  # Unrealistic intial value to later define the reference bus for each are when relative angles are used
        else:
            reference_angle = {i: 0.0 for i in set([bus[1] for bus in block_area_angle])} # The reference angle is zero for each area so the absolute angles as in the '_angles.txt' file are used
        if verbose and relative_angles: print("The reference buses are:")
        if reference_buses is not None and relative_angles:
            for bus in block_area_angle:
                if bus[0] in reference_buses: 
                    reference_angle[bus[1]] = float(bus[2])  # If the block of this bus is in the reference_buses list, then use the angle of this bus as the reference for its area
                    if verbose: print(" "+bus[0]+", for area",bus[1])
        elif -100 in reference_angle.values(): # If the user does not specify all reference buses, asign them automatically as the first bus of each area in the file
            for bus in block_area_angle:
                if reference_angle[bus[1]] == -100: 
                    reference_angle[bus[1]] = float(bus[2])  # If the reference angle of this block's area is not defined yet, then use the angle of the current block as the reference
                    if verbose: print(" "+bus[0]+", for area",bus[1])
    
    if rotate_edge and rotate_node and verbose:
        print("Rotating both the edge and node matrices results in matrices with unchanged eigenvalues; this also applies to the open and closed-loop matrices.")

    # 1) Read the topology matrix and extract block names
    if topology is not None:
        Ytopology = np.loadtxt(topology, skiprows=1, comments=["#", "%", "!"])
        # nameA-1 nameA-2 nameB-1 nameB-2 ... x nameA-1 nameA-2 nameB-1 ...
        # 0 means no interconnection, 1 means connection between the edges: diagonals are single-sided / shunt
        with open(topology, 'r') as f:
            block_names_Y = f.readline().strip('\n').split("\t")
    elif node_blocks is not None:
        if len(node_blocks) != 1:
            raise ValueError("The topology must be specified when there is more than one node block.")
        # If the user specifies the node blocks but not the topology, then it is assumed that the interconnection is full
        block_names_Y = [block[:-2] for block in node_blocks for _ in (0, 1)]
        block_names_Y = [block+"-1" if i%2==0 else block+"-2" for i, block in enumerate(block_names_Y)]
        Ytopology =  np.identity(2)
    else:
        raise ValueError("The topology or node_blocks must be specified. Check the function documentation by typing help(stability_analysis)")

    # 2) Read the admittance files based on the topology file
    admittances = []  # List containing the admittance objects
    # Create the undirected graph: adjacent matrix but diagonals can be 1
    g = Graph(len(Ytopology))
    for row, name in enumerate(block_names_Y):
        for col, edge in enumerate(Ytopology[row]):
            if int(edge) == 1: g.addEdge(row, col)
    # Obtain the connected components of the graph
    cc = g.connectedComponents()  # List of lists with blocks positions connected
    # For every group of connected buses, read and add the admittance matrix
    for buses in cc:
        involved_blocks=[block_names_Y[bus] for bus in buses]
        if node_blocks is not None:
            node = True if set(involved_blocks).issubset(node_blocks) else False  # If all the blocks in the file are in the node_blocks list, then it is a node matrix
        else:
            node = None  # If the user does not specify the node blocks, AC/DC matrices and single-port AC or DC components are node matrices by default
        admittances.append(read_admittance(path=results_folder, involved_blocks=involved_blocks, file_root=file_root, node=node))

    # 3) Update bus_names to be ordered as the variables in the individual admittance matrices and build the node matrix
    node_matrix = []  # Create the node matrix with the active components (block diagonal)
    node_variables = []  # List of variable names = the current/voltage vectors
    edge_aux_matrix = []  # Auxiliary edge matrix (block diagonal)
    edge_aux_variables = []  # The order of this aux matrix is different from the order of the node matrix
    for y in admittances:
        # print(y.blocks,y.y_type,"- Node:",y.node)
        if y.node:
            # Define the nodal matrix that sets the order of the electrical variables
            for var in y.vars: node_variables.append(var)
            node_matrix.append(y)
        else:
            # Update the aux edge matrix and its variables
            for var in y.vars: edge_aux_variables.append(var)
            edge_aux_matrix.append(y)

    edge_ordering = []  # List to re-sort the edge matrix acording to the node matrix variables
    for var in node_variables: edge_ordering.append(edge_aux_variables.index(var))
    # print("Node vars \n",node_variables,"\nEdge vars \n",edge_aux_variables)
    # print("\nSorted edge variables \n",sorted(edge_aux_variables,key=node_variables.index))

    # 4) Create the node and edge matrices with the frequency domain data
    frequencies = admittances[0].f  # Retreive the frequency vector
    # Create the auxiliary edge matrix (different order than node matrix), useful to check the network topology and scan
    Yedge_aux = np.zeros((len(frequencies),len(node_variables),len(node_variables)),dtype='cdouble')  # Or dtype='csingle'
    y_edge_idx = 0
    for yedge in edge_aux_matrix:
        # # Eliminate too small elements: related to PSCAD accuracy when the topology is not enforced
        # for col in range(np.size(yedge.y, 1)):
        #     for row in range(np.size(yedge.y, 2)):
        #         if max(abs(yedge.y[:, row, col])) < 1e-6:  # The threshold is system and time-step dependent
        #             yedge.y[:, row, col] = np.zeros(yedge.y[:, row, col].shape)
        Yedge_aux[:, y_edge_idx:y_edge_idx + len(yedge.vars), y_edge_idx:y_edge_idx + len(yedge.vars)] = yedge.y
        y_edge_idx = y_edge_idx + len(yedge.vars)  # Update the matrix index for the next admittance block
    
    # Sort the edge matrix acording to the node matrix variables
    Yedge = Yedge_aux[:,:,edge_ordering]  # Sort the columns
    Yedge = Yedge[:,edge_ordering,:]  # Sort the rows

    # 5) Build the block-diagonal node admittance and define the rotation matrix (if needed)
    Ynode = np.zeros((len(frequencies), len(node_variables), len(node_variables)),dtype='cdouble')  # Or dtype='csingle'
    if rotate_edge or rotate_node: T = np.zeros((len(node_variables), len(node_variables)), dtype='double') # Block-diagonal rotation matrix
    y_node_idx = 0 # Block index for the node matrix
    for ynode in node_matrix:
        Ynode[:, y_node_idx:y_node_idx+len(ynode.vars), y_node_idx:y_node_idx+len(ynode.vars)] = ynode.y  # Block-by-block construction of the node matrix        
        # Define the rotation matrix for the current node block if needed
        if rotate_edge or rotate_node: 
            sub_block_idx = 0 # Block index for the node matrix sub-blocks
            for block in ynode.blocks:
                if ynode.blocks_info[block]["type"] == "AC": # If the current sub-block is AC then define the AC rotation
                    for each_block in block_area_angle: # Iterate over all blocks to find the angle of the current sub-matrix
                        if each_block[0] == block[:-2]: # If the names match
                            theta = (each_block[2] - reference_angle[each_block[1]]) # Rotation angle = block's bus angle - reference angle for the area of this block 
                            # print('Nodal submatrix',idx,"between:",y_node_idx,"and",y_node_idx+len(ynode.vars),"with angle",round(theta,4),"for block",block,"between",str(y_node_idx+sub_block_idx),"and",str(y_node_idx+sub_block_idx+2))
                            T[y_node_idx+sub_block_idx:y_node_idx+sub_block_idx+2, y_node_idx+sub_block_idx:y_node_idx+sub_block_idx+2] = np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])
                    sub_block_idx = sub_block_idx + 2
                
                elif ynode.blocks_info[block]["type"] == "DC": # DC blocks need no rotation
                    # print("Node submatrix",idx,"between:",y_node_idx,"and",y_node_idx+len(ynode.vars),"DC type for block",block,"between",str(y_node_idx+sub_block_idx),"and",str(y_node_idx+sub_block_idx+1))
                    T[y_node_idx+sub_block_idx:y_node_idx+sub_block_idx+1, y_node_idx+sub_block_idx:y_node_idx+sub_block_idx+1] = 1.0
                    sub_block_idx = sub_block_idx + 1
                    
        y_node_idx = y_node_idx + len(ynode.vars)  # Update the matrix index for the next admittance block

    # Rotate the system matrices: T is orthogonal so T^-1 = T'; Note that rotating both matrices does not alter the eigenvalues of their product or sum!
    if rotate_node: Ynode = T.transpose() @ Ynode @ T # Equivalent to np.matmul(T.transpose(),np.matmul(Ynode,T))
    if rotate_edge: Yedge = T.transpose() @ Yedge @ T 

    # Apply a conversion of the matrices to the per-unit system if the base values are provided by the user, i.e. Y_pu = Ib^-1 @ Y_SI @ Vb
    if Ibase and Vbase:
        Ib = [] # Base currents in the order of the node variables
        Vb = [] # Base voltages in the order of the node variables
        variables = [var[:-len(var.split("_")[-1])-1] for var in node_variables] # Remove the ending, i.e. "dc", "d" and "q" to retrieve the node names
        if all(variable in Ibase.keys() for variable in variables) and all(variable in Vbase.keys() for variable in variables):
            for variable in variables:
                Ib.append(Ibase[variable])
                Vb.append(Vbase[variable])
            Vb = np.diag(Vb) # Diagonal matrix with the voltage base at each node
            Ib_inv = np.diag([1/Ibase_i for Ibase_i in Ib]) # Diagonal matrix with the inverse of the current base quantities at each node
            if verbose: print(" The matrices are converted to per unit.")
            Yedge = Ib_inv @ Yedge @ Vb # Base conversion of the edge matrix
            Ynode = Ib_inv @ Ynode @ Vb # Base conversion of the node matrix
        else:
            print(" The keys in the base quantities provided by the user do not match the node variables. The matrices are not converted to per unit.")
    
    # Sparsity plot for verification at the lowest frequency
    if save_results:
        sparsity_plot(Yedge_aux[0,:,:], title='Auxiliary edge admittance matrix at '+format(frequencies[0], '.2f')+' Hz', results_folder=results_folder, file_name=file_root+"_Edge_aux",  variables=edge_aux_variables)
        sparsity_plot(Yedge[0,:,:], title='Edge admittance matrix at '+format(frequencies[0], '.2f')+' Hz', results_folder=results_folder, file_name=file_root+"_Edge",  variables=node_variables)
        sparsity_plot(Ynode[0,:,:], title='Node admittance matrix at '+format(frequencies[0], '.2f')+' Hz', results_folder=results_folder, file_name=file_root+"_Node",  variables=node_variables)

    # 6) Perform stability analysis
    if run_nyquist or run_nyquist_det or run_small_gain:
        Zedge = np.linalg.inv(Yedge)
        L = np.matmul(Zedge,Ynode)  # Loop gain matrix

    # Stability via eigenvalue loci
    if run_nyquist:
        GNC_results = nyquist(L, frequencies, results_folder, file_root, verbose=verbose, check_conditioning=check_conditioning, condition_number_th=condition_number_th, make_plot=make_plot,
                              indentations=indentations, run_sensitivity=run_GNC_sensitivity, Z=Zedge, Y=Ynode if normalize_GNC_sensitivity else None, bus_names=node_variables, run_sigma=run_sigma,
                              modal_estimation=modal_estimation_nyquist, extra_poles=extra_poles, order_maxima=order_maxima, samples_fitting=samples_fitting, save_pickle=save_pickle, save_results=save_results)
    
    # Stability via determinant
    if run_nyquist_det:
        GNC_det_results = nyquist_det(L,frequencies, results_folder, file_root, verbose=verbose, offset=0.0, draw_arrows=True, show_plot=False,
                                      make_plot=make_plot, indentations=indentations, save_pickle=save_pickle, save_results=save_results, run_sigma=run_sigma) 

    # Oscillatory frequencies and bus participation factors based on the closed-loop impedance: the admittance is provided to avoid inversion so Z_closedloop=False
    if run_EVD or run_PMD:
        EVD_results = EVD(Yedge+Ynode, frequencies, node_variables, results_folder, file_root, Z_closedloop=False, PFs=run_EVD_PFs, PFs_extended=run_EVD_PFs_extended,
                          run_PMD=run_PMD, modal_estimation=modal_estimation_EVD, extra_poles=extra_poles, order_maxima=order_maxima, samples_fitting=samples_fitting,
                          verbose=verbose, make_plot=make_plot, save_pickle=save_pickle, save_results=save_results, run_sigma=run_sigma, PMD_zeta_threshold=PMD_zeta_threshold)

    # Save the admittance matrices
    if save_results and save_Y:
        results = [Yedge[:, row, col] for row in range(len(node_variables)) for col in range(len(node_variables))]
        results.insert(0, frequencies)
        results = tuple(results)
        np.savetxt(results_folder + '\\' + file_root + '_Y_edge.txt', np.stack(results, axis=1), delimiter='\t',
                header="f\t" + "\t".join(node_variables), comments='')

        results = [Ynode[:, row, col] for row in range(len(node_variables)) for col in range(len(node_variables))]
        results.insert(0, frequencies)
        results = tuple(results)
        np.savetxt(results_folder + '\\' + file_root + '_Y_node.txt', np.stack(results, axis=1), delimiter='\t',
                header="f\t" + "\t".join(node_variables), comments='')

    # Save the minor loop gain matrix
    if (run_nyquist or run_nyquist_det) and save_loop_gain:
        results = [L[:, row, col] for row in range(len(node_variables)) for col in range(len(node_variables))]
        results.insert(0, frequencies)
        # elements = [str(row) + "-" + str(col) for row in range(len(node_variables)) for col in range(len(node_variables))]
        results = tuple(results)
        np.savetxt(results_folder+'\\'+file_root+'_Minor_loop_gain.txt',np.stack(results, axis=1),delimiter='\t',
                header="f\t" + "\t".join(node_variables), comments='')

    # Compute the passivity index and gains of the system matrices
    if run_passivity:
        passivity(G=Ynode,frequencies=frequencies,results_folder=results_folder,filename=file_root+"_Ynode", variables=node_variables, make_plot=make_plot, save_pickle=save_pickle, save_results=save_results)
        passivity(G=Yedge, frequencies=frequencies, results_folder=results_folder, filename=file_root + "_Yedge", make_plot=make_plot, save_pickle=save_pickle, save_results=save_results)
        passivity_index = passivity(G=Ynode+Yedge, frequencies=frequencies, results_folder=results_folder, filename=file_root+"_Ynode_+_Yedge", make_plot=make_plot, save_pickle=save_pickle, save_results=save_results)
    if run_small_gain:
        small_gain_index = small_gain(G2=Ynode, G1=Zedge, frequencies=frequencies, results_folder=results_folder, filename=file_root,
                                      variables=node_variables, make_plot=make_plot, save_pickle=save_pickle, save_results=save_results)

    return dict(nyquist=GNC_results if run_nyquist else None, stability_nyquist_det=GNC_det_results if run_nyquist_det else None, EVD=EVD_results if run_EVD else None,
                passivity_index=passivity_index if run_passivity else None, small_gain_index=small_gain_index if run_small_gain else None)

def passivity(G, frequencies, results_folder=None, filename='passivity', variables=None, Yedge=None, make_plot=True, save_pickle=False, save_results=True):
    # The passivity index is computed as half of the minimum eigenvalue of the matrix plus its conjugate transpose
    # min{eig(A + A')}/2
    # A passive system has its Nyquist plot in the RHP: increased stability margins if connected to a passive system
    passivity_index = np.real(np.min(np.linalg.eig(G + G.swapaxes(-1, -2).conj())[0], axis=1))/2  # min{eig(A + A')}/2
    # The eigenvalues of a Hermitian matrix are always real but floating-point arithmetic renders a small complex part
    # and thus we take the real value at the end to get rid of the spurious numerical artifact.

    if (results_folder is not None) and (filename is not None):
        if not path.exists(results_folder): makedirs(results_folder)  # Create results folder if it does not exist
        
        # Plot the passivity index over the frequency range
        if variables is None and make_plot:
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
            ax.plot(frequencies, passivity_index, color='blue', linestyle='solid', linewidth=2.0,label=r"$\mathbf{Y}_{node}$")
            ax.set_xscale("log")
            ax.minorticks_on()
            ax.grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
            # ax.grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
            ax.set_ylabel(r'min $\{ \lambda (\mathbf{G} + \mathbf{G}^H) \}/2$')
            ax.set_title('Passivity evaluation for ' + str(len(frequencies)) + ' frequencies')
            ax.set_xlim([frequencies[0], frequencies[-1]])
            ax.set_xlabel('Frequency [Hz]')
            if Yedge is not None:
                passivity_index_Yedge = np.real(np.min(np.linalg.eig(Yedge+Yedge.swapaxes(-1, -2).conj())[0], axis=1))/2
                ax.plot(frequencies, passivity_index_Yedge, color='red', linestyle='solid', linewidth=2.0,label=r"$\mathbf{Y}_{edge}$")
                ax.plot(frequencies, passivity_index + passivity_index_Yedge, color='green', linestyle='dashed', linewidth=2.0,label=r'min $\{ \lambda (\mathbf{Y}_{node} + \mathbf{Y}_{node}^H) \}/2$+min $\{ \lambda (\mathbf{Y}_{edge} + \mathbf{Y}_{edge}^H) \}/2$')
                ax.plot(frequencies, np.real(np.min(np.linalg.eig(G + G.swapaxes(-1, -2).conj() + Yedge + Yedge.swapaxes(-1, -2).conj())[0],axis=1))/2, color='black', linestyle='dotted',linewidth=2.0, label=r"$\mathbf{Y}_{node}+\mathbf{Y}_{edge}$")

                ax.legend(loc='upper left', fancybox=True, shadow=True, ncol=2)
        elif make_plot:
            # Find the position of the block diagonal matrices as they are surrounded by zeros
            indices = []  # Tuple of start and end indeces of each matrix
            start_index = 0
            for index in range(G.shape[1]-1):
                if G[0, index, start_index] == 0:  # If the next element is zero, then boundary of block matrix if defined
                    indices.append((start_index, index-1))  # There is a block matrix between start_index and index
                    start_index = index  # Initialize start index of next matrix
            indices.append((start_index, G.shape[1] - 1))  # Last block matrix

            # Plot the passivity of the different matrices and the whole matrix
            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
            for block_pos in range(len(indices)):
                start_idx = indices[block_pos][0]
                end_idx = indices[block_pos][1] + 1
                A = G[:, start_idx:end_idx, start_idx:end_idx]
                passivity_index_block = np.real(np.min(np.linalg.eig(A + A.swapaxes(-1, -2).conj())[0], axis=1))/2
                ax[0].plot(frequencies, passivity_index_block, linestyle='solid', linewidth=2.0,label=", ".join(variables[start_idx:end_idx]))

            ax[0].set_xscale("log")
            ax[0].minorticks_on()
            ax[0].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
            # ax[0].grid(visible=False)
            # ax[0].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
            ax[0].set_ylabel(r'min $\{ \lambda (\mathbf{G}_i + \mathbf{G}_i^H) \}/2$')
            ax[0].set_title('Passivity analysis for ' + str(len(frequencies)) + ' frequencies')
            ax[0].set_xlim([frequencies[0], frequencies[-1]])
            ax[0].legend(loc='upper right', ncol=int(np.floor(np.sqrt(len(indices)))), prop={'size': 6}) # fancybox=True, shadow=True, 

            ax[1].plot(frequencies,passivity_index,color='blue',linestyle='solid',linewidth=2.0,label="Complete matrix")
            ax[1].set_xscale("log")
            ax[1].minorticks_on()
            ax[1].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
            # ax[1].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
            ax[0].set_ylabel(r'min $\{ \lambda (\mathbf{G} + \mathbf{G}^H) \}/2$')
            ax[1].legend(loc='upper left', fancybox=True, shadow=True, ncol=1)
            ax[1].set_xlim([frequencies[0], frequencies[-1]])
            ax[1].set_xlabel('Frequency [Hz]')
        
        if make_plot:
            fig.savefig(results_folder+'\\'+filename + "_passivity.pdf", format="pdf", bbox_inches="tight")
            if save_pickle:
                with open(results_folder + '\\' + filename + "_passivity.pickle", 'wb') as f: pickle.dump(fig, f)
            plt.close(fig)
        
        if save_results:
            np.savetxt(results_folder+'\\'+filename+'_passivity.txt', np.stack((frequencies, passivity_index), axis=1), delimiter='\t', header="f\t" + "Passivity_index", comments='')  

    return passivity_index

def nyquist(L, frequencies, results_folder=None, filename='nyquist', verbose=True, check_conditioning=False, condition_number_th=0.01/5e-9, make_plot=True, show_plot=False, indentations =[], save_pickle=False, save_results=True,
            run_sensitivity=False, Z=None, Y=None, bus_names=None, unstable_frequency=False, modal_estimation=False, extra_poles=0, order_maxima=4, samples_fitting=12, verbose_modal_estimation=False, run_sigma=False):
    # Generalized Nyquist Criteria (GNC) for stability analysis: graphically determine the number of unstable closed-loop poles from the encirclements of the critical point by the loci of the open-loop gain matrix
    if verbose: print("Performing Nyquist stability assessment based on the eigenvalues of L")
    if not path.exists(results_folder): makedirs(results_folder)  # Create results folder if it does not exist

    # 1) Compute the eigenvalues of the loop-gain at every frequency
    eigenvalues, right_eigenvectors = np.linalg.eig(L)
    # Compute the condition number to discard doubtful data: threshold based on input error
    if check_conditioning: condition_number = np.linalg.cond(L)  # Relative error output <= cond_num * relative error input

    # Sorting of eigenvalues for a continuous eigenloci in the Nyquist plot based on minimum changes between frequencies
    eigenvalues_sorted = np.empty(eigenvalues.shape, dtype='cdouble')  # Initialization of sorted eigenvalues
    eigenvalues_sorted[0,:] = eigenvalues[0,:]  # First frequency as reference for the sorting
    for idx in range(1,eigenvalues.shape[0]):
        if check_conditioning:
            if condition_number[idx] > condition_number_th:
                eigenvalues_sorted[idx, :] = eigenvalues_sorted[idx-1, :] # Replace by previous well-conditioned values
        else:
            eig_1 = eigenvalues_sorted[idx - 1, :]  # Previous eigenvalues
            eig_2 = eigenvalues[idx, :]  # Current eigenvalues
            # Create matrix of distances between eigenvalues
            x, y = np.meshgrid(np.real(eig_1), np.real(eig_2), indexing='ij')
            d_real = np.abs(x - y)
            x, y = np.meshgrid(np.imag(eig_1), np.imag(eig_2), indexing='ij')
            d_imag = np.abs(x - y)
            d_abs = np.sqrt(np.square(d_real) + np.square(d_imag))  # Absolute distance by element-wise operations
            # Solve the linear sum assignment problem to find the minimum variation and thus the correct order
            col_ind = linear_sum_assignment(d_abs)[1]  # The absolute distance is the cost matrix
            eigenvalues_sorted[idx, :] = eigenvalues[idx, col_ind]  # Sort the eigenvalues
            right_eigenvectors[idx, :] = right_eigenvectors[idx][:,col_ind]  # Sort the eigenvectors

    # 2) Eigenloci plot and count clockwise and counter-clockwise encirclements of (-1,0j) for each eigenvalue
    # Compute the coordinates of the eigenloci for the GNC aplication and plotting
    x = np.real(eigenvalues_sorted)  # Real axis
    y = np.imag(eigenvalues_sorted)  # Imaginary axis
    
    # Consider only indentations strictly in the frequency range
    valid_indent = (indentations > frequencies[0]) & (indentations < frequencies[-1])
    inds = np.asarray(indentations)[valid_indent]

    j = np.searchsorted(frequencies, inds, side='left') # Candidate neighbor frequencies
    is_match = frequencies[j] == inds
    # Indices to the left of or at the indentations, and to the right of the indentations:
    # - If match:       (match_idx, match_idx + 1)
    # - Else (no match): (j-1, j)
    left_idx  = np.where(is_match, j,     j - 1)
    right_idx = np.where(is_match, j + 1, j)
    idx_pairs = np.stack([left_idx, right_idx], axis=1)
    idx_indentations = np.ravel(idx_pairs).tolist() # Save the two frequency indexes closest to each indentation frequency

    to_drop = np.unique(np.ravel(idx_pairs)) if idx_pairs.size else np.array([], dtype=int)
    keep_idx = np.ones(frequencies.size, dtype=bool)
    if to_drop.size: keep_idx[to_drop] = False # Boolean array to keep the points excluding those around the indentations
    
    x_indent = x.copy() # Make the values around the indentation NaN for plotting purposes
    y_indent = y.copy()
    x_indent[~keep_idx] = np.nan 
    y_indent[~keep_idx] = np.nan
    
    if len(indentations)>0 and verbose:
        for idx in range(0,len(idx_indentations),2):
            print(f"GNC indentation at {indentations[idx//2]} Hz performed between {frequencies[idx_indentations[idx]]} and {frequencies[idx_indentations[idx+1]]} Hz")
    
    cw = []  # List of clockwise crossings for each locus
    ccw = []  # List of counter-clockwise crossings for each locus
    unstable_loci = []  # List of unstable loci indeces (those that encircle the critical point)
    if make_plot:
        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(6, 7))  # Create the figure and get the colors cycle
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for col in plt.rcParams['axes.prop_cycle'].by_key()['color']: colors.append(col)  # Triplicate the color cycle
        for col in plt.rcParams['axes.prop_cycle'].by_key()['color']: colors.append(col)  # in case of many eigenvalues
    # Loop over the sorted eigenvalues for ploting the locus and count the crossings
    for idx in range(eigenvalues_sorted.shape[1]):
        # Plot the eigenvalue locus avoiding the indentations
        if make_plot:
            ax[0].plot(x_indent[:,idx], y_indent[:,idx], color=colors[idx], linestyle='solid', linewidth=1.5, label=r'$\lambda_{'+format(idx+1,'.0f')+r'}$')
            ax[0].plot(x_indent[:,idx], -y_indent[:,idx], color=colors[idx], linestyle='solid', linewidth=1.5, label='_nolegend_')
            ax[1].plot(x_indent[:,idx],y_indent[:,idx], color=colors[idx],linestyle='solid',linewidth=1.5,label=r'$\lambda_{' + format(idx+1,'.0f')+r'}$')
            ax[1].plot(x_indent[:,idx], -y_indent[:,idx], color=colors[idx], linestyle='solid', linewidth=1.5, label='_nolegend_')
        # Count the number of real axis crossings to the left of (-1, 0j) by this eigenvalue locus
        cwi = 0  # Initialize the counters
        ccwi = 0
        for j in range(1,eigenvalues_sorted.shape[0]):
            # Only consider clockwise crossings of the real axis beyond (-1,0j)
            if y[j - 1, idx] < 0 < y[j, idx] and (x[j,idx] < -1) and j not in idx_indentations:  # x[j-1,idx] < -1 or
                # Check that the (-1,0j) is to the right of the line between (x1,y1) and (x2,y2)
                # If the cross product of vectors (x2-x1, y2-y1) and (-1-x1, 0-y1) is < 0, then (-1,0) is to the right
                if (x[j,idx] - x[j-1,idx])*(0 - y[j-1,idx]) - (y[j,idx] - y[j-1,idx])*(-1 - x[j-1,idx]) < 0:
                    cwi += 1
                    if verbose: print("Real axis CW crossing at",round(0.5*(frequencies[j] + frequencies[j-1]),4),"Hz by lambda =",str(idx+1))
                    if make_plot and show_plot:
                        fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(6, 7))
                        ax1.plot([x[j-1,idx],x[j,idx]],[y[j-1,idx],y[j,idx]], color='red', linestyle='solid', linewidth=2.0, label='_nolegend_')
                        ax1.scatter(x[j - 1, idx], y[j - 1, idx], color='green', label=str(frequencies[j-1]))
                        ax1.scatter(x[j, idx], y[j, idx], color='blue', label=str(frequencies[j]))
                        ax1.scatter(-1, 0, marker="+", c='black', label=r'$( -1, 0j )$')
                        ax1.legend(loc='upper right', ncol=1)
                        ax1.minorticks_on()
                        ax1.grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
                        if show_plot: plt.show()  # Visualize the plot 
                        with open(results_folder + '\\' + filename + "_GNC_lambda_"+str(idx)+"_cw_"+str(cwi)+".pickle", 'wb') as f: pickle.dump(fig1, f)
                        plt.close(fig1)

            # Only consider counter-clockwise crossings of the real axis beyond (-1,0j)
            elif y[j - 1, idx] > 0 > y[j, idx] and (x[j-1,idx] < -1 or x[j,idx] < -1) and j not in idx_indentations:
                if (x[j,idx] - x[j-1,idx])*(0 - y[j-1,idx]) - (y[j,idx] - y[j-1,idx])*(-1 - x[j-1,idx]) > 0:
                    # If the cross product of vectors (x2-x1, y2-y1) and (-1-x1, 0-y1) is > 0, then (-1,0) is to the left
                    # if (x[j,idx] - x[j-1,idx])*(0 - y[j-1,idx]) - (y[j,idx] - y[j-1,idx])*(-1 - x[j-1,idx]) > 0:
                    ccwi += 1
                    if verbose: print("Real axis CCW crossing at",round(0.5*(frequencies[j] + frequencies[j-1]),4),"Hz by lambda =",str(idx+1))
                    if make_plot and show_plot:
                        fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(6, 7))
                        ax1.plot([x[j-1,idx],x[j,idx]],[y[j-1,idx],y[j,idx]], color='red', linestyle='solid', linewidth=2.0, label='_nolegend_')
                        ax1.scatter(x[j - 1, idx], y[j - 1, idx], color='green', label=str(frequencies[j-1]))
                        ax1.scatter(x[j, idx], y[j, idx], color='blue', label=str(frequencies[j]))
                        ax1.scatter(-1, 0, marker="+", c='black', label=r'$( -1, 0j )$')
                        ax1.legend(loc='upper right', ncol=1)
                        ax1.minorticks_on()
                        ax1.grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
                        if show_plot: plt.show()  # Visualize the plot
                        with open(results_folder + '\\' + filename + "_GNC_lambda_"+str(idx)+"_ccw_"+str(ccwi)+".pickle", 'wb') as f: pickle.dump(fig1, f)
                        plt.close(fig1)
        
        if cwi - ccwi > 0:
            if verbose: print("CW encirclements for lambda =",str(idx+1)+":",cwi-ccwi)
            unstable_loci.append(idx) # Add the locus index to the list of unstable loci

        cw.append(cwi)  # Add the counters to the list
        ccw.append(ccwi)

    # print("CC: ",cw,"\nCCW: ",ccw)
    N = sum(cw) - sum(ccw)  # Net number of clockwise encirclements
    if N > 0:
        stable_system = False
        if verbose: print("\n GNC stability assessment: UNSTABLE closed-loop system \n")
    elif N < 0:
        stable_system = False
        if verbose: print("\n GNC stability assessment: UNSTABLE subsystem \n")
    else:
        stable_system = True
        if verbose: print("\n GNC stability assessment: STABLE closed-loop system if subsystems are stable \n")

    # Plot the unit circle and the critical point
    if make_plot:
        th = np.linspace(-np.pi * 1.01, np.pi * 1.01, 314)
        ax[0].plot(np.cos(th), np.sin(th), color='black', linestyle='dotted', linewidth=1.0, label='Unit circle')
        ax[0].scatter(-1, 0, marker="+", c='blue', label=r'$( -1, 0j )$')
        ax[0].minorticks_on()
        ax[0].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
        # ax[0].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
        ax[0].set_title('Eigenloci between '+format(frequencies[0],'.1f')+' and '+format(frequencies[-1],'.1f') + ' Hz')
        ax[0].set_xlim([np.min(x,axis=None), np.max(x,axis=None)])
        ax[0].set_ylim([np.min(np.concatenate((-y,y)),axis=None), np.max(np.concatenate((-y,y)),axis=None)])
        ax[0].set_xlabel('Real axis')
        ax[0].set_ylabel('Imaginary axis')
        ax[0].legend(loc='upper right', ncol=4, prop={'size': 6})

        ax[1].plot(np.cos(th), np.sin(th), color='black', linestyle='dotted', linewidth=1.0, label='Unit circle')
        ax[1].scatter(-1, 0,s=4*rcParams['lines.markersize'] ** 2, marker="+", c='blue', label=r'$( -1, 0j )$')
        ax[1].set_xlim([-2.0, 2.0])
        ax[1].set_ylim([-2.0, 2.0])
        ax[1].set_xlabel('Real axis')
        ax[1].set_ylabel('Imaginary axis')
        ax[1].minorticks_on()
        ax[1].grid(visible=False) # which='major', color='k', linestyle='-', linewidth=0.5
        fig.savefig(results_folder + '\\' + filename + "_GNC.pdf", format="pdf", bbox_inches="tight")
        if save_pickle:
            with open(results_folder + '\\' + filename + "_GNC.pickle", 'wb') as f: pickle.dump(fig, f)
        if show_plot: plt.show()  # Visualize the plot
        plt.close(fig)
        bode_plot(Y=1/(1+eigenvalues_sorted),  frequencies=frequencies, results_folder=results_folder, file_name=filename+"_inv(1+L)", style="solid",
                  title=r"Bode plot of $1/(1+\lambda(L))$ over "+str(len(frequencies))+' frequencies', legend=[str(idx+1) for idx in range(eigenvalues_sorted.shape[1])])
    
    # Unstable closed-loop poles estimation by analyzing the maxima of 1/(1+L)
    if modal_estimation or unstable_frequency:
        unstable_modes = {}  # Dict of identified unstable oscillatory modes per locus (key is locus, modes are stored in a list)
        mode_samples = samples_fitting//2  # Number of samples around the unstable frequency to consider for mode estimation
        for locus in unstable_loci:
            unstable_modes[locus] = [] # Initialize the list of unstable modes for this locus
            unstable_freqs, unstable_dampings = unstable_frequency(eigenvalues_sorted[:,locus], frequencies, results_folder=results_folder, filename=filename+"_inv(1+L"+str(locus+1)+")", order_maxima=order_maxima, make_plot=make_plot)
            N_locus = abs(cw[locus] - ccw[locus])  # Net number of encirclements for this locus
            if len(unstable_freqs) == 0 and N_locus>0: # No unstable frequency identified by the method: consider the closest to the critical point
                idx_closest = argrelmax(-np.abs(eigenvalues_sorted[:,locus] + 1.0), order=order_maxima)[0]  # Index with the smallest distance to the critical point
                unstable_freqs = [frequencies[i] for i in idx_closest[0:N_locus]] # If no unstable frequency is identified, then use the closest points to the critical point
                unstable_dampings = [-0.0]*N_locus  # Consider zero damping
            # Consider only the lowest damped modes up to the number of encirclements by this locus
            idx_zeta_sort = list(np.argsort(unstable_dampings)) # Sort the dampings in ascending order: most unstable first
            unstable_freqs = [unstable_freqs[i] for i in idx_zeta_sort[0:N_locus]] # Sort up to the number of encirclements
            unstable_dampings = [unstable_dampings[i] for i in idx_zeta_sort[0:N_locus]] 
            for freq_idx, freq in enumerate(unstable_freqs):
                if modal_estimation:
                    idx_mode = np.argmin(np.abs(frequencies - freq)) # Find the index of the frequency closest to the unstable frequency
                    idx_min = idx_mode-mode_samples-1 if idx_mode-mode_samples-1 >= 0 else 0 # Minimum index for mode estimation larger than zero
                    idx_max = idx_mode+mode_samples-1 if idx_mode+mode_samples-1 <= len(frequencies) else len(frequencies) # Maximum index for mode estimation within array bounds
                    mode_parameters = mode_estimation(1/(1+eigenvalues_sorted[idx_min:idx_max,locus]), 2*np.pi*np.array(frequencies[idx_min:idx_max]), zeta0=unstable_dampings[freq_idx], omega0=2*np.pi*freq, extra_poles=extra_poles, verbose=verbose_modal_estimation)
                    sigma = mode_parameters[0]
                    omega = mode_parameters[1]
                    unstable_modes[locus].append(sigma+1j*omega) # Add the identified unstable mode to the list of unstable modes for this locus
                else:
                    # Assuming a low-pass second-order system, we can estimate the unstable mode from the frequency and damping ratio
                    wn = 2*np.pi*freq/np.sqrt(1-2*unstable_dampings[freq_idx]**2) if 2*np.abs(unstable_dampings[freq_idx])**2 < 1 else 2*np.pi*freq  # Natural frequency of the unstable mode in rad/s
                    sigma = -wn*unstable_dampings[freq_idx]  # Real part of the unstable mode
                    omega = wn*np.sqrt(1-unstable_dampings[freq_idx]**2) if 1>unstable_dampings[freq_idx] else wn  # Imaginary part of the unstable mode
                    unstable_modes[locus].append(sigma+1j*omega)
                if verbose: print(f"Unstable mode at {freq:.2f} Hz with z = {-sigma/np.sqrt(sigma**2 + omega**2):.4e} from pole: {sigma:.4e} +/- {np.abs(omega):.4e}j rad/s")

    if (make_plot or save_results) and run_sensitivity:
        idx_eigen_closest = np.argmin(np.min(np.abs(eigenvalues_sorted + 1.0), axis=0))  # Index of the eigen-locus with the smallest distance to the critical point
        selected_loci = [idx_eigen_closest] if stable_system else unstable_loci
        diag_sensitivity = loci_sensitivity(right_eigenvectors, np.linalg.inv(right_eigenvectors), frequencies, results_folder=results_folder, bus_names=bus_names,
                                            selected_loci=selected_loci, loci=eigenvalues_sorted, Z=Z, normalize=True if Y is not None else False, Y=Y,
                                            filename=filename+"_GNC_sens_of_locus", make_plot=make_plot, save_pickle=save_pickle, save_results=save_results)

    # Save the eigenloci
    if save_results:
        loci = [eigenvalues_sorted[:, idx] for idx in range(eigenvalues_sorted.shape[1])]
        loci.insert(0, frequencies)
        loci = tuple(loci)
        np.savetxt(results_folder + '\\' + filename + '_GNC.txt', np.stack(loci, axis=1), delimiter='\t',
                header="f\t"+"\t".join(["lambda_"+format(idx+1,'.0f') for idx in range(eigenvalues_sorted.shape[1])]), comments='')

    if run_sigma:
        sigmas = np.linalg.svd(np.identity(L.shape[1]) + L, compute_uv=False)
        bode_plot(np.min(sigmas, axis=1), frequencies, results_folder, filename+"_GNC_sigma", title='Minimum singular value of $I + L(j\omega)$ over '+str(len(frequencies))+' frequencies',
                  legend=["\sigma_{min}"], style="solid", save_pickle=save_pickle, save_data=save_results)

    return dict(stability = stable_system, unstable_loci = unstable_loci, unstable_modes = unstable_modes if modal_estimation or unstable_frequency else {}, diag_sensitivity = diag_sensitivity if run_sensitivity else None, sigmas = sigmas if run_sigma else None)

def small_gain(G2, frequencies,  G1=None, results_folder=None, filename='small_gain', variables=None, make_plot=True, save_pickle=False, save_results=True):
    # Applies a conservative version of the small-gain theorem as |L| = |G1*G2| <= |G1|*|G2| < 1
    S2 = np.linalg.svd(G2, compute_uv=False)
    S2_max = np.max(S2, axis=1)
    if G1 is None: G1 = np.eye(G2.shape[1])[None, :, :].repeat(G2.shape[0], axis=0)  # If G1 is not provided, consider it as an identity matrix
    S1 = np.linalg.svd(G1, compute_uv=False)
    S1_max = np.max(S1, axis=1)
    S1_max_times_S2_max = np.multiply(S1_max, S2_max)
    S12 = np.linalg.svd(np.matmul(G1,G2), compute_uv=False)
    S12_max = np.max(S12, axis=1)

    if not path.exists(results_folder): makedirs(results_folder)  # Create results folder if it does not exist

    if make_plot:
        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
        ax[0].plot(frequencies, 1.0/S1_max, color='blue', linestyle='solid', linewidth=2.0, label=r"1 / max $\sigma (\mathbf{G}_1)$")
        if variables is not None:
            # Find the position of the block diagonal matrices as they are surrounded by zeros
            indices = []  # Tuple of start and end indeces of each matrix
            start_index = 0
            for index in range(G2.shape[1] - 1):
                if G2[0, index, start_index] == 0:  # If the next element is zero, then boundary of block matrix if defined
                    indices.append((start_index, index - 1))  # There is a block matrix between start_index and index
                    start_index = index  # Initialize start index of next matrix
            indices.append((start_index, G2.shape[1] - 1))  # Last block matrix

            # Plot the maximum singular values of the different block matrices and the whole matrix
            for block_pos in range(len(indices)):
                start_idx = indices[block_pos][0]
                end_idx = indices[block_pos][1] + 1
                S2_block = np.max(np.linalg.svd(G2[:,start_idx:end_idx,start_idx:end_idx], compute_uv=False), axis=1)
                ax[0].plot(frequencies, S2_block, linestyle='solid',linewidth=2.0, label=", ".join(variables[start_idx:end_idx]))
        ax[0].plot(frequencies, S2_max, color='red', linestyle='dashed', linewidth=2.0,label=r"max $\sigma (\mathbf{G}_2)$")

        # Setings for upper plot
        ax[0].set_xscale("log")
        ax[0].set_yscale("log")
        ax[0].minorticks_on()
        ax[0].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
        # ax[0].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
        ax[0].set_ylabel(r'max $\sigma ( \cdot )$')
        ax[0].set_title('Singular value analysis over ' + str(len(frequencies)) + ' frequencies')
        ax[0].set_xlim([frequencies[0], frequencies[-1]])
        # ax[0].set_xlabel('Frequency [Hz]')
        ax[0].legend(loc='best', fancybox=True, shadow=True, ncol=1, prop={'size': 6})

        ax[1].plot(frequencies, S12_max, color='black', linestyle='solid', linewidth=2.0, label=r"max $\sigma (\mathbf{G}_1  \mathbf{G}_2)$")
        ax[1].plot(frequencies, S1_max_times_S2_max, color='green', linestyle='dashed',
                linewidth=2.0, label=r"max $\sigma (\mathbf{G}_1) \cdot $ max $\sigma (\mathbf{G}_2)$")
        ax[1].plot([frequencies[0], frequencies[-1]],[1, 1], color='grey', linestyle='dotted',linewidth=2.0, label='_nolegend_')
        ax[1].set_xscale("log")
        ax[1].set_yscale("log")
        ax[1].minorticks_on()
        ax[1].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
        # ax[1].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
        ax[1].set_ylabel('Unitless')
        ax[1].set_xlim([frequencies[0], frequencies[-1]])
        ax[1].set_xlabel('Frequency [Hz]')
        ax[1].legend(loc='best', fancybox=True, shadow=True, ncol=1, prop={'size': 6})
        
        fig.savefig(results_folder + '\\' + filename + "_gain.pdf", format="pdf", bbox_inches="tight")
        if save_pickle:
            with open(results_folder + '\\' + filename + "_gain.pickle", 'wb') as f: pickle.dump(fig, f)
        plt.close(fig)

    if save_results:
        np.savetxt(results_folder + '\\' + filename + '_gain.txt',
                    np.stack((frequencies, S1_max, S2_max, S12_max), axis=1),
                    delimiter='\t', header="f\t" + "max_sigma_G1\t" + "max_sigma_G2\t" + "max_sigma_G1_G2", comments='')
        
    return S1_max_times_S2_max

def EVD(G, frequencies, bus_names=None, results_folder=None, filename='EVD', verbose=True, Z_closedloop=True, make_plot=True, save_pickle=False, save_results=True,
        PFs=True, PFs_extended=False, run_PMD=False, modal_estimation=False, order_maxima=4, extra_poles=0, samples_fitting=12, run_sigma=False, PMD_zeta_threshold=0.707):
    if bus_names is None: bus_names = [str(bus+1) for bus in range(G.shape[1])]  # Sorted numbers if names not provided
    if not path.exists(results_folder): makedirs(results_folder)  # Create results folder if it does not exist

    # 1) Eigenvalue decomposition over the frequency
    eigenvalues, right_eigenvectors = np.linalg.eig(G)

    if not Z_closedloop: eigenvalues = 1.0 / eigenvalues  # Eigenvalues of G^-1 are the inverse of the eigenvalues of G
    
    # Sorting of eigenvalues for a continuous plot based on minimum changes between adjacent frequencies
    eigenvalues_sorted = np.empty(eigenvalues.shape, dtype='cdouble')  # Initialization of sorted eigenvalues
    eigenvalues_sorted[0, :] = eigenvalues[0, :]  # First frequency as reference for the sorting
    for idx in range(1, eigenvalues.shape[0]):
        eig_1 = eigenvalues_sorted[idx - 1, :]  # Previous eigenvalues
        eig_2 = eigenvalues[idx, :]  # Current eigenvalues
        # Create matrix of distances between eigenvalues
        x, y = np.meshgrid(np.real(eig_1), np.real(eig_2), indexing='ij')
        d_real = np.abs(x - y)
        x, y = np.meshgrid(np.imag(eig_1), np.imag(eig_2), indexing='ij')
        d_imag = np.abs(x - y)
        d_abs = np.sqrt(np.square(d_real) + np.square(d_imag))  # Absolute distance by element-wise operations
        # Solve the linear sum assignment problem to find the minimum variation and thus the correct order
        col_ind = linear_sum_assignment(d_abs)[1]  # The absolute distance is the cost matrix
        eigenvalues_sorted[idx, :] = eigenvalues[idx, col_ind]  # Sort the eigenvalues
        right_eigenvectors[idx, :] = right_eigenvectors[idx][:,col_ind]  # Sort the eigenvectors

    left_eigenvectors = np.linalg.inv(right_eigenvectors)

    # 2) Oscillation modes identification
    lambda_re = np.real(eigenvalues_sorted)  # Real part
    lambda_imag = np.imag(eigenvalues_sorted)  # Imaginary part
    lambda_abs = np.abs(eigenvalues_sorted)  # Absolute value (magnitude)

    # Oscillation mode identification based on the magnitude peaks of the closed-loop impedance matrix
    idx_lambda_max = np.argmax(lambda_abs,axis=0)  # Frequency index of the maximum magnitude of each eigenvalue
    idx_lambda_max_max = np.argmax([lambda_abs[idx_lambda_max[idx],idx] for idx in range(eigenvalues.shape[1])])  # Critical mode = the highest mag peak
    freq_idx = idx_lambda_max[idx_lambda_max_max]  # Main oscillation frequency index
    if verbose: print("The main oscillation frequency is around",round(frequencies[freq_idx],2),"Hz based on the magnitude of eigenvalue",idx_lambda_max_max+1,"=",np.round(eigenvalues_sorted[idx_lambda_max[idx_lambda_max_max],idx_lambda_max_max], 5))
    idx_lambda_envelope = np.argmax(lambda_abs, axis=1)  # Index of the maximum magnitude eigenvalue at each frequency

    if run_PMD or modal_estimation:
        # Apply the Positive Mode Damping (PMD) criterion, more information here: https://doi.org/10.1016/j.ijepes.2023.108957
        lambda_envelope = np.take_along_axis(eigenvalues_sorted, idx_lambda_envelope[:,None], axis=1)  # Eigenvalues with the largest magnitude at each frequency
        critical_points = argrelmax(np.abs(lambda_envelope), order=order_maxima)[0]
        if verbose: print("Critical frequencies:",frequencies[critical_points])
        PMD_indexes = [] # List of PMD criterion indexes: index > 0 -> potential unstable mode
        modes = [] # Initialize the list of unstable modes
        mode_samples = samples_fitting//2  # Number of samples around the unstable frequency to consider for modal estimation
        for critical_point in critical_points:
            if modal_estimation:
                # Perfom mode estimation around the critical frequencies by fitting a low-order rational function
                idx_min = critical_point-mode_samples-1 if critical_point-mode_samples-1 >= 0 else 0 # Minimum index for mode estimation larger than zero
                idx_max = critical_point+mode_samples-1 if critical_point+mode_samples-1 <= len(frequencies) else len(frequencies) # Maximum index for mode estimation 
                mode_parameters = mode_estimation(np.squeeze(lambda_envelope)[idx_min:idx_max], 2*np.pi*np.array(frequencies[idx_min:idx_max]), zeta0=0.05, omega0=2*np.pi*frequencies[critical_point], extra_poles=extra_poles)
                sigma = mode_parameters[0]
                omega = mode_parameters[1]
                modes.append(sigma+1j*omega) # Add the identified mode to the list of modes
                zeta = -sigma/np.sqrt(sigma**2 + omega**2) # Damping ratio of the identified mode
                if verbose: print(f"Mode at {frequencies[critical_point]:.2f} Hz with z = {zeta:.4e} from fitted pole: {sigma:.4e} +/- {np.abs(omega):.4e}j rad/s")
            else:
                zeta = 0.0 # Define a low damping ratio so the PMD criterion is applied below if run_PMD is True but no mode estimation is performed

            # Run the PMD only if the estimated damping ratio magnitude is below a threshold or if there is no modal estimation
            if run_PMD and np.abs(zeta) < PMD_zeta_threshold:
                # Run the PMD: three-point formula to approximate the derivate of the imaginary part around the local maxima of the mangnitude
                d1 = frequencies[critical_point] - frequencies[critical_point-1]
                d2 = frequencies[critical_point+1] - frequencies[critical_point]
                c_minus = -d2 / ( d1*(d1+d2) )
                c_zero  =  (d2 - d1) / ( d1*d2 )
                c_plus  =  d1 / ( d2*(d1+d2) )
                dImag_df = c_minus*np.imag(lambda_envelope[critical_point-1]) + c_zero*np.imag(lambda_envelope[critical_point]) + c_plus*np.imag(lambda_envelope[critical_point+1])
                PMD_index = dImag_df*np.real(lambda_envelope[critical_point]) # PMD criterion index: if > 0 then this is a potential unstable mode
                PMD_indexes.append(PMD_index)
                if PMD_index > 0 and verbose:
                    print("Unstable mode by the PMD criterion at",round(frequencies[critical_point], 2),"Hz for eigenvalue", np.round(lambda_envelope[critical_point], 5))

    # 3) Compute the bus participation factors (PFs) of the critical eigenvalue at the oscillation frequency
    # Controllability (right eigenvectors) and observability (transpose of left eigenvectors)
    Obs = right_eigenvectors[freq_idx, :]
    Cont = np.transpose(left_eigenvectors[freq_idx, :])
    # PF[frequency, row = bus, column = mode]
    if PFs:
        PF = right_eigenvectors * left_eigenvectors.transpose(0,2,1)  # Element-wise product of right eigenvectors and the transposed left eigenvectors
        PF_envelope = np.take_along_axis(PF, idx_lambda_envelope[:,None,None], axis=2)[:, :, 0]  # PF of the maximum magnitude eigenvalue at each frequency
    PF_freq_idx = Obs * Cont  # Element-wise product
    PF_mode = PF_freq_idx[:,idx_lambda_max_max]  # Select the target mode

    # The controllability, observability and PF of the critical mode at each bus
    if verbose:
        len_longest_bus_name = max([len(bus) for bus in bus_names])
        print("Bus"+(len_longest_bus_name-3)*" "+"\t","Cont.\t","Obs.\t","PF")  # Header
        for idx, bus in enumerate(bus_names):
            print(bus+(len_longest_bus_name-len(bus))*" "+"\t",
                  f"{np.abs(Cont[idx,idx_lambda_max_max]) / np.sum(np.abs(Cont[:,idx_lambda_max_max])):.4f}"+"\t",
                  f"{np.abs(Obs[idx,idx_lambda_max_max]) / np.sum(np.abs(Obs[:,idx_lambda_max_max])):.4f}"+"\t",
                  f"{np.abs(PF_mode[idx]) / np.sum(np.abs(PF_mode)):.4f}")

    # 4) Plot the eigenvalues over frequency
    if make_plot:
        fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(6, 8))  # Create the figure and get the colors cycle
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for col in plt.rcParams['axes.prop_cycle'].by_key()['color']: colors.append(col)  # Triplicate colour cycle
        for col in plt.rcParams['axes.prop_cycle'].by_key()['color']: colors.append(col)  # in case of many eigenvalues
        # Loop over the sorted eigenvalues and plot them over the frequency range
        for idx in range(eigenvalues_sorted.shape[1]):
            # Plot the eigenvalue locus
            ax[0].plot(frequencies, lambda_abs[:,idx], color=colors[idx], linestyle='solid', linewidth=1.5,
                    label=r'$\lambda_{' + format(idx + 1, '.0f') + r'}$')
            ax[1].plot(frequencies, lambda_re[:, idx], color=colors[idx], linestyle='solid', linewidth=1.5,
                    label=r'$\lambda_{' + format(idx + 1, '.0f') + r'}$')
            ax[2].plot(frequencies, lambda_imag[:, idx], color=colors[idx], linestyle='solid', linewidth=1.5,
                    label=r'$\lambda_{' + format(idx + 1, '.0f') + r'}$')
        # Figure settings and save to pdf
        ax[0].minorticks_on()
        ax[0].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
        # ax[0].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
        ax[0].set_title('Eigenvalue decomposition between ' + format(frequencies[0], '.1f') + ' and ' + format(frequencies[-1], '.1f') + ' Hz')
        ax[0].set_xlim([frequencies[0], frequencies[-1]])
        ax[0].set_ylim([np.min(lambda_abs, axis=None), np.max(lambda_abs, axis=None)])
        ax[0].set_ylabel('Magnitude')
        ax[0].set_xscale("log")
        ax[0].set_yscale("log")
        ax[0].legend(loc='lower right', ncol=int(np.ceil(np.sqrt(eigenvalues.shape[1]))), prop={'size': 6})

        ax[1].minorticks_on()
        ax[1].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
        # ax[1].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
        ax[1].set_xlim([frequencies[0], frequencies[-1]])
        ax[1].set_ylim([np.min(lambda_re, axis=None), np.max(lambda_re, axis=None)])
        ax[1].set_ylabel('Real part')
        ax[1].set_xscale("log")

        ax[2].minorticks_on()
        ax[2].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
        # ax[2].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
        ax[2].set_xlim([frequencies[0], frequencies[-1]])
        ax[2].set_ylim([np.min(lambda_imag, axis=None), np.max(lambda_imag, axis=None)])
        ax[2].set_xlabel('Frequency [Hz]')
        ax[2].set_ylabel('Imaginary part')
        ax[2].set_xscale("log")

        # plt.show()  # Visualize the plot interactively
        fig.savefig(results_folder + '\\' + filename + "_EVD.pdf", format="pdf", bbox_inches="tight")

        if save_pickle:
            with open(results_folder + '\\' + filename + "_EVD.pickle", 'wb') as f: pickle.dump(plt.gcf(), f)
        plt.close(fig)

        bode_plot(eigenvalues_sorted, frequencies, results_folder, filename+"_EVD_Bode", title='Modal impedance: EVD of the impedance matrix over '+str(len(frequencies))+' frequencies',
                  legend=[format(idx+1,'.0f') for idx in range(eigenvalues_sorted.shape[1])], style="solid", save_pickle=save_pickle)
        if PFs:
            bode_plot(PF_envelope, frequencies, results_folder, filename+"_EVD_max_PFs", title='Sensitivity of the largest modal impedance w.r.t. its diagonal elements',
                       legend=bus_names, style="solid", save_pickle=save_pickle, linear_mag=True)
        if PFs_extended:
            # Sensitivity of the critical modal impedance to each matrix element
            loci_sensitivity(right_eigenvectors, left_eigenvectors, frequencies, loci=eigenvalues_sorted, selected_loci=[idx_lambda_max_max], bus_names=bus_names, wrt_all_elements=PFs_extended,
                             results_folder=results_folder, filename=filename+"_EVD_PFs_extended_locus", make_plot=make_plot, save_pickle=save_pickle, save_results=save_results, normalize=False)   

     # Save the EVD and PFs of the envelope into a text file
    if save_results:
        np.savetxt(results_folder + '\\' + filename + '_EVD.txt', np.column_stack((frequencies, eigenvalues_sorted)), delimiter='\t', header="Frequency [Hz]\t" + "\t".join([str(i+1) for i in range(len(bus_names))]), comments='')
        if PFs:
            np.savetxt(results_folder + '\\' + filename + '_EVD_max_PFs.txt', np.column_stack((frequencies, PF_envelope)), delimiter='\t', header="Frequency [Hz]\t" + "\t".join(bus_names), comments='')
            np.savetxt(results_folder + '\\' + filename + '_EVD_max_controllability.txt', np.column_stack((frequencies, np.take_along_axis(right_eigenvectors, idx_lambda_envelope[:,None,None], axis=2)[:, :, 0])), delimiter='\t', header="Frequency [Hz]\t" + "\t".join(bus_names), comments='')
            np.savetxt(results_folder + '\\' + filename + '_EVD_max_observability.txt', np.column_stack((frequencies, np.take_along_axis(left_eigenvectors.transpose(0,2,1), idx_lambda_envelope[:,None,None], axis=2)[:, :, 0])), delimiter='\t', header="Frequency [Hz]\t" + "\t".join(bus_names), comments='')
    
    if run_sigma:
        sigmas = np.linalg.svd(G, compute_uv=False)
        bode_plot(np.max(sigmas, axis=1) if Z_closedloop else 1/np.min(sigmas, axis=1), frequencies, results_folder, filename+"_sigma_CL", title='Maximum singular value of the closed-loop dynamics over '+str(len(frequencies))+' frequencies',
                  legend=["\sigma_{max}"], style="solid", save_pickle=save_pickle, save_data=save_results)

    return dict(stability = not any(np.array(PMD_indexes) > 0) if run_PMD else None, modes = modes if modal_estimation else [], sigmas = sigmas if run_sigma else None, modal_impedances = eigenvalues_sorted, PFs = PF_envelope if PFs else PF_mode)

def nyquist_det(L, frequencies, results_folder=None, filename='nyquist_det', verbose=True, offset=0.0, draw_arrows=True, make_plot=True, show_plot=False, f0=50.0,
                indentations=[], run_sigma=False, save_pickle=False, save_results=True):
    # Stability assessment based on the determinant of I + L
    if verbose: print("Performing Nyquist stability assessment based on det(I + L) +",offset)
    if not path.exists(results_folder): makedirs(results_folder)  # Create results folder if it does not exist
    
    # Compute the determinant
    det = np.linalg.det(np.identity(L.shape[1]) + L) + offset  
    x = np.real(det)
    y = np.imag(det)

    # Consider only indentations strictly in the frequency range
    valid_indent = (indentations > frequencies[0]) & (indentations < frequencies[-1])
    inds = np.asarray(indentations)[valid_indent]
    
    j = np.searchsorted(frequencies, inds, side='left') # Candidate neighbor frequencies
    is_match = frequencies[j] == inds
    # Indices to the left of or at the indentations, and to the right of the indentations:
    # - If match:       (match_idx, match_idx + 1)
    # - Else (no match): (j-1, j)
    left_idx  = np.where(is_match, j,     j - 1)
    right_idx = np.where(is_match, j + 1, j)

    idx_pairs = np.stack([left_idx, right_idx], axis=1)
    idx_indentations = np.ravel(idx_pairs).tolist() # Save the two frequency indexes closest to each indentation frequency

    to_drop = np.unique(np.ravel(idx_pairs)) if idx_pairs.size else np.array([], dtype=int)
    keep_idx = np.ones(frequencies.size, dtype=bool)
    if to_drop.size: keep_idx[to_drop] = False # Boolean array to keep the points excluding those around the indentations
    
    x_indent = x.copy()
    y_indent = y.copy()
    x_indent[~keep_idx] = np.nan # Make the values around the indentation NaN to avoid plotting lines across them
    y_indent[~keep_idx] = np.nan
    
    if len(indentations)>0 and verbose:
        for idx in range(0,len(idx_indentations),2):
            print(f"GNC indentation at {indentations[idx//2]} Hz performed between {frequencies[idx_indentations[idx]]} and {frequencies[idx_indentations[idx+1]]} Hz")

    # Count the number of crossings of the positive vertical axis at (offset,0j) by det
    cwi = 0  # Initialize the counters
    ccwi = 0
    for j in range(1,len(frequencies)):
        # Clockwise crossings of the positive vertical axis at (offset,0j) avoiding the indentations
        if x[j - 1] < offset < x[j] and j not in idx_indentations:
            # Check that the (offset,0j) is to the right of the line between (x1,y1) and (x2,y2)
            # If the cross product of vectors (x2-x1, y2-y1) and (offset-x1, 0-y1) is < 0, then (offset,0) is to the right
            if (x[j] - x[j-1])*(0 - y[j-1]) - (y[j] - y[j-1])*(offset - x[j-1]) < 0:
                cwi += 1
                if show_plot:
                    print("Vertical line CW crossing at",round(0.5*(frequencies[j] + frequencies[j-1]),4),"Hz")
                    fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(6, 7))
                    ax1.plot([x[j-1],x[j]],[y[j-1],y[j]], color='red', linestyle='solid', linewidth=2.0, label='_nolegend_')
                    ax1.scatter(x[j-1], y[j-1], color='green', label=str(frequencies[j-1]))
                    ax1.scatter(x[j], y[j], color='blue', label=str(frequencies[j]))
                    ax1.scatter(offset, 0, marker="+", c='black', label=r'$( -1, 0j )$')
                    ax1.legend(loc='upper right', ncol=1)
                    ax1.minorticks_on()
                    ax1.grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
                    with open(results_folder + '\\' + filename + "_det_cw_"+str(cwi)+".pickle", 'wb') as f: pickle.dump(fig1, f)
                    if show_plot: plt.show() 
                    plt.close(fig1)

        # Counter-clockwise crossings of the positive vertical axis at (offset,0j) avoiding the indentations
        elif x[j - 1] > offset > x[j] and j not in idx_indentations:
            if (x[j] - x[j-1])*(0 - y[j-1]) - (y[j] - y[j-1])*(offset - x[j-1]) > 0:
                # If the cross product of vectors (x2-x1, y2-y1) and (offset-x1, 0-y1) is > 0, then (offset,0) is to the left
                ccwi += 1
                if show_plot:
                    print("Vertical line CCW crossing at",round(0.5*(frequencies[j] + frequencies[j-1]),4),"Hz")
                    fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(6, 7))
                    ax1.plot([x[j-1],x[j]],[y[j-1],y[j]], color='red', linestyle='solid', linewidth=2.0, label='_nolegend_')
                    ax1.scatter(x[j-1], y[j-1], color='green', label=str(frequencies[j-1]))
                    ax1.scatter(x[j], y[j], color='blue', label=str(frequencies[j]))
                    ax1.scatter(offset, 0, marker="+", c='black', label=r'$( -1, 0j )$')
                    ax1.legend(loc='upper right', ncol=1)
                    ax1.minorticks_on()
                    ax1.grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
                    with open(results_folder + '\\' + filename + "_det_ccw_"+str(ccwi)+".pickle", 'wb') as f: pickle.dump(fig1, f)
                    if show_plot: plt.show() 
                    plt.close(fig1)

    N = cwi - ccwi  # Net number of clockwise encirclements
    if N > 0:
        stable_system = False
        if verbose: print("\n GNC stability assessment: UNSTABLE closed-loop system \n")
    elif N < 0:
        stable_system = False
        if verbose: print("\n GNC stability assessment: UNSTABLE subsystem \n")
    else:
        stable_system = True
        if verbose: print("\n GNC stability assessment: STABLE closed-loop system if subsystems are stable \n")

    if make_plot:
        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(6, 7))  # Create the figure
        ax[0].scatter(offset, 0, marker="+", c='blue', label=r'$( '+str(round(offset,0))+', 0j )$') # Plot the critical point
        ax[0].plot(x_indent, y_indent, color='red', linestyle='solid', linewidth=1.5)
        if draw_arrows:
            da = int(np.log(x_indent.shape[0] + 1))  # decimate the number of arrows
            ax[0].quiver(x_indent[::da], y_indent[::da], np.gradient(x_indent)[::da], np.gradient(y_indent)[::da], angles='xy', scale_units='xy', scale=da, linewidth=1, edgecolor='black', facecolor='green')
        id0 = np.argmin(np.abs(frequencies - f0))
        ax[0].text(x[id0], y[id0], str(round(frequencies[id0], 2)) + ' Hz', fontsize=12, color='blue', ha='right', va='bottom')
        ax[0].text(x[-1], y[-1], str(round(frequencies[-1], 2)) + ' Hz', fontsize=12, color='blue', ha='right', va='bottom')
        ax[0].text(x[0], y[0], str(round(frequencies[0], 2)) + ' Hz', fontsize=12, color='blue', ha='right', va='bottom')

        ax[0].minorticks_on()
        ax[0].grid(visible=True, which='major', color='k', alpha=0.5, linestyle='-', linewidth=0.5)
        ax[0].set_title(str(offset)+r' + det[I + L(s)] between '+format(frequencies[0], '.1f')+' and ' + format(frequencies[-1], '.1f') + ' Hz')
        ax[0].set_xlim([np.min(x, axis=None), np.max(x, axis=None)])
        ax[0].set_ylim([np.min(y, axis=None), np.max(y, axis=None)])
        ax[0].set_xlabel('Real axis')
        ax[0].set_ylabel('Imaginary axis')
        ax[0].legend(loc='best', ncol=1)

        ax[1].plot(x_indent, y_indent, color='red', linestyle='solid', linewidth=2.0)
        # if draw_arrows:
        #     ax[1].quiver(x_indent, y_indent, np.gradient(x_indent), np.gradient(y_indent), angles='xy', scale_units='xy', scale=1, linewidth=1, edgecolor='black', facecolor='green')
        ax[1].scatter(offset, 0, s=4 * rcParams['lines.markersize'] ** 2, marker="+", c='blue', label=r'$( ' + str(offset) + ', 0j )$')
        id0 = np.argmin(np.abs(x + y*1j - offset)) # Closest point to the critical point
        ax[1].text(x[id0], y[id0], str(round(frequencies[id0], 2)) + ' Hz', fontsize=12, color='blue', ha='right', va='bottom')
        ax[1].set_xlim([-1.0-1.2*np.abs(x[id0]), 1.0+1.2*np.abs(x[id0])])
        ax[1].set_ylim([-1.0-1.2*np.abs(y[id0]), 1.0+1.2*np.abs(y[id0])])
        ax[1].minorticks_on()
        ax[1].grid(visible=True, which='major', color='k', alpha=0.5, linestyle='-', linewidth=0.5)
        ax[1].set_xlabel('Real axis')
        ax[1].set_ylabel('Imaginary axis')

        fig.savefig(results_folder + '\\' + filename + "_GNC_det.pdf", format="pdf", bbox_inches="tight")
        if save_pickle:
            with open(results_folder + '\\' + filename + "_GNC_det.pickle", 'wb') as f: pickle.dump(fig, f)
        if show_plot: plt.show()  # Visualize the plot interactively
        plt.close(fig)
        bode_plot(det,  frequencies, results_folder, file_name=filename + "_GNC_det_Bode", title='Bode plot of '+str(offset)+r' + det[$I + L(j \omega)$] over '+str(len(frequencies))+' frequencies', style="solid", save_pickle=save_pickle, legend=None)

    # Save the results
    if save_results:
        np.savetxt(results_folder + '\\' + filename + '_GNC_det.txt', np.stack((frequencies,det), axis=-1), delimiter='\t',
                header="Frequency [Hz]\t"+str(offset)+"+det[I + L(s)]", comments='')
    
    if run_sigma:
        sigmas = np.linalg.svd(np.identity(L.shape[1]) + L, compute_uv=False)
        bode_plot(np.min(sigmas, axis=1), frequencies, results_folder, filename+"_GNC_sigma", title='Minimum singular value of $I + L(j\omega)$ over '+str(len(frequencies))+' frequencies',
                  legend=["\sigma_{min}"], style="solid", save_pickle=save_pickle, save_data=save_results)
        
    return dict(stable=stable_system, sigmas= sigmas if run_sigma else None, net_crossings=N)

def loci_sensitivity(right_eigenvectors, left_eigenvectors, frequencies, results_folder=None, filename='loci_sensitivity', Z=None, selected_loci=[], bus_names=[],
                     normalize=False, loci=None, Y=None, wrt_all_elements=False, make_plot=True, save_pickle=False, save_results=True):
    # Compute different sensitivities of the eigenloci of a given matrix for the frequencies of interest.
    # 1) The most basic calculation is with respect to changes in the diagonal elements of the original matrix, i.e. below diag_sensitivity[freq_idx, diag_element, locus] gives the sensitivity of the locus to the diag_element of the matrix for the frequency at freq_idx.
    # 2) If the Z matrix is provided, the sensitivity of the selected open-loop (L=Z*Y) loci with respect to changes in the elements of Y is computed by applying the chain rule
    # 3) If the Z matrix is not provided but wrt_all_elements is True, then the sensitivity of the selected loci with respect to all the elements of the original matrix is computed. This is useful to extend the closed-loop participation factors beyond the diagonal elements.
   
    if wrt_all_elements: print("Plotting or saving the sensitivity with respect to changes in all the matrix elements can be time and memory intensive.")
    
    diag_sensitivity = right_eigenvectors * left_eigenvectors.transpose(0,2,1) # Element-wise product: sensitivity of the eigenvalues (loci) to changes in the diagonal elements of the matrix
    
    if len(selected_loci) == 0:
        loci_range = range(right_eigenvectors.shape[1]) # Loop over all the loci if selected_loci is not provided and there are multiple loci
    elif right_eigenvectors.shape[1] == 1:
        loci_range = [0] # Only one locus, so only one iteration
    else:
        loci_range = selected_loci # Loop over the selected loci

    if len(bus_names) == 0: bus_names = [str(bus+1) for bus in range(right_eigenvectors.shape[1])]  # Sorted numbers if names not provided

    if normalize and loci is not None:
        for i in range(diag_sensitivity.shape[2]):
            diag_sensitivity[:, :, i] = diag_sensitivity[:, :, i] / np.abs(loci[:, None, i])  # Normalized sensitivity by the magnitude of each locus at each frequency

    for locus in loci_range:
        if make_plot:
            bode_plot(diag_sensitivity[:, :, locus], frequencies, results_folder=results_folder, file_name=filename+"_"+str(locus+1)+"_wrt_diag",
                      title='Sensitivity of the locus ' + str(locus+1) + ' with respect to the matrix diagonal elements', style="solid", save_pickle=save_pickle,
                      linear_mag=True, legend=bus_names)

        if save_results:
            np.savetxt(results_folder+'\\'+filename+"_"+str(locus+1)+'_wrt_diag.txt', np.column_stack((frequencies, diag_sensitivity[:, :, locus])), delimiter='\t',
                       header="Frequency [Hz]\t" + "\t".join(bus_names), comments='')
    
        if Z is not None:
            # Compute the sensitivity of the open-loop locus to changes in each (i, j) element of Y by chain rule  
            B = np.transpose(right_eigenvectors[:, :, locus][:, :, None] @ left_eigenvectors[:, locus, :][:, None, :], axes=(0,2,1)) 
            S = Z.transpose(0,2,1) @ B 

            if normalize and loci is not None:
                locus_abs = np.abs(loci[:, locus])
                S = S / locus_abs[:,None,None] # Normalized sensitivity w.r.t. the magnitude of the locus

            if normalize and Y is not None:
                S = S * np.abs(Y)  # Normalized sensitivity w.r.t. the magnitude of each element of Y

            if make_plot:
                if wrt_all_elements:
                    bode_plot(S, frequencies, results_folder=results_folder, file_name=filename+"_"+str(locus+1)+"_wrt_Y",style="solid", save_pickle=save_pickle,
                              title='Sensitivity of locus ' + str(locus+1) + ' with respect to the elements of Y', linear_mag=True, legend=bus_names)
                bode_plot(np.diagonal(S,axis1=1,axis2=2), frequencies, results_folder=results_folder, file_name=filename+"_"+str(locus+1)+"_wrt_Y_diag",
                          title='Sensitivity of locus ' + str(locus+1) + ' with respect to the diagonal elements of Y', style="solid",
                          save_pickle=save_pickle, linear_mag=True, legend=bus_names)
            if save_results:
                if wrt_all_elements:
                    np.savetxt(results_folder+'\\'+filename+"_"+str(locus+1)+"_wrt_Y.txt", np.column_stack((frequencies, S.reshape(S.shape[:-2] + (-1,), order='C'))), delimiter='\t',
                               header="Frequency [Hz]\t" + "\t".join(bus_names), comments='')

                np.savetxt(results_folder+'\\'+filename+"_"+str(locus+1)+"_wrt_Y_diag.txt", np.column_stack((frequencies, np.diagonal(S,axis1=1,axis2=2))), delimiter='\t',
                           header="Frequency [Hz]\t" + "\t".join(bus_names), comments='')
        
        if wrt_all_elements:
            # Compute the matrix product of right and left eigenvectors for the selected locus to get its sensitivity to changes in all the elements of the original matrix
            if Z is not None:
                PFs = B.transpose(0,2,1) # Re-use the matrix product
            else:
                PFs = right_eigenvectors[:, :, locus][:, :, None] @ left_eigenvectors[:, locus, :][:, None, :] # Compute the matrix product

            if make_plot:
                bode_plot(PFs, frequencies, results_folder=results_folder, file_name=filename+"_"+str(locus+1),style="solid", save_pickle=save_pickle,
                          title='Sensitivity of locus ' + str(locus+1) + ' with respect to all matrix elements', linear_mag=True, legend=bus_names)
            if save_results:
                np.savetxt(results_folder+'\\'+filename+"_"+str(locus+1)+".txt", np.column_stack((frequencies, PFs.reshape(PFs.shape[:-2] + (-1,), order='C'))), delimiter='\t',
                           header="Frequency [Hz]\t" + "\t".join(bus_names), comments='')

    # Plotting and saving the sensitivity of all the loci to changes in all the diagonal elements of the matrix
    if len(selected_loci) == 0:
        warn("Plotting or saving the sensitivity of all the loci to changes in all the diagonal elements can be time and memory intensive. Consider selecting a specific loci.")
        if make_plot:
            bode_plot(diag_sensitivity, frequencies, results_folder=results_folder, file_name=filename+"_wrt_diag",
                    title='Sensitivity of the loci with respect to the matrix diagonal elements', style="solid", save_pickle=save_pickle,
                    linear_mag=True, legend=bus_names)
        if save_results:
            np.savetxt(results_folder+'\\'+filename+'_wrt_diag.txt', np.column_stack((frequencies, diag_sensitivity if diag_sensitivity.ndim < 3 else diag_sensitivity.reshape(diag_sensitivity.shape[:-2] + (-1,), order='C'))),
                       delimiter='\t', header="Frequency [Hz]\t" + "\t".join(bus_names), comments='')

    return diag_sensitivity

def unstable_frequency(locus, frequencies, results_folder=None, filename='unstable_frequency', order_maxima=3, make_plot=True, save_pickle=False, open_loop=True, max_zeta_abs=0.707, diff_check=True):
    G = 1/(1+locus) if open_loop else locus  # Built the function to be analised containing the closed-loop poles
    G_mag = np.abs(G)
    G_ph = np.unwrap(np.angle(G)) # Unwrapped phase to avoid discontinuities
    critical_points = argrelmax(G_mag, order=order_maxima)[0] # Peaks = potential unstable modes
    unstable_freqs = [] # List to store the frequencies of the unstable modes
    unstable_dampings = [] # List to store the approximate damping ratios of the unstable modes
    for critical_point in critical_points:
        # Three-point formula to approximate the derivate of the phase around the local maxima of the magnitude
        d1 = frequencies[critical_point] - frequencies[critical_point-1]
        d2 = frequencies[critical_point+1] - frequencies[critical_point]
        c_minus = -d2 / ( d1*(d1+d2) )
        c_zero  =  (d2 - d1) / ( d1*d2 )
        c_plus  =  d1 / ( d2*(d1+d2) )
        dtheta_df = c_minus*G_ph[critical_point-1] + c_zero*G_ph[critical_point] + c_plus*G_ph[critical_point+1]
        zeta = -1/(dtheta_df*frequencies[critical_point]) # Damping ratio approximation based on the phase derivative

        # Optional: check that the first order approximations also hold around the critical point
        dtheta_df_left = (G_ph[critical_point] - G_ph[critical_point-1]) / (frequencies[critical_point] - frequencies[critical_point-1])
        dtheta_df_right = (G_ph[critical_point+1] - G_ph[critical_point]) / (frequencies[critical_point+1] - frequencies[critical_point])
        if diff_check and (np.sign(dtheta_df_left) != np.sign(dtheta_df) or np.sign(dtheta_df_right) != np.sign(dtheta_df)):
            # Show a warning if the three derivates do not match
            print(f" Warning: Uncertain phase shift analysis at {frequencies[critical_point]} Hz. Try increasing the frequency resolution.")
        elif abs(zeta) < max_zeta_abs and not np.signbit(dtheta_df): # Only consider unstable modes with a positive phase derivative and a low damping ratio (otherwise the second-order approximation might not be valid)
            unstable_freqs.append(frequencies[critical_point]) # Unstable frequency if the phase derivative is positive
            unstable_dampings.append(zeta) # Estimated damping ratio using the phase derivative

    if make_plot and results_folder is not None:
        fig_bode, ax_bode = bode_plot(Y=G, frequencies=frequencies, results_folder=None, style="solid", legend=None, return_plot=True,
                                      title=r"Unstable mode identification: $1/(1+\lambda)$ over "+str(len(frequencies))+' frequencies')
        for freq in unstable_freqs:
            ax_bode[0].axvline(x=freq, color='red', linestyle=':', linewidth=1, label='_nolegend_')
            ax_bode[0].text(freq, 0.20, str(round(freq,2)), color='r', ha='right', va='bottom', rotation=90, transform=ax_bode[0].get_xaxis_transform())
            ax_bode[1].axvline(x=freq, color='red', linestyle=':', linewidth=1, label='_nolegend_')
            ax_bode[1].text(freq, 0.20, str(round(freq,2)), color='r', ha='right', va='bottom', rotation=90, transform=ax_bode[1].get_xaxis_transform())

        fig_bode.savefig(results_folder + '\\' + filename + ".pdf", format="pdf", bbox_inches="tight")
        if save_pickle:
            with open(results_folder + '\\' + filename + ".pickle", 'wb') as f: pickle.dump(fig_bode, f)
        plt.close(fig_bode)
        
    return unstable_freqs, unstable_dampings
             
def mode_estimation(G, omegas, zeta0=0.05, omega0=None, extra_poles=0, weight=None, reg=1e-3, enable_d=False, verbose=False):
    if omega0 is None: omega0 = np.sum(omegas * np.abs(G)) / np.sum(np.abs(G)) # Initial guess of the frequency 
    p0 = -zeta0*omega0 + 1j*omega0 # Initial guess of the complex conjugate pole 
    r0 = G[np.argmin(np.abs(omegas - omega0))] * (1j*omega0 - p0) # Initial residue guess
    d0 = (G[0] + G[-1])/2 # Initial guess of the direct term based on the average of the first and last points of G (baseline estimation)
    # "parameters0" containes the real and imaginary parts of each poles and residues, as well as the direct term at the last entry
    parameters0 = np.array([p0.real, p0.imag, r0.real, r0.imag]) # Complex conjugate mode first
    # Add as many extra poles as specified by extra_poles, each with its own initial guess of pole and residue
    poles_guess = -np.logspace(np.log10(omegas[0]), np.log10(omegas[-1]), extra_poles) # Initial guess of the pole
    for k in range(extra_poles):
        pk = poles_guess[k] # Initial guess of the pole
        rk = 0.0 # Initial guess of the residue
        parameters0 = np.append(parameters0, [pk.real, pk.imag, rk.real, rk.imag]) # Real and imaginary parts of the pole and residue
    parameters0 = np.append(parameters0, d0.real) # Add the direct term
    
    # Define the frequency-wise weights for the least-squares optimization
    if weight is None:
        dw = np.gradient(omegas)
        weight = dw / np.mean(dw)
    # Perform least-squares optimization
    res = least_squares(lsq_residuals, parameters0, args=(omegas, G, extra_poles, weight, reg, enable_d) ) 
    parameters_opt = res.x
    if not enable_d: parameters_opt[-1] = 0.0 # Set the direct term to zero if it is not enabled in the optimization

    if verbose:
        print(f"Modal estimation results with a final cost of {res.cost:.4e} and success status {res.success}:")
        for modes in range(extra_poles+1):
            p = parameters_opt[4*modes] + 1j*parameters_opt[4*modes+1]
            r = parameters_opt[4*modes+2] + 1j*parameters_opt[4*modes+3]
            print(f" Pole {modes+1}: {p.real:.4e} + {p.imag:.4e}j rad/s, and residue {r.real:.4e} + {r.imag:.4e}j")
        if enable_d: print(f" Direct term: {parameters_opt[-1]:.4e}\n")

    return parameters_opt

def lsq_residuals(parameters, omegas, G, extra_poles=1, weight=None, reg=1e-3, enable_d=False):
    # Simple function to solve the least-squares optimization problem for rational fitting in the mode_estimation function.
    # It computes the residuals between the given G and the estimated G based on the given pole and residue parameters.
    # The first four elements in "parameters" correspond to the real and imaginary parts of the complex conjugate pole and its residue,
    # The rest of the "parameters" entries correspond to the real and imaginary parts of each extra pole and residue, and the direct term at the end.
    # Frequency-wise weights can be applied to the residuals to prioritize certain frequency ranges in the optimization.
    # An optional regularization term "reg" can be used to promote small absolute real parts of the poles.
    # The "enable_d" flag can be used to include (True) or exclude (False) the direct term in the optimization.

    s = 1j * omegas
    # Second order system with at least one pair of complex conjugate poles
    sigma, omega = parameters[0], parameters[1]
    r0 = parameters[2] + 1j*parameters[3]
    G_est = r0/(s-(sigma+1j*omega)) + np.conj(r0)/(s-(sigma-1j*omega)) # Complex conjugate mode

    # Optional extra poles and collection of real parts for optional regularization
    p_re = [parameters[0]]  # Sigma of oscillatory mode
    idx = 4 # Four parameters already defined
    for k in range(extra_poles):
        p = parameters[idx] + 1j*parameters[idx+1]
        r = parameters[idx+2] + 1j*parameters[idx+3]
        p_re.append(parameters[idx]) # Real part of the pole for regularization
        G_est += r / (s - p) # Rational function with the given poles and residues at s
        idx += 4
    G_est += parameters[-1] if enable_d else 0.0 # Add the direct term

    if weight is None:
        weight = 1.0
    else:
        weight = np.sqrt(weight).repeat(2)
    e = weight * np.hstack([(G - G_est).real,  (G - G_est).imag]) # Error vector

    p_re = np.array(p_re) # Pole real-part regularization (optional)
    reg_term = np.sqrt(reg) * p_re

    return np.hstack([e, reg_term])

nyquist.__doc__ = """
Stability assessment based the Generalized Nyquist Criteria (GNC): eigenvalue decomponsition (EVD) of the open-loop (minor-loop) matrix over the frequency.

The GNC can be stated as follows considering a contour along the imaginary axis and around the whole Right-Half Plane (RHP) avoiding open-loop poles:
N: Net number of clockwise encirclements by the open-loop eigenloci =  clockwise -  counter-clockwise
P: Number of RHP poles of the open-loop system
Z: number of RHP poles of the closed-loop system
Argument principle states N = Z - P over the closed contour. Therefore, if Z = N + P > 0, then the closed-loop system has RHP poles and it is unstable.
Assuming standalone-stable subsystems means P = 0, and thus N = 0 implies stability while N > 0 implies instability.
If the subsystems are standalone unstable, P > 0, then possibly N < 0.
The interested reader is referred to S. Skogestad and I. Postlethwaite, "Multivariable Feedback Control: Analysis and Design", Wiley, 2005 for a more detailed explanation.

The function computes the eigenvalues of L at every frequency, plots the eigenloci and counts the number of clockwise and counter-clockwise encirclements of (-1,0j).
The EVD of L is saved as filename_GNC.txt and its plot is saved as filename_GNC.pdf. In addition, the Bode plot 1/(1+locus) is generated to aid in the determination of the unstable frequencies.
Identification of unstable frequencies is performed when unstable_frequency=True based on the local maxima of the magnitude of 1/(1+locus) and the sign of its phase derivative at said frequencies.
The number of reported unstable frequencies is equal to the number of net encirclements by each eigenlocus. To minimize numerical errors the candidate frequencies are first sorted by their approximate damping ratio.
The indentations argument specifies frequencies at which indentations in the Nyquist contour are performed so as to avoid open-loop poles in the imaginary axis.
Lastly, the sensitivity of the critical loci with respect to changes in each element of Y is computed by applying the chain rule if Z in L=Z*Y is provided.

Required arguments
        L                   (numpy ndarray of complex double) Minor loop gain (transfer matrix) for different frequencies.
        frequencies         (numpy array) Frequencies over which L is computed [Hz].
        results_folder      (str) Absolute path where the results are to be stored. If it does not exist, it is created.
        filename            (str) Name root of the results output files.        
      
Optional arguments
        verbose                 (bool) Bool flag to show detailed GNC application information, such as the number of counter clock-wise (CCW) and clock-wise (CC) encirclements of the critical point.
        check_conditioning      (bool) Bool flag to discard values with poor numerical conditioning of L.
        indentations            (list of double) Frequencies [Hz] at which indentations around open-loop poles are performed.
        condition_number_th     (double) Condition number threshold of L above which the data is ignored.
                                This threshold can be set based on the expected input error and maximum acceptable ouTput error.
                                For example, for a relative output error <= 0.01 considering a relative input error of 5e-9, the condition number threshold can be set to 0.01/5e-9 (default value).
        save_pickle             (bool) Bool flag to save the generated plots as pickle objects in addition to pdf files. Default = False.
        save_results            (bool) Bool flag to save the results in a text file. Default = True.
        run_sensitivity         (bool) Bool flag to run the sensitivity analysis of the critical loci with respect to changes in each element of Y with L=Z*Y. Default = False.
        Z                       (numpy ndarray of complex double) If provided, the sensitivity of the critical loci with respect to changes in each element of Y is computed by applying the chain rule. Default = None.
                                The critical loci are selected as those showing encirclements of (-1,0j) or that closest to the critical point (-1,0j).
        Y                       (numpy ndarray of complex double) If provided together with Z, the sensitivity is normalized. Default = None.
        bus_names               (list of str) List of bus names to be used in the sensitivity analysis. Default = empty list, which results in the use of sorted numbers as bus names.
        unstable_frequency      (bool) Bool flag to run the unstable frequency identification based on the local maxima of the magnitude of 1/(1+locus) and the sign of its phase derivative. Default = False.
        order_maxima            (int) Points on each side of each local maximum used for the comparison and maxima identification in the unstable frequency identification, i.e., the 'order' argument of the 'argrelmax' function. Default = 4.
        modal_estimation        (bool) Bool flag to run the modal estimation of the dominant modes based on least-squares rational fitting around the unstable modes. Default = False.
        verbose_modal_estimation (bool) Bool flag to show detailed information about the modal estimation results, such as the estimated poles and residues. Default = False.
        extra_poles             (int) Number of extra poles to be added in the least-squares rational fitting for modal estimation. Default = 0, which results in a second-order fit with one pair of complex conjugate poles.
        samples_fitting         (int) Number of frequency samples to be used in the least-squares rational fitting for modal estimation. Default = 12, which results in a fitting around the unstable modes with 6 points on each side.
        run_sigma               (bool) Bool flag to run the sigma analysis for the identification of the frequencies at which the system is most prone to de-stabilization. Default = False.

Returns
        Bool flag indicating closed-loop or interconnected stability: True means stable.

"""

nyquist_det.__doc__ = """
Stability assessment based on based the Generalized Nyquist Criteria (GNC) by counting the encirclements of the critical point by the determinant of I + L over the frequency.

Theorem 4.14 in S. Skogestad and I. Postlethwaite, "Multivariable Feedback Control: Analysis and Design", Wiley (2005), assuming standalone stable subsystems,
i.e. L(s) does not have any open-loop unstable poles, the stability conditions are:
a) Zero net number of clockwise encirclements of (0,j0) by det[I+L(s)] as s travels the imaginary axis from 0 to +j*infinity avoiding the pure imaginary poles of L
b) No crossings of the orgin by by det[I + L(s)] as s travels the imaginary axis from 0 to +j*infinity
Since only real systems are considered in practice, it implies that det[I + L(+j*infinity)] settles on the real axis and just a large enough frequency can be used to approximate det[I + L(+j*infinity)].

The interested reader is referred to S. Skogestad and I. Postlethwaite, "Multivariable Feedback Control: Analysis and Design", Wiley, 2005 for a more detailed explanation.

The function computes the eigenvalues of L at every frequency, plots the eigenloci and counts the number of clockwise and counter-clockwise encirclements of (-1,0j).
The EVD of L is saved as filename_GNC.txt and its plot is saved as filename_GNC.pdf.
The indentations argument specifies frequencies at which indentations in the Nyquist contour are performed so as to avoid open-loop poles in the imaginary axis.

Required arguments
        L                   (numpy ndarray of complex double) Minor loop gain (transfer matrix) for different frequencies.
        frequencies         (numpy array) Frequencies over which L is computed [Hz].
        results_folder      (str) Absolute path where the results are to be stored. If it does not exist, it is created.
        filename            (str) Name root of the results output files.        
      
Optional arguments
        indentations        (list of double) Frequencies [Hz] at which indentations around open-loop poles are performed.
        verbose             Bool flag to show detailed GNC application information, such as the number of counter clock-wise (CCW) and clock-wise (CC) encirclements of the critical point. Default = True.
        draw_arrows         Bool flag to draw arrows on the direction of the Nyquist plot. Default = True.
        offset              (double) The offset parameter can be used to shift the critical point on the real axis instead of (0,j0). Default = 0.
                            For example, offset = -1.0 defines the critical point as (-1,j0) as when applying the GNC via the eigenvalue loci of L.
        show_plot           (bool) to show the plot interactively. Default = False.
        save_pickle         Bool flag to save the generated plots as pickle objects in addition to pdf files. Default = False.
        save_results        (bool) Bool flag to save the results in a text file. Default = True.

Returns
        Bool flag indicating closed-loop or interconnected stability: True means stable.

"""

stability_analysis.__doc__ = """
Performs a small-signal analysis based on the scanned matrices of the network defined by the topology file over the frequency and stores the results in the specified results folder.
Firstly, the function reads the individual admittance matrices from the specified topology file. Then it builds the edge and node (block-diagonal) matrices, including the necessary rotations to a common reference frame.
Finally, it calculates different metrics used to assess the small-signal dynamics of the system:
- Generalized Nyquist Criteria (GNC) based on the eigenvalue decomposition (EVD) and the determinant of the open-loop (minor-loop) matrix L computed as the product of the edge and node matrices.
- Eigenvalue decomposition (EVD) of the closed-loop impedance matrix for oscillation modes identification.
- Bus participation factors computation at the frequency of largest current amplification.
- Small-gain theorem based on the singular value decomposition (SVD) of the edge and node matrices:: used to identify risk-free frequencies.
- Passivity assessment based on the minimum singular value of the Hermitian part of the different system matrix: used to identify risk-free frequencies.
The user can choose which of the above analyses to perform by setting the corresponding boolean flags.
The results are stored in the specified results folder as text files and pdf plots based on the boolean arguments save_results and make_plot.

Required arguments
        topology            (str) Absolute path to the topology file defining the system the subsystems interconnections.
        results_folder      (str) Absolute path where the matrices and powerflow files are stored.
        file_root           (str) Name root of the files used in the analysis

Optional arguments
        indentations            (list of double) Frequencies [Hz] at which indentations around open-loop poles are performed in the GNC.
        node_blocks             (list of strings) List of strings where each entry is "BlockName-side" corresponding to the node matrix components. Default = None, which results in the automated identification of the blocks as per the read_admittance function.
        check_conditioning      (bool) Bool flag to discard values with poor numerical conditioning of the system matrices.
        condition_number_th     (double) Condition number threshold of the system matrices above which the data is ignored.
                                This threshold can be set based on the expected input error and maximum acceptable ouTput error.
        make_plot               (bool) Bool flag to enable/disable the generation of pdf plot files.
        save_pickle             (bool) Bool flag to save the generated plots as pickle objects in addition to pdf files. Default = False.
        save_results            (bool) Bool flag to save the results in a text file. Default = True.
        save_Y                  (bool) Bool flag to save the system admittance matrices in text files. Default = True.
        save_loop_gain          (bool) Bool flag to save the loop gain matrix as a text file. Default = True.
        verbose                 (bool) Bool flag to show detailed analysis information, such as the crossings in the GNC and the participation factors of the dominant mode. Default = True.
        reference_buses         (list of str) List of AC bus names (e.g. Z-tool scan block names) which define the reference angles for each area. Default = None, which results in the use of the first block of each area as reference.
        relative_angles         (bool) Bool flag to use the relative angles for the matrix rotations. If set to True, the relative angle between each bus and the reference buses is used; otherwise the angles in the text file are used directly. Default = True.
        run_nyquist             (bool) Bool flag to run the Generalized Nyquist Criteria (GNC) based on the eigenvalues of the open-loop matrix L. Default = True.
        run_nyquist_det         (bool) Bool flag to run the determinant-based Nyquist stability assessment. Default = False.
        run_EVD                 (bool) Bool flag to run the eigenvalue decomposition (EVD) of the closed-loop system for oscillation modes identification and bus participation factors computation. Default = True.
        run_EVD_PFs             (bool) Bool flag to compute the bus participation factors (PFs) of the largest magnitude modal impedance magnitude. Default = True.
        run_EVD_PFs_extended    (bool) Bool flag to compute the sensitivity of the critical locus to changes all elements of the original matrix. This is the extension of the bus PFs beyond the diagonal elements. Default = False.
        run_passivity           (bool) Bool flag to run the passivity assessment of the system matrices. Default = True.
        run_small_gain          (bool) Bool flag to run the small-gain theorem based on the system matrices. Default = True.
        run_GNC_sensitivity     (bool) Bool flag to run the sensitivity analysis of the critical loci when applying the GNC with respect to the diagonal elements of L as well as with respect to Y with L=Z*Y. Default = False.
        normalize_GNC_sensitivity (bool) Bool flag to normalize the sensitivity of the critical loci in the GNC with respect to changes in each element of Y with L=Z*Y. The normalization is with respect to each admittance magnitude and locus magnitude. Default = False.
        order_maxima            (int) Points on each side of each local maximum used for the comparison and maxima identification in the PMD criterion and unstable frequency identification, i.e., the 'order' argument of the 'argrelmax' function. Default = 4.
        modal_estimation_nyquist(bool) Bool flag to run the modal estimation of the unstable modes based on least-squares rational fitting around the peaks of 1/(1+lambda) where lambda are the critical eigenvalues of the open-loop matrix. Default = False.
        modal_estimation_EVD    (bool) Bool flag to run the modal estimation of the dominant closed-loop modes based on least-squares rational fitting around the peaks of the closed-loop modal impedance. Default = False.
        extra_poles             (int) Number of extra poles to be added in the least-squares rational fitting for modal estimation. Default = 0, which results in a second-order fit with one pair of complex conjugate poles.
        samples_fitting         (int) Number of frequency samples to be used in the least-squares rational fitting for modal estimation. Default = 12, which results in a fitting around the unstable mode with 6 points on each side.
        run_sigma               (bool) Bool flag to run the sigma analysis for the identification of the frequencies at which the system is most sensitive to perturbations. Default = False.
        Ibase                   (dict) Dictionary of base currents per node used to per-unitize the system matrices. The node name, i.e. scan block name, is used as key. Default = None, which results in no per-unitization.
        Vbase                   (dict) Dictionary of base voltages per node used to per-unitize the system matrices. The node name, i.e. scan block name, is used as key. Default = None, which results in no per-unitization.
        PMD_zeta_threshold      (float) Damping ratio threshold for the PMD criterion. When modal_estimation is True only modes with a damping ratio below this value are checked for stability via the PMD criterion. Default = 0.25,

Returns
        Dictionary of dynamic analysis results. The keys correspond to the different methods, e.g., "nyquist", "nyquist_det", "EVD", "passivity_index", "small_gain_index", etc., depending on the analyses performed.
        The values are the corresponding results, such as boolean flags for stability assessments, arrays of modal frequencies and damping ratios, participation factors, passivity indices, etc.
        
"""

passivity.__doc__ = """
Passivity assessment based on the minimum eigenvalue value of the Hermitian part of the system matrices over the frequency.
The passivity index is defined as the minimum eigenvalue of the Hermitian part of G, i.e., 0.5*(G + G^H), where G^H is the conjugate transpose of G.
If the system under evaluation is stable with a positive passivity index across all frequencies, then no instability can arise when interconnected to another passive system.
Required arguments
        G                   (numpy ndarray of complex double) System matrix at different frequencies.
        frequencies         (numpy array) Frequencies over which the matrix G is evaluated [Hz].
        results_folder      (str)Absolute path where the results are to be stored.
        filename            (str) Name root of the results output files.
        
Optional arguments
        variables           (list of str) Names of the block-diagonal matrices in G for block-wise analysis. Default = None.
        make_plot           (bool) Bool flag to enable/disable the generation of pdf plot files.
        save_pickle         (bool) Bool flag to save the generated plots as pickle objects in addition to pdf files. Default = False.
        save_results        (bool) Bool flag to save the results in a text file. Default = True.

Returns
        (numpy array) Passivity index value over the frequency.
"""

EVD.__doc__ = """
Eigenvalue decomposition (EVD) of the closed-loop system impedance matrix over the frequency for oscillation modes identification and bus participation factors (PF) computation.
The function is based on the developements presented in https://doi.org/10.1109/TPWRD.2004.834856 and https://doi.org/10.1016/j.ijepes.2023.108957
If the provided matrix G is the closed-loop admittance matrix, set Z_closedloop = True to avoid needless matrix inversion.

Required arguments
        G                   (numpy ndarray of complex double) Closed-loop system matrix at different frequencies.
        frequencies         (numpy array) Frequencies over which the matrix G is evaluated [Hz].
        results_folder      (str) Absolute path where the results are to be stored.
        filename            (str) Name root of the results output files.

Optional arguments
        bus_names          (list of str) Names of the buses in the system for PF computation. Default = None.
        verbose            (bool) Bool flag to show detailed analysis information.
        Z_closedloop       (bool) Bool flag indicating if G is the closed-loop impedance matrix. This can be used to avoid the inversion of G = Ynode + Yedge. Default = True.
        make_plot          (bool) Bool flag to enable/disable the generation of pdf plot files.
        save_pickle        (bool) Bool flag to save the generated plots as pickle objects in addition to pdf files. Default = False.
        save_results       (bool) Bool flag to save the results in a text file. Default = True.
        PFs                (bool) Bool flag to compute the bus participation factors (PFs) of the largest magnitude modal impedance magnitude. Default = True.
        PFs_extended       (bool) Bool flag to compute the sensitivity of the critical locus to changes the elements of the original matrix. This is the extension of the bus PFs beyond the diagonal elements. Default = False.
        run_PMD            (bool) Bool flag to run the positive mode damping (PMD) criterion. Default = False. Find more information on the PMD on the following paper
                           Luis Orellana, et al. "Study of black-box models and participation factors for the Positive-Mode Damping stability criterion",2023 https://doi.org/10.1016/j.ijepes.2023.108957
        order_maxima       (int) Points on each side of each local maximum used for the comparison and maxima identification in the PMD criterion, i.e., the 'order' argument of the 'argrelmax' function. Default = 4.
        modal_estimation   (bool) Bool flag to run the modal estimation of the dominant modes based on least-squares rational fitting around the mode. Default = False.
        extra_poles        (int) Number of extra poles to be added in the least-squares rational fitting for modal estimation. Default = 0, which results in a second-order fit with one pair of complex conjugate poles.
        samples_fitting    (int) Number of frequency samples to be used in the least-squares rational fitting for modal estimation. Default = 12, which results in a fitting around the unstable mode with 6 points on each side.
        run_sigma          (bool) Bool flag to run the sigma analysis for the identification of the frequencies at which the system is most sensitive to perturbations. Default = False.
        PMD_zeta_threshold (float) Damping ratio threshold for the PMD criterion when modal_estimation is True: only modes with a damping ratio below this value are checked for stability via the PMD criterion. Default = 0.707.

Returns
        Dictionary including main results such as stability assessment by the PMD criterion (True means stable), the modal impedances and participation factors, the fitted modes if modal_estimation is True, and other relevant metrics.
                            
"""

small_gain.__doc__ = """
Applies a conservative version of the small-gain theorem verifying if |L| = |G1*G2| <= |G1|*|G2| < 1 holds by plotting
the maximum singular value of G1, G2 and G1*G2 over the frequency as well as the unitary gain line.
If the unitary gain line is not crossed by the maximum singular value of G1*G2, then no instability can arise at said frequency.
To visually check this, the plot of 1/|G1| is compared with that of |G2|, since if |G2| < 1/|G1|, then |G1*G2| < 1.
Therefore, the Bode plot of |G2| should be below 1/|G1| to guarantee stability.
If G2 is block-diagonal, the plot of the maximum singular value of each diagonal block in G2 is also computed.

Required arguments
        G2                  (numpy ndarray of complex double) System matrix at different frequencies. Possibly block-diagonal.
        frequencies         (numpy array) Frequencies over which the matrix G is evaluated [Hz].
        results_folder      (str) Absolute path where the results are to be stored.
        filename            (str) Name root of the results output files.
 
Optional arguments
        G1                  (numpy ndarray of complex double) System matrix at different frequencies.
        variables           (list of str) Names of the block-diagonal matrices in G2 for block-wise analysis. Default = None.
        make_plot           (bool) Bool flag to enable/disable the generation of pdf plot files.
        save_pickle         (bool) Bool flag to save the generated plots as pickle objects in addition to pdf files. Default = False.
        save_results        (bool) Bool flag to save the results in a text file. Default = True.

Returns
        (numpy array) Array of the frequency-wise product of the maximum singular value of G1 and the maximum singular value of G2.
"""

loci_sensitivity.__doc__ = """
Computation of the sensitivity of the eigenloci of a given matrix for the frequencies of interest. The results and their interpretation depend on the original matrix and arguments provided to the function:
1) The most basic sensitivity calculation is with respect to changes in the diagonal elements of the original matrix, i.e. each entry [freq_idx, diag_element, locus] gives the sensitivity of the locus to the diag_element of the matrix for the frequency at freq_idx.
2) If the Z matrix is provided, the sensitivity of the selected open-loop (L=Z*Y) loci with respect to changes in the elements of Y is computed by applying the chain rule.
3) If the Z matrix is not provided but wrt_all_elements is True, then the sensitivity of the selected loci with respect to all the elements of the original matrix is computed. This is useful to extend the closed-loop participation factors beyond the diagonal elements.
Note that plotting or saving the sensitivity of all the loci to changes in all the diagonal elements can be time and memory intensive. Consider specifying selected_loci before calling the function.
Lastly, the results can be normalized with respect to the magnitude of each locus and/or the magnitude of the elements of Y by providing the corresponding arguments.

Required arguments
        right_eigenvectors  (numpy ndarray of complex double) Right eigenvectors of the matrix for the frequencies of interest.
        left_eigenvectors   (numpy ndarray of complex double) Left eigenvectors of the matrix for the frequencies of interest.
        frequencies         (numpy array or list) Frequencies [Hz].
        results_folder      (str) Absolute path where the results are to be stored.
        filename            (str) Name root of the results output files.


Optional arguments
        Z                   (numpy ndarray of complex double) If provided, the sensitivity of the open-loop locus to changes in each (i, j) element of Y is computed by applying the chain rule. Default = None.
        selected_loci       (list of int) List of loci indexes to compute the sensitivity for. Default = [], which results in the computation of the sensitivity for all the loci with respect to the diagonal elements.
        bus_names           (list of str) Names of the buses in the system for labeling the sensitivity results. Default = [] which generates numeric labels.
        normalize           (bool) Bool flag to normalize the sensitivity of each locus with respect to the locus magnitude. Default = False (no normalization).
        loci                (numpy array) Array of shape (number of frequencies, number of loci) with all eigenloci for each frequency. It is used for normalization if normalize = True. Default = None (no locus-based normalization).
        Y                   (numpy ndarray of complex double) If provided, the sensitivity of the open-loop locus to changes in Y is normalized with respect to the Y matrix elements. Default = None (no normalization with respect to Y).
        wrt_all_elements    (bool) Bool flag to compute the sensitivity of the selected loci with respect to all the elements of the original matrix. This is useful to extend the closed-loop participation factors beyond the diagonal elements. Default = False.
        make_plot           (bool) Bool flag to enable/disable the generation of pdf plot files.
        save_pickle         (bool) Bool flag to save the generated plots as pickle objects in addition to pdf files. Default = False.
        save_results        (bool) Bool flag to save the results in a text file. Default = True.

"""

mode_estimation.__doc__ = """
Mode estimation based on least-squares optimization to fit a rational function to the provided frequency response data.
The rational fit contains a complex-conjugate pole pair, an additional number of extra poles defined by the user plus a direct term.
The function relies on the lsq_residuals function to compute the residuals of the least-squares optimization.
This lsq_residuals includes an optional regularization term on the real part of the poles to penalize too large real-parts poles.
See the book "System Identification: Theory for the User" by Lennart Ljung for more information on the rational approximation of frequency response data.

Required arguments
        G           (numpy ndarray of complex double) Frequency response data to be fitted.
        omegas      (numpy array) Frequencies [rad/s] corresponding to the frequency response data in G.
Optional arguments
        zeta0       (double) Initial guess of the damping ratio of the complex-conjugate pole pair. Default = 0.05.
        omega0      (double) Initial guess of the natural frequency of the complex-conjugate pole pair. Default = None, which uses the weighted average of the frequencies.
        extra_poles (int) Number of extra poles to be added to the rational fit. Default = 1.
        weight      (numpy array) Frequency-wise weights for the least-squares optimization. Default = None, which uses the normalized frequency spacing as weights.
        reg         (double) Regularization parameter for penalizing too large real parts of the poles. Default = 1e-3.
        enable_d    (bool) Bool flag to include the direct term in the optimization. Default = False.
        verbose     (bool) Bool flag to show the fitting results: the estimated poles and residues. Default = True.
Returns
        parameters_opt (numpy array) Array containing the optimized parameters of the rational fit.
                        The first four entries correspond to the real and imaginary parts of the complex-conjugate pole pair and its residue. Next, the real and imaginary parts of the extra poles and their residues, and lastly the direct term.
"""