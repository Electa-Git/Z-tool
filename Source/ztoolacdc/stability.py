"""
This program contains several functions for frequency-domain stability analysis including:
    1) Generalized Nyquist Criterion (GNC) application to determine system stability: via eigenvalue decomposition and via the determinant
    2) Eigenvalue Decomposition (EVD) of the closed-loop matrix to determine oscillatory modes and bus participation factors
    3) Passivity index (for the application of the passivity theorem) and singular value decomposition (small-gain theorem) of target matrices
    4) A main stability_analysis function to apply all the previously described to a specific system

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

__all__ = ['stability_analysis','passivity','nyquist','small_gain','EVD','nyquist_det']

import numpy as np
from scipy.optimize import linear_sum_assignment
import matplotlib.pyplot as plt
from .read_admittance import read_admittance
from os import path, makedirs
import pickle

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

def stability_analysis(topology=None, results_folder=None, file_root=None, check_conditioning=False, condition_number_th=10e6):
    # 0) Read the terminal angle information
    block_area_angle = []  # List for each block containing a list as [bus/block name, area_id, terminal angle in rad]
    with open(results_folder+r'\\'+file_root+'_angles.txt', 'r') as file:
        next(file)  # First line contains the header
        for line in file:
            content = [str(line.strip().split("\t")[0]),int(line.strip().split("\t")[1]),float(line.strip().split("\t")[2])]
            block_area_angle.append(content)
    areas = set([bus[1] for bus in block_area_angle])
    # Define a reference bus for each area
    reference_angle = {i: -100.0 for i in areas}  # Unrealistic intial values
    for bus in block_area_angle:
        if len(areas) == 1:
            reference_angle[0] = float(bus[2])
            reference_angle[1] = float(bus[2])
        else:
            if reference_angle[bus[1]] == -100.0:
                reference_angle[bus[1]] = float(bus[2])  # If no angle is already defined, use the first one from that area
                # print("The reference for area",bus[1],'is bus',bus[0])

    # This function loads and builds the edge and nodal admittance matrices and for stability analysis
    # 1) Read the topology matrix and extract block names
    Ytopology = np.loadtxt(topology, skiprows=1, comments=["#", "%", "!"])
    # nameA-1 nameA-2 nameB-1 nameB-2 ... x nameA-1 nameA-2 nameB-1 ...
    # 0 means no interconnection, 1 means connection between the edges: diagonals are single-sided / shunt
    with open(topology, 'r') as f:
        block_names_Y = f.readline().strip('\n').split("\t")

    # 2) Read the admittance files based on the topology file
    admittances = []  # List containing the admittance objects
    bus_names = []  # List of list with block (bus) names
    # Create the undirected graph: adjacent matrix but diagonals can be 1
    g = Graph(len(Ytopology))
    for row, name in enumerate(block_names_Y):
        for col, edge in enumerate(Ytopology[row]):
            if int(edge) == 1: g.addEdge(row, col)
    # Obtain the connected components of the graph
    cc = g.connectedComponents()  # List of lists with blocks positions connected
    # For every group of connected buses, read and add the admittance matrix
    for buses in cc:
        bus_names.append([block_names_Y[bus] for bus in buses])
        admittances.append(read_admittance(path=results_folder, involved_blocks=bus_names[-1], file_root=file_root))

    # 3) Update bus_names to be ordered as the variables in the individual admittance matrices and build the node matrix
    node_matrix = []  # Create the node matrix with the active components (block diagonal)
    node_variables = []  # List of variable names = the current/voltage vectors
    edge_aux_matrix = []  # Auxiliary edge matrix (block diagonal)
    edge_aux_variables = []  # The order of this aux matrix is different from the order of the node matrix
    for idx, y in enumerate(admittances):
        bus_names[idx] = y.blocks
        # print(bus_names[idx],y.y_type,"- Node:",y.node)
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
    for idx, yedge in enumerate(edge_aux_matrix):
        # # Eliminate too small elements: related to PSCAD accuracy when the topology is not enforced
        # for col in range(np.size(yedge.y, 1)):
        #     for row in range(np.size(yedge.y, 2)):
        #         if max(abs(yedge.y[:, row, col])) < 1e-6:  # The threshold is system and time-step dependent
        #             yedge.y[:, row, col] = np.zeros(yedge.y[:, row, col].shape)
        Yedge_aux[:, y_edge_idx:y_edge_idx + len(yedge.vars), y_edge_idx:y_edge_idx + len(yedge.vars)] = yedge.y
        y_edge_idx = y_edge_idx + len(yedge.vars)  # Update the matrix index for the next admittance block

    # Sparsity plot for verification at the lowest frequency
    np.seterr(divide='ignore')
    plt.imshow(20 * np.log10(np.abs(Yedge_aux[0, :, :])), cmap='spring', interpolation='nearest')
    plt.colorbar()
    plt.xticks(ticks=np.arange(0, len(edge_aux_variables), step=1), labels=[])
    plt.yticks(ticks=np.arange(0, len(edge_aux_variables), step=1), labels=edge_aux_variables)
    plt.grid(visible=True, which='minor', alpha=0.3, color='k', linestyle='-', linewidth=0.5)
    plt.title('Auxiliary edge admittance matrix at '+format(frequencies[0], '.1f')+' Hz')
    plt.savefig(results_folder + '\\' + file_root + "_Edge_aux.pdf", format="pdf", bbox_inches="tight")
    with open(results_folder + '\\' + file_root + "_Edge_aux.pickle", 'wb') as f:
        pickle.dump(plt.gcf(), f)
    plt.close()

    # Re-sort the edge matrix acording to the node matrix variables
    Yedge = Yedge_aux[:,:,edge_ordering]  # Sort the columns
    Yedge = Yedge[:,edge_ordering,:]  # Sort the rows

    # Plot the sorted edge admittance for verification at the lowest frequency
    plt.imshow(20 * np.log10(np.abs(Yedge[0, :, :])), cmap='spring', interpolation='nearest')
    plt.colorbar()
    plt.xticks(ticks=np.arange(0, len(node_variables), step=1), labels=[])
    plt.yticks(ticks=np.arange(0, len(node_variables), step=1), labels=node_variables)
    #  plt.minorticks_on()
    plt.grid(visible=True, which='minor', alpha=0.3, color='k', linestyle='-', linewidth=0.5)
    plt.title('Edge admittance matrix at '+format(frequencies[0], '.1f')+' Hz')
    plt.savefig(results_folder + '\\' + file_root + "_Edge.pdf", format="pdf", bbox_inches="tight")
    with open(results_folder + '\\' + file_root + "_Edge.pickle", 'wb') as f:
        pickle.dump(plt.gcf(), f)
    plt.close()

    # Node admittance matrix and associated edge matrix
    Ynode = np.zeros((len(frequencies), len(node_variables), len(node_variables)),dtype='cdouble')  # Or dtype='csingle'
    y_node_idx = 0
    for idx, ynode in enumerate(node_matrix):
        # Rotate the matrix accordingly to the terminal angle if it involves AC variables
        if len(ynode.vars) >= 2:
            for block in block_area_angle:
                # Iterate over the different blocks and use the angle of the current matrix
                if block[0] in [bus[:-2] for bus in ynode.blocks]:
                    # The reference angle for the area of this block - bus angle of this block
                    theta = (reference_angle[block[1]] - block[2])
                    # print("Angle of",round(theta,4),"in the admittance involving",ynode.blocks)

        # Define the rotation matrix for each component: AC-side 2x2, AC/DC 3x3
        if ynode.y_type == "AC" and len(ynode.vars) == 2:
            T = np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])
            T_1 = np.array([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])
        elif ynode.y_type == "ACDC" and len(ynode.vars) == 3:
            T = np.array([[1, 0, 0],[0, np.cos(theta), -np.sin(theta)],[0, np.sin(theta), np.cos(theta)]])
            T_1 = np.array([[1, 0, 0], [0, np.cos(theta), np.sin(theta)], [0, -np.sin(theta), np.cos(theta)]])
        else:
            T = np.identity(len(ynode.vas))  # DC-side analysis only: does not need frame alignment
            T_1 = np.identity(len(ynode.vars))

        # print('Matrix',idx,", nodal matrix start idx:",y_node_idx,"and end:",y_node_idx+len(ynode.vars))
        Ynode[:, y_node_idx:y_node_idx+len(ynode.vars), y_node_idx:y_node_idx+len(ynode.vars)] = T@ynode.y@T_1
        # Equivalent to np.matmul(T,np.matmul(ynode.y,T_1))
        y_node_idx = y_node_idx + len(ynode.vars)  # Update the matrix index for the next admittance block

    # Scatter verification of the node admittance for the lowest frequency
    plt.imshow(20*np.log10(np.abs(Ynode[0, :, :])), cmap='spring', interpolation='nearest')
    plt.colorbar()
    plt.xticks(ticks=np.arange(0, len(node_variables), step=1), labels=[])
    plt.yticks(ticks=np.arange(0, len(node_variables), step=1), labels=node_variables)
    #  plt.minorticks_on()
    plt.grid(visible=True, which='minor', alpha=0.3, color='k', linestyle='-', linewidth=0.5)
    plt.title('Node admittance matrix at '+format(frequencies[0], '.1f')+' Hz')
    plt.savefig(results_folder + '\\' + file_root + "_Node.pdf", format="pdf", bbox_inches="tight")
    with open(results_folder + '\\' + file_root + "_Node.pickle", 'wb') as f: pickle.dump(plt.gcf(), f)
    plt.close()
    np.seterr(divide='warn')

    # Perform stability analysis
    L = np.matmul(np.linalg.inv(Yedge),Ynode)  # Loop gain matrix
    # Stability via eigenvalue loci
    nyquist(L=L,frequencies=frequencies,results_folder=results_folder,filename=file_root,verbose=True,check_conditioning=check_conditioning,condition_number_th=condition_number_th)
    # Stability via determinant
    nyquist_det(L=L,frequencies=frequencies,results_folder=results_folder,filename=file_root, verbose=False, offset=0.0, draw_arrows=True, show_plot=False)

    # Oscillatory frequencies and participation factors based on the closed-loop impedance matrix (bus)
    EVD(np.linalg.inv(Yedge+Ynode), frequencies, node_variables, results_folder, file_root)

    # Save the admittance matrices
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
    results = [L[:, row, col] for row in range(len(node_variables)) for col in range(len(node_variables))]
    results.insert(0, frequencies)
    # elements = [str(row) + "-" + str(col) for row in range(len(node_variables)) for col in range(len(node_variables))]
    results = tuple(results)
    np.savetxt(results_folder+'\\'+file_root+'_Minor_loop_gain.txt',np.stack(results, axis=1),delimiter='\t',
               header="f\t" + "\t".join(node_variables), comments='')

    # Compute the passivity and SVD of the system matrices
    # passivity(G=Yedge_aux, frequencies=frequencies, results_folder=results_folder, filename=file_root + "_Yedge_aux",variables=edge_aux_variables)
    passivity(G=Ynode,frequencies=frequencies,results_folder=results_folder,filename=file_root+"_Ynode",variables=node_variables)
    passivity(G=Yedge, frequencies=frequencies, results_folder=results_folder, filename=file_root + "_Yedge")
    passivity(G=Ynode+Yedge, frequencies=frequencies, results_folder=results_folder, filename=file_root+"_Ynode_+_Yedge")
    small_gain(np.linalg.inv(Yedge), Ynode, frequencies, results_folder, file_root, variables=node_variables)

def passivity(G, frequencies, results_folder=None, filename=None, variables=None, Yedge=None):
    # The passivity index is computed as half of the minimum eigenvalue of the matrix plus its conjugate transpose
    # min{eig(A + A')}/2
    # A passive system has its Nyquist plot in the RHP: increased stability margins if connected to a passive system
    passivity_index = np.real(np.min(np.linalg.eig(G + G.swapaxes(-1, -2).conj())[0], axis=1))/2  # min{eig(A + A')}/2
    # The eigenvalues of a Hermitian matrix are always real but floating-point arithmetic renders a small complex part
    # and thus we take the real value at the end to get rid of the spurious numerical artifact.

    if (results_folder is not None) and (filename is not None):
        if not path.exists(results_folder): makedirs(results_folder)  # Create results folder if it does not exist
        
        # Plot the passivity index over the frequency range
        if variables is None:
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
            ax.plot(frequencies, passivity_index, color='blue', linestyle='solid', linewidth=2.0,label=r"$\mathbf{Y}_{node}$")
            ax.set_xscale("log")
            ax.minorticks_on()
            ax.grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
            ax.grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
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
        else:
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
            ax[0].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
            ax[0].set_ylabel(r'min $\{ \lambda (\mathbf{G}_i + \mathbf{G}_i^H) \}/2$')
            ax[0].set_title('Passivity analysis for ' + str(len(frequencies)) + ' frequencies')
            ax[0].set_xlim([frequencies[0], frequencies[-1]])
            ax[0].legend(loc='upper left', fancybox=True, shadow=True, ncol=1)

            ax[1].plot(frequencies,passivity_index,color='blue',linestyle='solid',linewidth=2.0,label="Complete matrix")
            ax[1].set_xscale("log")
            ax[1].minorticks_on()
            ax[1].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
            ax[1].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
            ax[0].set_ylabel(r'min $\{ \lambda (\mathbf{G} + \mathbf{G}^H) \}/2$')
            ax[1].legend(loc='upper left', fancybox=True, shadow=True, ncol=1)
            ax[1].set_xlim([frequencies[0], frequencies[-1]])
            ax[1].set_xlabel('Frequency [Hz]')
        
        fig.savefig(results_folder+'\\'+filename + "_passivity.pdf", format="pdf", bbox_inches="tight")
        with open(results_folder + '\\' + filename + "_passivity.pickle", 'wb') as f:
            pickle.dump(fig, f)
        np.savetxt(results_folder+'\\'+filename+'_passivity.txt', np.stack((frequencies, passivity_index), axis=1), delimiter='\t',
                   header="f\t" + "Passivity_index", comments='')
        plt.close(fig)

    return passivity_index

def nyquist(L, frequencies, results_folder=None, filename=None, verbose=True, check_conditioning=False, condition_number_th=0.01/5e-9):
    # Stability assessment via the Generalized Nyquist Criteria (GNC) over the open-loop frequency response

    if not path.exists(results_folder): makedirs(results_folder)  # Create results folder if it does not exist

    # 1) Compute the eigenvalues of the loop-gain at every frequency
    eigenvalues = np.linalg.eig(L)[0]
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

    # 2) Eigenloci plot and count clockwise and counter-clockwise encirclements of (-1,0j) for each eigenvalue
    # Compute the coordinates of the eigenloci for the GNC aplication and plotting
    x = np.real(eigenvalues_sorted)  # Real axis
    y = np.imag(eigenvalues_sorted)  # Imaginary axis
    cw = []  # List of clockwise encirclements for each locus
    ccw = []  # List of counter-clockwise encirclements for each locus
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(6, 7))  # Create the figure and get the colors cycle
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for col in plt.rcParams['axes.prop_cycle'].by_key()['color']: colors.append(col)  # Triplicate the color cycle
    for col in plt.rcParams['axes.prop_cycle'].by_key()['color']: colors.append(col)  # in case of many eigenvalues
    # Loop over the sorted eigenvalues for ploting the locus and count the encirclements
    for idx in range(eigenvalues_sorted.shape[1]):
        # Plot the eigenvalue locus
        ax[0].plot(x[:,idx], y[:,idx], color=colors[idx], linestyle='solid', linewidth=2.0, label=r'$\lambda_{'+format(idx+1,'.0f')+r'}$')
        ax[0].plot(x[:,idx], -y[:,idx], color=colors[idx], linestyle='solid', linewidth=2.0, label='_nolegend_')
        ax[1].plot(x[:,idx],y[:,idx], color=colors[idx],linestyle='solid',linewidth=2.0,label=r'$\lambda_{' + format(idx+1,'.0f')+r'}$')
        ax[1].plot(x[:,idx], -y[:,idx], color=colors[idx], linestyle='solid', linewidth=2.0, label='_nolegend_')

        # Count the number of (-1, 0j) encirclements by this eigenvalue locus
        cwi = 0  # Initialize the counters
        ccwi = 0
        for j in range(1,eigenvalues_sorted.shape[0]):
            # Only consider clockwise crossings of the imaginary axis beyond (-1,0j)
            if y[j - 1, idx] < 0 < y[j, idx] and (x[j,idx] < -1):  # x[j-1,idx] < -1 or
                # Check that the (-1,0j) is to the right of the line between (x1,y1) and (x2,y2)
                # If the cross product of vectors (x2-x1, y2-y1) and (-1-x1, 0-y1) is < 0, then (-1,0) is to the right
                if (x[j,idx] - x[j-1,idx])*(0 - y[j-1,idx]) - (y[j,idx] - y[j-1,idx])*(-1 - x[j-1,idx]) < 0:
                    cwi += 1
                    if verbose: print("CW crossing around ",round(0.5*(frequencies[j] + frequencies[j-1]),4)," Hz by lambda =",str(idx+1))
                    fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(6, 7))
                    ax1.plot([x[j-1,idx],x[j,idx]],[y[j-1,idx],y[j,idx]], color='red', linestyle='solid', linewidth=2.0, label='_nolegend_')
                    ax1.scatter(x[j - 1, idx], y[j - 1, idx], color='green', label=str(frequencies[j-1]))
                    ax1.scatter(x[j, idx], y[j, idx], color='blue', label=str(frequencies[j]))
                    ax1.scatter(-1, 0, marker="+", c='black', label=r'$( -1, 0j )$')
                    ax1.legend(loc='upper right', ncol=1)
                    ax1.minorticks_on()
                    ax1.grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
                    ax1.grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
                    # fig1.show()  # Visualize the plot interactively
                    with open(results_folder + '\\' + filename + "_GNC_lambda_"+str(idx)+"_cw_"+str(cwi)+".pickle", 'wb') as f: pickle.dump(fig1, f)
                    plt.close(fig1)
                    # fig1.clf()

            # Only consider counter-clockwise crossings of the imaginary axis beyond (-1,0j)
            elif y[j - 1, idx] > 0 > y[j, idx] and (x[j-1,idx] < -1 or x[j,idx] < -1):
                if (x[j,idx] - x[j-1,idx])*(0 - y[j-1,idx]) - (y[j,idx] - y[j-1,idx])*(-1 - x[j-1,idx]) > 0:
                    # If the cross product of vectors (x2-x1, y2-y1) and (-1-x1, 0-y1) is > 0, then (-1,0) is to the left
                    # if (x[j,idx] - x[j-1,idx])*(0 - y[j-1,idx]) - (y[j,idx] - y[j-1,idx])*(-1 - x[j-1,idx]) > 0:
                    ccwi += 1
                    if verbose: print("CCW crossing around ",round(0.5*(frequencies[j] + frequencies[j-1]),4)," Hz by lambda =",str(idx+1))
                    fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(6, 7))
                    ax1.plot([x[j-1,idx],x[j,idx]],[y[j-1,idx],y[j,idx]], color='red', linestyle='solid', linewidth=2.0, label='_nolegend_')
                    ax1.scatter(x[j - 1, idx], y[j - 1, idx], color='green', label=str(frequencies[j-1]))
                    ax1.scatter(x[j, idx], y[j, idx], color='blue', label=str(frequencies[j]))
                    ax1.scatter(-1, 0, marker="+", c='black', label=r'$( -1, 0j )$')
                    ax1.legend(loc='upper right', ncol=1)
                    ax1.minorticks_on()
                    ax1.grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
                    ax1.grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
                    # fig1.show()  # Visualize the plot interactively
                    with open(results_folder + '\\' + filename + "_GNC_lambda_"+str(idx)+"_ccw_"+str(ccwi)+".pickle", 'wb') as f: pickle.dump(fig1, f)
                    plt.close(fig1)

        cw.append(cwi)  # Add the counters to the list
        ccw.append(ccwi)
    # print("CC: ",cw,"\nCCW: ",ccw)
    N = sum(cw) - sum(ccw)  # Net number of clockwise encirclements
    if N > 0:
        stable_system = False
        if verbose: print("\n GNC stability assessment: unstable closed-loop system \n")
    elif N < 0:
        stable_system = False
        if verbose: print("\n GNC stability assessment: unstable subsystem \n")
    else:
        stable_system = True
        if verbose: print("\n GNC stability assessment: stable closed-loop system if subsystems are stable \n")

    # Plot the unit circle and the critical point
    th = np.linspace(-np.pi * 1.01, np.pi * 1.01, 314)
    ax[0].plot(np.cos(th), np.sin(th), color='black', linestyle='dotted', linewidth=2.0, label='Unit circle')
    ax[0].scatter(-1, 0, marker="+", c='blue', label=r'$( -1, 0j )$')
    ax[0].minorticks_on()
    ax[0].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
    ax[0].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
    ax[0].set_title('Eigenloci between '+format(frequencies[0],'.1f')+' and '+format(frequencies[-1],'.1f') + ' Hz')
    ax[0].set_xlim([np.min(x,axis=None), np.max(x,axis=None)])
    ax[0].set_ylim([np.min(np.concatenate((-y,y)),axis=None), np.max(np.concatenate((-y,y)),axis=None)])
    ax[0].set_xlabel('Real axis')
    ax[0].set_ylabel('Imaginary axis')
    ax[0].legend(loc='upper right', ncol=4)

    ax[1].plot(np.cos(th), np.sin(th), color='black', linestyle='dotted', linewidth=2.0, label='Unit circle')
    ax[1].scatter(-1, 0,s=4*rcParams['lines.markersize'] ** 2, marker="+", c='blue', label=r'$( -1, 0j )$')
    ax[1].set_xlim([-2.0, 2.0])
    ax[1].set_ylim([-2.0, 2.0])
    # ax[1].minorticks_on()
    ax[1].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
    ax[1].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
    # plt.show()  # Visualize the plot interactively
    fig.savefig(results_folder + '\\' + filename + "_GNC.pdf", format="pdf", bbox_inches="tight")
    with open(results_folder + '\\' + filename + "_GNC.pickle", 'wb') as f:
        pickle.dump(fig, f)
    plt.close(fig)

    # Save the eigenloci
    loci = [eigenvalues_sorted[:, idx] for idx in range(eigenvalues_sorted.shape[1])]
    loci.insert(0, frequencies)
    loci = tuple(loci)
    np.savetxt(results_folder + '\\' + filename + '_GNC.txt', np.stack(loci, axis=1), delimiter='\t',
               header="f\t"+"\t".join(["lambda_"+format(idx+1,'.0f') for idx in range(eigenvalues_sorted.shape[1])]), comments='')

    return stable_system

def small_gain(G1, G2, frequencies, results_folder=None, filename=None, variables=None):
    # Applies a conservative version of the small-gain theorem as |L| = |G1*G2| <= |G1|*|G2| < 1

    S1 = np.linalg.svd(G1, compute_uv=False)
    S2 = np.linalg.svd(G2, compute_uv=False)

    if not path.exists(results_folder): makedirs(results_folder)  # Create results folder if it does not exist

    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
    ax[0].plot(frequencies, 1.0/np.max(S1, axis=1), color='blue', linestyle='solid', linewidth=2.0, label=r"1 / max $\sigma (\mathbf{G}_1)$")
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
    ax[0].plot(frequencies, np.max(S2, axis=1), color='red', linestyle='dashed', linewidth=2.0,label=r"max $\sigma (\mathbf{G}_2)$")

    # Setings for upper plot
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].minorticks_on()
    ax[0].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
    ax[0].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
    ax[0].set_ylabel(r'max $\sigma ( \cdot )$')
    ax[0].set_title('Singular value analysis over ' + str(len(frequencies)) + ' frequencies')
    ax[0].set_xlim([frequencies[0], frequencies[-1]])
    # ax[0].set_xlabel('Frequency [Hz]')
    ax[0].legend(loc='best', fancybox=True, shadow=True, ncol=1)

    S12 = np.linalg.svd(np.matmul(G1,G2), compute_uv=False)
    ax[1].plot(frequencies, np.max(S12, axis=1), color='black', linestyle='solid', linewidth=2.0, label=r"max $\sigma (\mathbf{G}_1  \mathbf{G}_2)$")
    ax[1].plot(frequencies, np.multiply(np.max(S1, axis=1), np.max(S2, axis=1)), color='green', linestyle='dashed',
            linewidth=2.0, label=r"max $\sigma (\mathbf{G}_1) \cdot $ max $\sigma (\mathbf{G}_2)$")
    ax[1].plot([frequencies[0], frequencies[-1]],[1, 1], color='grey', linestyle='dotted',linewidth=2.0, label='_nolegend_')
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[1].minorticks_on()
    ax[1].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
    ax[1].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
    ax[1].set_ylabel('Unitless')
    ax[1].set_xlim([frequencies[0], frequencies[-1]])
    ax[1].set_xlabel('Frequency [Hz]')
    ax[1].legend(loc='best', fancybox=True, shadow=True, ncol=1)
    fig.savefig(results_folder + '\\' + filename + "_gain.pdf", format="pdf", bbox_inches="tight")

    with open(results_folder + '\\' + filename + "_gain.pickle", 'wb') as f:
        pickle.dump(fig, f)
    np.savetxt(results_folder + '\\' + filename + '_gain.txt',
                np.stack((frequencies, np.max(S1, axis=1), np.max(S2, axis=1), np.max(S12, axis=1)), axis=1),
                delimiter='\t', header="f\t" + "max_sigma_G1\t" + "max_sigma_G2\t" + "max_sigma_G1_G2", comments='')
    plt.close(fig)

def EVD(G, frequencies, bus_names=None, results_folder=None, filename=None, verbose=True):
    if bus_names is None: bus_names = [str(bus+1) for bus in range(G.shape[1])]  # Sorted numbers if names not provided
    if not path.exists(results_folder): makedirs(results_folder)  # Create results folder if it does not exist

    # 1) Eigenvalue decomposition over the frequency
    eigenvalues, right_eigenvectors = np.linalg.eig(G)

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
    lambda_imag = np.real(eigenvalues_sorted)  # Imaginary part
    lambda_abs = np.abs(eigenvalues_sorted)  # Absolute value (magnitude)

    # 2.1) Oscillation mode identification based on the magnitude peaks of G = closed-loop impedance matrix
    idx_lambda_max = np.argmax(lambda_abs,axis=0)  # Frequency index of the maximum magnitude of each eigenvalue
    idx_lambda_max_max = np.argmax([lambda_abs[idx_lambda_max[idx],idx] for idx in range(eigenvalues.shape[1])])  # Critical mode = the highest mag peak
    freq_idx = idx_lambda_max[idx_lambda_max_max]  # Oscillation frequency index; or also freq_indices[idx_lambda_min]
    if verbose: print("The main oscillation frequency is",round(frequencies[freq_idx],2),"Hz based on the magnitude of eigenvalue",idx_lambda_max_max+1,"=",round(eigenvalues_sorted[idx_lambda_max[idx_lambda_max_max],idx_lambda_max_max], 5))

    # # 2.2) Based on the minimum real part at imaginary zero-crossing (used in the Positive Net Damping criterion)
    # sign_changes = np.diff(np.sign(lambda_imag),axis=0)
    # critical_points = np.nonzero(sign_changes)  # Frequency of sign change from + to - or viceversa of the imaginary parts
    # freq_indices = critical_points[0]
    # lambda_indices = critical_points[1]
    # # Frequency index of the minimum real part of each eigenvalue over the critical frequencies
    # critical_lambdas = [lambda_re[freq_indices[idx], lambda_indices[idx]] for idx in range(len(freq_indices))]  # Real part at the critical frequencies
    # idx_lambda_min = np.argmin(critical_lambdas)
    # print("According to the PND, the critical mode is at",round(frequencies[freq_indices[idx_lambda_min]],3),"Hz based on the minimum real part of eigenvalue",lambda_indices[idx_lambda_min]+1)
    # print(" Eigenvalue",lambda_indices[idx_lambda_min]+1, "=",round(eigenvalues_sorted[freq_indices[idx_lambda_min],lambda_indices[idx_lambda_min]],5))
    # print(" The slope of the imaginary part is",round(lambda_imag[freq_indices[idx_lambda_min],lambda_indices[idx_lambda_min]]-lambda_imag[freq_indices[idx_lambda_min]-1,lambda_indices[idx_lambda_min]],5))
    # freq_idx = freq_indices[idx_lambda_min] # Oscillation frequency index or also round(critical_lambdas[idx_lambda_min],5)

    # 3) Compute the bus participation factors (PFs) of the critical eigenvalue at the oscillation frequency
    # Controllability (right eigenvectors) and observability (transpose of left eigenvectors)
    Obs = right_eigenvectors[freq_idx, :]
    Cont = np.transpose(left_eigenvectors[freq_idx, :])
    # PF[row = bus, column = mode]
    PF = right_eigenvectors[freq_idx, :] * np.transpose(left_eigenvectors[freq_idx, :])  # Element-wise product
    PF_mode = PF[:,idx_lambda_max_max]  # Select the target mode
    # The controllability, observability and PF of the critical mode at each bus
    if verbose: print("Bus"+(max([len(bus) for bus in bus_names])-3)*" "+"\t","Cont.\t","Obs.\t","PF")  # Header
    for idx, bus in enumerate(bus_names):
        if verbose: 
            print(bus+"\t", f"{np.abs(Cont[idx,idx_lambda_max_max]) / np.sum(np.abs(Cont[:,idx_lambda_max_max])):.4f}"+"\t",
              f"{np.abs(Obs[idx,idx_lambda_max_max]) / np.sum(np.abs(Obs[:,idx_lambda_max_max])):.4f}"+"\t",
              f"{np.abs(PF_mode[idx]) / np.sum(np.abs(PF_mode)):.4f}")

    # 4) Plot the eigenvalues over frequency
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(6, 8))  # Create the figure and get the colors cycle
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for col in plt.rcParams['axes.prop_cycle'].by_key()['color']: colors.append(col)  # Triplicate colour cycle
    for col in plt.rcParams['axes.prop_cycle'].by_key()['color']: colors.append(col)  # in case of many eigenvalues
    # Loop over the sorted eigenvalues and plot them over the frequency range
    for idx in range(eigenvalues_sorted.shape[1]):
        # Plot the eigenvalue locus
        ax[0].plot(frequencies, lambda_abs[:,idx], color=colors[idx], linestyle='solid', linewidth=2.0,
                   label=r'$\lambda_{' + format(idx + 1, '.0f') + r'}$')
        ax[1].plot(frequencies, lambda_re[:, idx], color=colors[idx], linestyle='solid', linewidth=2.0,
                   label=r'$\lambda_{' + format(idx + 1, '.0f') + r'}$')
        ax[2].plot(frequencies, lambda_imag[:, idx], color=colors[idx], linestyle='solid', linewidth=2.0,
                   label=r'$\lambda_{' + format(idx + 1, '.0f') + r'}$')

    # Figure settings and save to pdf
    ax[0].minorticks_on()
    ax[0].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
    ax[0].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
    ax[0].set_title('Eigenvalue decomposition between ' + format(frequencies[0], '.1f') + ' and ' + format(frequencies[-1], '.1f') + ' Hz')
    ax[0].set_xlim([frequencies[0], frequencies[-1]])
    ax[0].set_ylim([np.min(lambda_abs, axis=None), np.max(lambda_abs, axis=None)])
    ax[0].set_ylabel('Magnitude')
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].legend(loc='lower right', ncol=4)

    ax[1].minorticks_on()
    ax[1].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
    ax[1].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
    ax[1].set_xlim([frequencies[0], frequencies[-1]])
    ax[1].set_ylim([np.min(lambda_re, axis=None), np.max(lambda_re, axis=None)])
    ax[1].set_ylabel('Real part')
    ax[1].set_xscale("log")

    ax[2].minorticks_on()
    ax[2].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
    ax[2].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
    ax[2].set_xlim([frequencies[0], frequencies[-1]])
    ax[2].set_ylim([np.min(lambda_imag, axis=None), np.max(lambda_imag, axis=None)])
    ax[2].set_xlabel('Frequency [Hz]')
    ax[2].set_ylabel('Imaginary part')
    ax[2].set_xscale("log")

    # plt.show()  # Visualize the plot interactively
    fig.savefig(results_folder + '\\' + filename + "_EVD.pdf", format="pdf", bbox_inches="tight")
    with open(results_folder + '\\' + filename + "_EVD.pickle", 'wb') as f:
        pickle.dump(plt.gcf(), f)
    plt.close(fig)

    # Save the EVD into a text file
    evd_results = [eigenvalues_sorted[:,idx] for idx in range(eigenvalues_sorted.shape[1])]
    evd_results.insert(0, frequencies)
    evd_results = tuple(evd_results)
    np.savetxt(results_folder + '\\' + filename + '_EVD.txt', np.stack(evd_results, axis=1), delimiter='\t',
               header="Frequency [Hz]\t" + "\t".join(bus_names), comments='')

def nyquist_det(L, frequencies, results_folder=None, filename=None, verbose=True, offset=0.0, draw_arrows=True, show_plot=False):
    # Stability assessment based on the determinant of I + L

    if not path.exists(results_folder): makedirs(results_folder)  # Create results folder if it does not exist
    
    det = np.linalg.det(np.identity(L.shape[1]) + L) + offset  
    x = np.real(det)
    y = np.imag(det)

    # Count the number of (offset, 0j) encirclements by det
    cwi = 0  # Initialize the counters
    ccwi = 0
    for j in range(1,len(frequencies)):
        # Only consider clockwise crossings of the imaginary axis beyond (offset,0j)
        if y[j - 1] < 0 < y[j] and (x[j] < offset):  # x[j-1] < offset or
            # Check that the (offset,0j) is to the right of the line between (x1,y1) and (x2,y2)
            # If the cross product of vectors (x2-x1, y2-y1) and (offset-x1, 0-y1) is < 0, then (offset,0) is to the right
            if (x[j] - x[j-1])*(0 - y[j-1]) - (y[j] - y[j-1])*(offset - x[j-1]) < 0:
                cwi += 1
                if verbose: print("CW crossing around ",round(0.5*(frequencies[j] + frequencies[j-1]),4)," Hz")
                fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(6, 7))
                ax1.plot([x[j-1],x[j]],[y[j-1],y[j]], color='red', linestyle='solid', linewidth=2.0, label='_nolegend_')
                ax1.scatter(x[j-1], y[j-1], color='green', label=str(frequencies[j-1]))
                ax1.scatter(x[j], y[j], color='blue', label=str(frequencies[j]))
                ax1.scatter(offset, 0, marker="+", c='black', label=r'$( -1, 0j )$')
                ax1.legend(loc='upper right', ncol=1)
                ax1.minorticks_on()
                ax1.grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
                ax1.grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
                if show_plot: fig1.show()  # Visualize the plot interactively
                with open(results_folder + '\\' + filename + "_det_cw_"+str(cwi)+".pickle", 'wb') as f: pickle.dump(fig1, f)
                plt.close(fig1)

        # Only consider counter-clockwise crossings of the imaginary axis beyond (offset,0j)
        elif y[j - 1] > 0 > y[j] and (x[j-1] < offset or x[j] < offset):
            if (x[j] - x[j-1])*(0 - y[j-1]) - (y[j] - y[j-1])*(offset - x[j-1]) > 0:
                # If the cross product of vectors (x2-x1, y2-y1) and (offset-x1, 0-y1) is > 0, then (offset,0) is to the left
                ccwi += 1
                if verbose: print("CCW crossing around ",round(0.5*(frequencies[j] + frequencies[j-1]),4)," Hz")
                fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(6, 7))
                ax1.plot([x[j-1],x[j]],[y[j-1],y[j]], color='red', linestyle='solid', linewidth=2.0, label='_nolegend_')
                ax1.scatter(x[j-1], y[j-1], color='green', label=str(frequencies[j-1]))
                ax1.scatter(x[j], y[j], color='blue', label=str(frequencies[j]))
                ax1.scatter(offset, 0, marker="+", c='black', label=r'$( -1, 0j )$')
                ax1.legend(loc='upper right', ncol=1)
                ax1.minorticks_on()
                ax1.grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
                ax1.grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
                if show_plot: fig1.show()  # Visualize the plot interactively
                with open(results_folder + '\\' + filename + "_det_ccw_"+str(ccwi)+".pickle", 'wb') as f: pickle.dump(fig1, f)
                plt.close(fig1)

    N = cwi - ccwi  # Net number of clockwise encirclements
    if N > 0:
        stable_system = False
        if verbose: print("\n GNC stability assessment: unstable closed-loop system \n")
    elif N < 0:
        stable_system = False
        if verbose: print("\n GNC stability assessment: unstable subsystem \n")
    else:
        stable_system = True
        if verbose: print("\n GNC stability assessment: stable closed-loop system if subsystems are stable \n")

    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(6, 7))  # Create the figure and get the colors cycle
    # Plot the critical point
    ax[0].scatter(offset, 0, marker="+", c='blue', label=r'$( '+str(round(offset,0))+', 0j )$')
    ax[0].plot(x, y, color='red', linestyle='solid', linewidth=2.0)
    if draw_arrows:
        da = int(np.log(L.shape[0] + 1))  # decimate the number of arrows
        ax[0].quiver(x[::da], y[::da], np.gradient(x)[::da], np.gradient(y)[::da], angles='xy', scale_units='xy',
                    scale=da, linewidth=1, edgecolor='black', facecolor='green')
    id0 = np.argmin(np.abs(frequencies - 50.0))
    ax[0].text(x[id0], y[id0], str(round(frequencies[id0], 2)) + ' Hz', fontsize=12, color='blue', ha='right', va='bottom')
    ax[0].text(x[-1], y[-1], str(round(frequencies[-1], 2)) + ' Hz', fontsize=12, color='blue', ha='right', va='bottom')
    ax[0].text(x[0], y[0], str(round(frequencies[0], 2)) + ' Hz', fontsize=12, color='blue', ha='right', va='bottom')

    ax[0].minorticks_on()
    ax[0].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
    ax[0].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
    ax[0].set_title(str(offset)+r' + det[I + L(s)] between '+format(frequencies[0], '.1f')+' and ' + format(frequencies[-1], '.1f') + ' Hz')
    ax[0].set_xlim([np.min(x, axis=None), np.max(x, axis=None)])
    ax[0].set_ylim([np.min(y, axis=None), np.max(y, axis=None)])
    ax[0].set_xlabel('Real axis')
    ax[0].set_ylabel('Imaginary axis')
    ax[0].legend(loc='best', ncol=1)

    ax[1].plot(x, y, color='red', linestyle='solid', linewidth=2.0)
    if draw_arrows:
        ax[1].quiver(x, y, np.gradient(x), np.gradient(y), angles='xy', scale_units='xy',
                    scale=1, linewidth=1, edgecolor='black', facecolor='green')
    ax[1].scatter(offset, 0, s=4 * rcParams['lines.markersize'] ** 2, marker="+", c='blue', label=r'$( ' + str(offset) + ', 0j )$')
    ax[1].set_xlim([-1.0+offset, 1.0+offset])
    ax[1].set_ylim([-1.0+offset, 1.0+offset])
    ax[1].minorticks_on()
    ax[1].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
    ax[1].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
    ax[1].set_xlabel('Real axis')
    ax[1].set_ylabel('Imaginary axis')

    fig.savefig(results_folder + '\\' + filename + "_det.pdf", format="pdf", bbox_inches="tight")
    with open(results_folder + '\\' + filename + "_det.pickle", 'wb') as f:
        pickle.dump(fig, f)
    if show_plot: plt.show()  # Visualize the plot interactively
    plt.close(fig)

    # Save the results
    np.savetxt(results_folder + '\\' + filename + '_det.txt', np.stack((frequencies,det), axis=-1), delimiter='\t',
               header="Frequency [Hz]\t"+str(offset)+"+det[I + L(s)]", comments='')
    
    return stable_system


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
The EVD of L is saved as filename_GNC.txt and its plot is saved as filename_GNC.pdf.

Required arguments
        L                   (numpy ndarray of complex double) Minor loop gain (transfer matrix) for different frequencies.
        frequencies         (numpy array) Frequencies over which L is computed [Hz].
        results_folder      Absolute path where the results are to be stored. If it does not exist, it is created.
        filename            Name root of the results output files.        
      
Optional arguments
        verbose             Bool flag to show detailed GNC application information, such as the number of counter clock-wise (CCW) and clock-wise (CC) encirclements of the critical point.
        check_conditioning  Bool flag to discard values with poor numerical conditioning of L.
        condition_number_th (double) Condition number threshold of L above which the data is ignored.
                            This threshold can be set based on the expected input error and maximum acceptable ourput error.
                            For example, for a relative output error <= 0.01 considering a relative input error of 5e-9, the condition number threshold can be set to 0.01/5e-9 (default value).

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

Required arguments
        L                   (numpy ndarray of complex double) Minor loop gain (transfer matrix) for different frequencies.
        frequencies         (numpy array) Frequencies over which L is computed [Hz].
        results_folder      Absolute path where the results are to be stored. If it does not exist, it is created.
        filename            Name root of the results output files.        
      
Optional arguments
        verbose             Bool flag to show detailed GNC application information, such as the number of counter clock-wise (CCW) and clock-wise (CC) encirclements of the critical point. Default = True.
        draw_arrows         Bool flag to draw arrows on the direction of the Nyquist plot. Default = True.
        offset              (double) The offset parameter can be used to shift the critical point on the real axis instead of (0,j0). Default = 0.
                            For example, offset = -1.0 defines the critical point as (-1,j0) as when applying the GNC via the eigenvalue loci of L.
        show_plot           Bool to show the plot interactively. Default = False.

Returns
        Bool flag indicating closed-loop or interconnected stability: True means stable.

"""