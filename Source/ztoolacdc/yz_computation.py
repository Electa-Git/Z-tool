"""
This program contains a function to compute multi-infeed admittance matrices both for AC and DC systems

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

__all__ = ['admittance','SISO_TF']

import matplotlib.pyplot as plt
import numpy as np  # Numerical python functions
from scipy import interpolate  # For interpolation
from matplotlib import rcParams  # Text's parameters for plots
from .plot_utils import bode_plot

rcParams['mathtext.fontset'] = 'cm'  # Font selection
rcParams['font.family'] = 'STIXGeneral'  # 'cmu serif'

def admittance(f_base=None, frequencies=None, freq_multi=None, fft_periods=1, 
               sides=None, dt=None, exploit_dq_sym=False, start_idx=None, 
               zblocks=None, results_folder=None, results_name='Y', 
               multiport=None, make_plot=True):
    
    # Normalize frequency input
    if freq_multi is not None:
        if frequencies is not None:
            raise ValueError("Only one of 'frequencies' or 'freq_multi' should be provided")
        is_multi_freq = True
        freq_multi_shape = freq_multi.shape
        frequencies_flat = freq_multi.reshape(np.prod(freq_multi.shape))
        num_sim = freq_multi_shape[1]
    else:
        if frequencies is None:
            raise ValueError("Either 'frequencies' or 'freq_multi' must be provided")
        is_multi_freq = False
        freq_multi_shape = None
        frequencies_flat = np.asarray(frequencies)
        num_sim = 1
    
    # Normalize zblocks to list
    if not isinstance(zblocks, list): zblocks = [zblocks]
    
    # Normalize sides to list
    if isinstance(sides, str):
        sides = [sides] * len(zblocks)
    elif not isinstance(sides, list):
        sides = list(sides)
    
    # Verify multiport requirement
    if multiport is None and len(zblocks) > 2:
        raise ValueError("multiport object is required for more than 2 blocks")
    
    # Determine scan type
    if multiport is not None:
        scan_type = multiport.scan_type
    else:
        # Infer from zblocks types
        if len(zblocks) == 1:
            scan_type = zblocks[0].type
        else:
            # Multiple blocks without multiport - infer as mixed
            types = set(block.type for block in zblocks)
            if types == {"AC"}:
                scan_type = "AC"
            elif types == {"DC"}:
                scan_type = "DC"
            elif types == {"AC", "DC"}:
                scan_type = "ACDC" # Mixed AC/DC subsystem
            else:
                raise ValueError(f"Cannot infer scan_type from blocks with types: {types}")
    
    # Determine matrix size NxN
    if multiport is not None:
        if exploit_dq_sym and multiport.scan_type == "AC":
            N = 2 * multiport.runs
        else:
            N = multiport.runs
    else:
        # Simple single-port scan or AC/DC subsystem
        if scan_type == "AC":
            N = 2  # d and q axes
        elif scan_type == "DC":
            N = 1
        elif scan_type == "ACDC":
            N = 3  # dc, d, q
        else:
            N = len(zblocks) * (2 if any(b.type == "AC" for b in zblocks) else 1)
    
    # Define perturbation types to iterate over
    if multiport is not None:
        # Use multiport's run list
        data_ending = ["_" + str(run) for run in multiport.runs_list]
    else:
        # Define based on scan type
        if scan_type == "AC":
            data_ending = ["_d", "_q"]
        elif scan_type == "DC":
            data_ending = ["_dc"]
        elif scan_type == "ACDC":
            data_ending = ["_dc", "_d", "_q"]
        else:
            # Generic multiport
            data_ending = []
            for block in zblocks:
                if block.type == "AC":
                    data_ending.extend(["_d", "_q"])
                else:
                    data_ending.append("_dc")
    
    # Initialization of the matrix and data storage
    num_freqs = len(frequencies_flat)
    L = int(fft_periods * 1 / f_base * 1.0 / dt) # FFT length
    deltaV = np.empty((num_freqs, N, N), dtype='cdouble')
    deltaI = np.empty((num_freqs, N, N), dtype='cdouble')
    Y = np.empty((num_freqs, N, N), dtype='cdouble')
    
    # Main admittance computation loop
    for freq_idx, frequency in enumerate(frequencies_flat):
        # Determine which simulation index to use
        fft_idx = int(round(frequency * fft_periods * 1 / f_base))
        if is_multi_freq:
            sim_idx = freq_idx % num_sim
        else:
            sim_idx = freq_idx
        
        # Fill in the perturbation data for this frequency row by row for each column
        for col, sim_type in enumerate(data_ending):
            current_row = 0
            for block_num, block in enumerate(zblocks):
                # AC scan block
                if block.type == "AC":
                    names = ['VDUTac:1', 'VDUTac:2', 'IDUTacA' + sides[block_num] + ':1', 'IDUTacA' + sides[block_num] + ':2']
                    
                    if multiport is not None and exploit_dq_sym and col < N // 2:
                        # Symmetry exploitation: only d-axis perturbation data exists and q-axis scan is derived via symmetry
                        for name_pos, name in enumerate(names):
                            # Retrieve and compute FFT of perturbation
                            delta = (block.perturbation_data[sim_idx][name + sim_type][start_idx:] - block.snapshot_data[name][start_idx:])
                            delta_FD = np.fft.rfft(delta, n=L, axis=0) * 2 / L
                            row = current_row + (name_pos % 2)  # Alternating d,q rows
                            if "V" in name:
                                deltaV[freq_idx, row, 2*col] = delta_FD[fft_idx]
                            else:
                                deltaI[freq_idx, row, 2*col] = delta_FD[fft_idx]
                        
                        # Apply dq symmetry to derive q-axis results from d-axis perturbations
                        deltaV[freq_idx, current_row,     2*col+1] = -deltaV[freq_idx, current_row+1, 2*col]
                        deltaV[freq_idx, current_row+1,   2*col+1] = deltaV[freq_idx, current_row, 2*col]
                        deltaI[freq_idx, current_row,     2*col+1] = -deltaI[freq_idx, current_row+1, 2*col]
                        deltaI[freq_idx, current_row+1,   2*col+1] = deltaI[freq_idx, current_row, 2*col]
                    
                    else:
                        # No symmetry exploitation: use actual d and q perturbation data
                        for name_pos, name in enumerate(names):
                            delta = (block.perturbation_data[sim_idx][name + sim_type][start_idx:] -  block.snapshot_data[name][start_idx:])
                            delta_FD = np.fft.rfft(delta, n=L, axis=0) * 2 / L
                            row = current_row + (name_pos % 2)
                            
                            if "V" in name:
                                deltaV[freq_idx, row, col] = delta_FD[fft_idx]
                            else:
                                deltaI[freq_idx, row, col] = delta_FD[fft_idx]
                    current_row += 2
                
                # DC scan block
                elif block.type == "DC":
                    names = ['VDUTdc', 'IDUTdcA' + sides[block_num]]
                    for name in names:
                        delta = (block.perturbation_data[sim_idx][name + sim_type][start_idx:] - 
                               block.snapshot_data[name][start_idx:])
                        delta_FD = np.fft.rfft(delta, n=L, axis=0) * 2 / L
                        
                        if "V" in name:
                            deltaV[freq_idx, current_row, col] = delta_FD[fft_idx]
                        else:
                            deltaI[freq_idx, current_row, col] = delta_FD[fft_idx]
                    
                    current_row += 1
        
        # Admittance computation
        if multiport is not None and multiport.enforce:
            # Enforce network connectivity constraints
            Yextended = np.copy(multiport.adj_matrix)
            np.fill_diagonal(Yextended, 1)
            if scan_type == "AC": Yextended = np.kron(Yextended, np.ones((2, 2), dtype=int))
            Y[freq_idx, :, :] = np.multiply(Yextended, np.matmul(deltaI[freq_idx, :, :], np.linalg.inv(deltaV[freq_idx, :, :])))
        else:
            Y[freq_idx, :, :] = np.matmul(deltaI[freq_idx, :, :],  np.linalg.inv(deltaV[freq_idx, :, :]))
    
    # Sort frequencies and admittance data
    sorting_idx = np.argsort(frequencies_flat)
    frequencies_sorted = frequencies_flat[sorting_idx]
    Y = Y[sorting_idx]
    
    # Save the results
    if results_folder is not None:
        filename = (results_name + '#Y_' + scan_type + '#' + '#'.join([zblocks[idx].name + "-" + sides[idx] for idx in range(len(sides))]))
        results = [Y[:, row, col] for row in range(N) for col in range(N)]
        results.insert(0, frequencies_sorted)

        header = ["f"]
        for block_num, block in enumerate(zblocks):
            if block.type == "AC":
                header.append(block.name + "-" + sides[block_num] + "_d")
                header.append(block.name + "-" + sides[block_num] + "_q")
            elif block.type == "DC":
                header.append(block.name + "-" + sides[block_num] + "_dc")

        np.savetxt(results_folder + '\\' + filename + '#.txt',np.stack(results, axis=1),delimiter='\t', header='\t'.join(header), comments='')
    
    # Plot the admittance
    if results_folder is not None and make_plot:
        if multiport is not None and multiport.enforce:
            Y_plot = Y * Yextended[:, :, np.newaxis].squeeze()  # Apply mask
        else:
            Y_plot = Y
        bode_plot(Y_plot, frequencies_sorted, results_folder, filename)
    
    return Y, frequencies_sorted

def SISO_TF(f_base=None, frequencies=None, fft_periods=1, dt=None, freq_multi=None, 
            start_idx=None, zblock=None, results_folder=None, results_name='SISO_TF', make_plot=True):
    
    # Normalize frequency input
    if freq_multi is not None:
        if frequencies is not None:
            raise ValueError("Only one of 'frequencies' or 'freq_multi' should be provided")
        is_multi_freq = True
        freq_multi_shape = freq_multi.shape
        frequencies_flat = freq_multi.reshape(np.prod(freq_multi.shape),order='F')
        num_sim = freq_multi_shape[1]
        freqs_per_sim = freq_multi_shape[0]

    else:
        if frequencies is None:
            raise ValueError("Either 'frequencies' or 'freq_multi' must be provided")
        is_multi_freq = False
        freq_multi_shape = None
        frequencies_flat = np.asarray(frequencies)
        num_sim = 1
        freqs_per_sim = 1

    # Small-signal sinusoidal steady state computation and rFFT (no target frequency-based FFT distinction)
    L = int(fft_periods * 1 / f_base * 1.0 / dt)  # For the FFT computation

    # For each frequency a transfer function is computed
    num_freqs = len(frequencies_flat)
    delta_input = np.empty(num_freqs, dtype='cdouble')  # Also dtype='csingle'
    delta_output= np.empty(num_freqs, dtype='cdouble')

    # Build the small-signal input and output variables, compute the FFT and retreive the OUT/IN at every frequency
    for name in ['inputTF', 'outputTF']:
        for sim_idx in range(num_freqs if not is_multi_freq else num_sim):
            delta = zblock.perturbation_data[sim_idx][name+"_TF"][start_idx:] - zblock.snapshot_data[name][start_idx:]
            delta_FD = np.fft.rfft(delta, n=L, axis=0) * 2 / L

            for freq_idx, frequency in enumerate(freq_multi[:,sim_idx].squeeze()) if is_multi_freq else enumerate([frequencies_flat[sim_idx]]):
                fft_idx = int(round(frequency * fft_periods * 1 / f_base))
                # Retrieve the response at the target frequency
                if "input" in name:
                    delta_input[freq_idx + sim_idx*freqs_per_sim] = delta_FD[fft_idx]
                else:
                    delta_output[freq_idx + sim_idx*freqs_per_sim] = delta_FD[fft_idx]
            
    Y = delta_output / delta_input

    if results_folder is not None:
        filename = results_name+'#SISO_TF_'+"#"+zblock.name
        sorting_idx = np.argsort(frequencies_flat)
        frequencies = frequencies_flat[sorting_idx]
        Y = Y[sorting_idx]
        results = [Y]
        results.insert(0, frequencies)
        results = tuple(results)
        header = "f\t" + zblock.name
        np.savetxt(results_folder+r'\\'+filename+'#.txt',np.stack(results, axis=1),  delimiter='\t', header=header,comments='')

    if results_folder is not None and make_plot:
        bode_plot(Y, frequencies, results_folder, filename, legend=[zblock.name])
    
    return Y, frequencies

admittance.__doc__ = """
Compute the admittance at every frequency based on the perturbation data obtained with the "frequency_sweep" function.

The function iterates over the time-domain simulation results for every frequency stored in "zblocks" and computes the real-side FFT of the perturbed
waveforms under sinusoidal steady-state. The small-signal voltages and currents at every frequency matrices are built and the admittance is simply
computed as Currents * inv(Voltages). Additional options allow the user to compute the dq-frame admittance with half of the perturbations needed
for a general NxN matrix at each frequency by assuming dq-frame symmetry. 
The function accepts several arguments:

Required
    f_base          Base frequency (determines frequency resolution) [Hz]
    frequencies     List of frequencies at which the perturbations are performed [Hz]
    scantype        To distinguish between the following subsystems
                        single AC bus: scantype = "AC" (default), e.g. three-phase system bus
                        single DC node: scantype = "DC", e.g. DC system node
                        single AC bus and single DC node: cantype = "ACDC", e.g. AC/DC converter
                        generic AC, DC and/or AC/DC multi-terminal: scantype = "Network"
    zblocks         Instance of the Scanblock class in frequency_sweep.py, or list thereof, storing information from the PSCAD Z-tool library blocks such as EMT simulation data.
    sides           Sides of the PSCAD Z-tool library blocks at which the admittance is to be computed.
                    sides is a single number for single AC bus or single DC node cases, while it is a list otherwise.
                    The list contains the side of each block as apearing in the same order as in zblocks.
    dt              Simulation sampling time used to compute the FFT [s]
    start_idx       Index of the time-domain waveforms corresponding to the sinusoidal steady-state after the injections.
    results_folder  Absolute path where the admittance is to be stored.

Optional
    fft_periods     Number of periods used to compute the FFT. Default = 1.
    results_name    Default = "Y".
    network         Instance of the Network class in frequency_sweep.py containing information related to the type of subnetwork, blocks involved and their sides,
                    topology or connectivity, perturbations performed at each EMT run, block identifiers, etc. 
    exploit_dq_sym  The dq-frame symmetry of AC subsystems is exploited when set to True so as to compute their admittance matrix with half of the number of simulations. Default = False.
    
"""

SISO_TF.__doc__ = """
Compute a single-input single-output transfer function based on the perturbation data obtained with the "frequency_sweep_TF" function.
The function iterates over the time-domain simulation results for every simulation stored in "zblock" and computes the real-side FFT of the perturbed
waveforms under sinusoidal steady-state. The small-signal input and output vectors at every frequency are built and the transfer function is simply
computed as Output / Input. The function also supports multi-frequency simulations.
The function accepts several arguments:
Required
    f_base          Base frequency (determines frequency resolution) [Hz]
    frequencies     List of frequencies at which the perturbations are performed [Hz]
    zblock          Instance of the Scanblock class in frequency_sweep.py, storing information from the PSCAD Z-tool library blocks such as EMT simulation data.
    dt              Simulation sampling time used to compute the FFT [s]
    start_idx       Index of the time-domain waveforms corresponding to the sinusoidal steady-state after the injections.
    results_folder  Absolute path where the transfer function is to be stored.
Optional
    fft_periods     Number of periods used to compute the FFT. Default = 1.
    freq_multi      2D array where each column contains the list of frequencies for each simulation. If provided, 'frequencies' is ignored.
    results_name    Default = "SISO_TF".
    make_plot       If True, Bode plots of the transfer function are generated and saved. Default = True.

"""