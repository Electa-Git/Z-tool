"""
This program contains a function to compute multi-infeed admittance matrices both for AC and DC systems

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

__all__ = ['admittance','admittance_multi_freq','admittance_generic','SISO_TF']

import matplotlib.pyplot as plt
import numpy as np  # Numerical python functions
from matplotlib import rcParams  # Text's parameters for plots

rcParams['mathtext.fontset'] = 'cm'  # Font selection
rcParams['font.family'] = 'STIXGeneral'  # 'cmu serif'

def admittance(f_base=None, frequencies=None, fft_periods=1, scantype="AC", sides=None, dt=None, exploit_dq_sym=False,
               start_idx=None, zblocks=None, results_folder=None, results_name='Y', network=None):
    # Small-signal sinusoidal steady state computation and rFFT (no target frequency-based FFT distinction)
    L = int(fft_periods * 1 / f_base * 1.0 / dt)  # For the FFT computation
    if scantype == "AC":
        # Compute the small-signal sinusoidal steady-state waveforms
        # For each simulation (freq) it contains a matrix: [Vd_d Vd_q; Vq_d Vq_q]
        deltaV = np.empty((len(frequencies), 2, 2), dtype='cdouble')  # Also dtype='csingle'
        deltaI = np.empty((len(frequencies), 2, 2), dtype='cdouble')  # Also dtype='csingle'
        Y = np.empty((len(frequencies), 2, 2), dtype='cdouble')  # Also dtype='csingle'
        names = ['VDUTac:1', 'VDUTac:2', 'IDUTacA' + sides + ':1', 'IDUTacA' + sides + ':2']
        row = {names[0]: 0, names[1]: 1, names[2]: 0, names[3]: 1}
        for sim, frequency in enumerate(frequencies):
            fft_idx = int(round(frequency * fft_periods * 1 / f_base))  # Index of the target FFT frequency
            for col, sim_type in enumerate(["_d", "_q"]):
                if exploit_dq_sym:
                    if col != 1:
                        for name in names:
                            delta = zblocks.perturbation_data[sim][name+sim_type][start_idx:] - zblocks.snapshot_data[name][start_idx:]
                            delta_FD = np.fft.rfft(delta, n=L, axis=0) * 2 / L
                            # Retrieve the response at the target frequency
                            if "V" in name:
                                deltaV[sim, row[name], col] = delta_FD[fft_idx]
                            else:
                                deltaI[sim, row[name], col] = delta_FD[fft_idx]
                        
                        # The data for the q-axis is mirrored w.r.t. the d-axis data
                        deltaV[sim, row[names[0]], 1] = - deltaV[sim, row[names[1]], 0]
                        deltaV[sim, row[names[1]], 1] =   deltaV[sim, row[names[0]], 0]
                        deltaI[sim, row[names[2]], 1] = - deltaI[sim, row[names[3]], 0]
                        deltaI[sim, row[names[3]], 1] =   deltaI[sim, row[names[2]], 0]

                else:       
                    for name in names:
                        delta = zblocks.perturbation_data[sim][name+sim_type][start_idx:] - zblocks.snapshot_data[name][start_idx:]
                        delta_FD = np.fft.rfft(delta, n=L, axis=0) * 2 / L
                        # freq_FD = np.fft.rfftfreq(L, d=dt)  # rFFT frequency points
                        # Retrieve the response at the target frequency
                        if "V" in name:
                            deltaV[sim, row[name], col] = delta_FD[fft_idx]
                        else:
                            deltaI[sim, row[name], col] = delta_FD[fft_idx]
            
            Y[sim,...] = np.matmul(deltaI[sim,...], np.linalg.inv(deltaV[sim,...]))


        if results_folder is not None:
            filename = results_name + '#Y_AC#' + zblocks.name + "-" + sides
            np.savetxt(results_folder + '\\' + filename+'#.txt',
                       np.stack((frequencies, Y[:, 0, 0], Y[:, 0, 1], Y[:, 1, 0], Y[:, 1, 1]), axis=1), delimiter='\t',
                       header="f\t"+zblocks.name+"-"+sides+"_d\t"+zblocks.name+"-"+sides+"_q", comments='')
            # "f\td-d\td-q\tq-d\tq-q"
            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
            ax[0].scatter(frequencies, 20 * np.log10(np.abs(Y[:, 0, 0])), marker='o', facecolors='none', edgecolors='b',
                          linewidths=1.5, label=r'$Y_{dd}$')
            ax[0].scatter(frequencies, 20 * np.log10(np.abs(Y[:, 0, 1])), marker='x', c='r', linewidths=1.5,
                          label=r'$Y_{dq}$')
            ax[0].scatter(frequencies, 20 * np.log10(np.abs(Y[:, 1, 0])), marker='+', c='m', linewidths=1.5,
                          label=r'$Y_{qd}$')
            ax[0].scatter(frequencies, 20 * np.log10(np.abs(Y[:, 1, 1])), marker='.', c='g', linewidths=1.5,
                          label=r'$Y_{qq}$')
            # ax[0].set_yscale("log")
            ax[0].set_xscale("log")
            ax[0].set_xlim([frequencies[0], frequencies[-1]])
            ax[0].minorticks_on()
            ax[0].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
            ax[0].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
            ax[0].set_ylabel('Magnitude [dB]')
            ax[0].set_title('DUT admittance ― ' + str(len(frequencies)) + ' scanned frequencies')
            ax[0].legend(loc='upper right', ncol=2)

            ax[1].scatter(frequencies, np.angle(Y[:, 0, 0], deg=True), marker='o', facecolors='none', edgecolors='b',
                          linewidths=1.5, label=r'$Y_{dd}$')
            ax[1].scatter(frequencies, np.angle(Y[:, 0, 1], deg=True), marker='x', c='r', linewidths=1.5,
                          label=r'$Y_{dq}$')
            ax[1].scatter(frequencies, np.angle(Y[:, 1, 0], deg=True), marker='+', c='m', linewidths=1.5,
                          label=r'$Y_{qd}$')
            ax[1].scatter(frequencies, np.angle(Y[:, 1, 1], deg=True), marker='.', c='g', linewidths=1.5,
                          label=r'$Y_{qq}$')
            ax[1].set_xscale("log")
            ax[1].set_ylim([-200, 200])
            ax[1].set_yticks([-180, -90, 0, 90, 180])
            ax[1].set_xlim([frequencies[0], frequencies[-1]])
            ax[1].minorticks_on()
            ax[1].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
            ax[1].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
            ax[1].set_ylabel('Phase [°]')
            ax[1].set_xlabel('Frequency [Hz]')
            ax[1].legend(loc='upper right', fancybox=True, shadow=True, ncol=2)
            fig.savefig(results_folder + '\\' + filename + ".pdf", format="pdf", bbox_inches="tight")
            plt.close(fig)

    elif scantype == "DC":
        # Use the actual steady-state waveforms (allows to remove periodicity)
        # For each simulation (freq) it contains a matrix: [Vd_d Vd_q; Vq_d Vq_q]
        deltaV = np.empty((len(frequencies),), dtype='cdouble')  # Also dtype='csingle'
        deltaI = np.empty((len(frequencies),), dtype='cdouble')  # Also dtype='csingle'
        names = ['VDUTdc', 'IDUTdcA' + sides]
        for sim, frequency in enumerate(frequencies):
            fft_idx = int(round(frequency * fft_periods * 1 / f_base))  # Index of the target FFT frequency
            for name in names:
                delta = zblocks.perturbation_data[sim][name + "_dc"][start_idx:] - zblocks.snapshot_data[name][start_idx:]
                delta_FD = np.fft.rfft(delta, n=L, axis=0) * 2 / L
                # Retrieve the response at the target frequency
                if "V" in name:
                    deltaV[sim] = delta_FD[fft_idx]
                else:
                    deltaI[sim] = delta_FD[fft_idx]
        Y = deltaI / deltaV

        if results_folder is not None:
            filename = results_name + '#Y_DC#' + zblocks.name + "-" + sides
            np.savetxt(results_folder + '\\' + filename + '#.txt', np.stack((frequencies, Y), axis=1), delimiter='\t', header="f\t"+zblocks.name+"-"+sides+"_dc", comments='')
            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
            ax[0].scatter(frequencies, 20 * np.log10(np.abs(Y)), marker='o', facecolors='none', edgecolors='b',
                          linewidths=1.5, label=r'$Y_{dc}$')
            # ax[0].set_yscale("log")
            ax[0].set_xscale("log")
            ax[0].set_xlim([frequencies[0], frequencies[-1]])
            ax[0].minorticks_on()
            ax[0].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
            ax[0].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
            ax[0].set_ylabel('Magnitude [dB]')
            ax[0].set_title('DUT admittance ― ' + str(len(frequencies)) + ' scanned frequencies')
            ax[0].legend(loc='upper right', fancybox=True, shadow=True, ncol=2)

            ax[1].scatter(frequencies, np.angle(Y, deg=True), marker='o', facecolors='none', edgecolors='b',
                          linewidths=1.5, label=r'$Y_{dc}$')
            ax[1].set_xscale("log")
            ax[1].set_ylim([-200, 200])
            ax[1].set_yticks([-180, -90, 0, 90, 180])
            ax[1].set_xlim([frequencies[0], frequencies[-1]])
            ax[1].minorticks_on()
            ax[1].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
            ax[1].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
            ax[1].set_ylabel('Phase [°]')
            ax[1].set_xlabel('Frequency [Hz]')
            ax[1].legend(loc='upper right', fancybox=True, shadow=True, ncol=2)
            fig.savefig(results_folder + '\\' + filename + ".pdf", format="pdf", bbox_inches="tight")
            plt.close(fig)

    elif scantype == "ACDC":
        # For each simulation (freq) it contains a 3x3 matrix: [Vdc_dc Vdc_d Vdc_q; Vd_dc Vd_d Vd_q; Vq_dc Vq_d Vq_q]
        deltaV = np.empty((len(frequencies), 3, 3), dtype='cdouble')  # Also dtype='csingle'
        deltaI = np.empty((len(frequencies), 3, 3), dtype='cdouble')
        Y = np.empty((len(frequencies), 3, 3), dtype='cdouble')
        for block_num, block in enumerate(zblocks):
            if block.type == "AC":
                namesAC = ['VDUTac:1', 'VDUTac:2', 'IDUTacA' + sides[block_num] + ':1', 'IDUTacA' + sides[block_num] + ':2']
            else:
                namesDC = ['VDUTdc', 'IDUTdcA' + sides[block_num]]
        names = namesDC + namesAC
        row = {names[0]: 0, names[1]: 0, names[2]: 1, names[3]: 2, names[4]: 1, names[5]: 2}
        for sim, frequency in enumerate(frequencies):
            fft_idx = int(round(frequency * fft_periods * 1 / f_base))  # Index of the target FFT frequency
            for col, sim_type in enumerate(["_dc", "_d", "_q"]):
                for block in zblocks:
                    # Select whether it is a DC or AC block
                    if block.type == "AC":
                        names = namesAC
                    else:
                        names = namesDC
                    for name in names:
                        delta = block.perturbation_data[sim][name + sim_type][start_idx:] - block.snapshot_data[name][start_idx:]
                        delta_FD = np.fft.rfft(delta, n=L, axis=0) * 2 / L
                        # Retrieve the response at the target frequency
                        if "V" in name:
                            deltaV[sim, row[name], col] = delta_FD[fft_idx]
                        else:
                            deltaI[sim, row[name], col] = delta_FD[fft_idx]
            Y[sim, ...] = np.matmul(deltaI[sim, ...], np.linalg.inv(deltaV[sim, ...]))

        if results_folder is not None:
            filename = results_name + '#Y_ACDC#' + zblocks[0].name+"-"+sides[0]+"#"+zblocks[1].name+"-"+sides[1]
            # header = "f\tdc-dc\tdc-d\tdc-q\td-dc\td-d\td-q\tq-dc\tq-d\tq-q"
            if zblocks[0].type == "AC":
                header = "f\t" + zblocks[1].name + "-" + sides[1] + "_dc\t" + zblocks[0].name + "-" + sides[0] + "_d\t" + zblocks[0].name + "-" + sides[0] + "_q"
            else:
                header = "f\t" + zblocks[0].name + "-" + sides[0] + "_dc\t" + zblocks[1].name + "-" + sides[1] + "_d\t" + zblocks[1].name + "-" + sides[1] + "_q"
            np.savetxt(results_folder + '\\' + filename+'#.txt', np.c_[frequencies, Y.reshape(Y.shape[0],-1)], delimiter='\t',header=header,comments='')
            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
            ax[0].scatter(frequencies, 20 * np.log10(np.abs(Y[:, 1, 1])), marker='o', facecolors='none', edgecolors='b',
                          linewidths=1.5, label=r'$Y_{dd}$')
            ax[0].scatter(frequencies, 20 * np.log10(np.abs(Y[:, 1, 2])), marker='x', c='r', linewidths=1.5,
                          label=r'$Y_{dq}$')
            ax[0].scatter(frequencies, 20 * np.log10(np.abs(Y[:, 2, 1])), marker='+', c='m', linewidths=1.5,
                          label=r'$Y_{qd}$')
            ax[0].scatter(frequencies, 20 * np.log10(np.abs(Y[:, 2, 2])), marker='.', c='g', linewidths=1.5,
                          label=r'$Y_{qq}$')
            # ax[0].set_yscale("log")
            ax[0].set_xscale("log")
            ax[0].set_xlim([frequencies[0], frequencies[-1]])
            ax[0].minorticks_on()
            ax[0].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
            ax[0].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
            ax[0].set_ylabel('Magnitude [dB]')
            ax[0].set_title('DUT admittance ― ' + str(len(frequencies)) + ' scanned frequencies')
            ax[0].legend(loc='upper right', ncol=2)

            ax[1].scatter(frequencies, np.angle(Y[:, 1, 1], deg=True), marker='o', facecolors='none', edgecolors='b',
                          linewidths=1.5, label=r'$Y_{dd}$')
            ax[1].scatter(frequencies, np.angle(Y[:, 1, 2], deg=True), marker='x', c='r', linewidths=1.5,
                          label=r'$Y_{dq}$')
            ax[1].scatter(frequencies, np.angle(Y[:, 2, 1], deg=True), marker='+', c='m', linewidths=1.5,
                          label=r'$Y_{qd}$')
            ax[1].scatter(frequencies, np.angle(Y[:, 2, 2], deg=True), marker='.', c='g', linewidths=1.5,
                          label=r'$Y_{qq}$')
            ax[1].set_xscale("log")
            ax[1].set_ylim([-200, 200])
            ax[1].set_yticks([-180, -90, 0, 90, 180])
            ax[1].set_xlim([frequencies[0], frequencies[-1]])
            ax[1].minorticks_on()
            ax[1].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
            ax[1].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
            ax[1].set_ylabel('Phase [°]')
            ax[1].set_xlabel('Frequency [Hz]')
            ax[1].legend(loc='upper right', ncol=2)
            filename_aux = filename.split("#")
            filename_aux[1] = "Y_AC"
            fig.savefig(results_folder + '\\' + "#".join(filename_aux) + ".pdf", format="pdf", bbox_inches="tight")
            plt.close(fig)

            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
            ax[0].scatter(frequencies, 20 * np.log10(np.abs(Y[:, 0, 0])), marker='o', facecolors='none', edgecolors='b',
                          linewidths=1.5, label=r'$Y_{dc}$')
            # ax[0].set_yscale("log")
            ax[0].set_xscale("log")
            ax[0].set_xlim([frequencies[0], frequencies[-1]])
            ax[0].minorticks_on()
            ax[0].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
            ax[0].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
            ax[0].set_ylabel('Magnitude [dB]')
            ax[0].set_title('DUT admittance ― ' + str(len(frequencies)) + ' scanned frequencies')
            ax[0].legend(loc='upper right', ncol=2)

            ax[1].scatter(frequencies, np.angle(Y[:, 0, 0], deg=True), marker='o', facecolors='none', edgecolors='b',
                          linewidths=1.5, label=r'$Y_{dc}$')
            ax[1].set_xscale("log")
            ax[1].set_ylim([-200, 200])
            ax[1].set_yticks([-180, -90, 0, 90, 180])
            ax[1].set_xlim([frequencies[0], frequencies[-1]])
            ax[1].minorticks_on()
            ax[1].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
            ax[1].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
            ax[1].set_ylabel('Phase [°]')
            ax[1].set_xlabel('Frequency [Hz]')
            ax[1].legend(loc='upper right', ncol=2)
            filename_aux[1] = "Y_DC"
            fig.savefig(results_folder + '\\' + "#".join(filename_aux) + ".pdf", format="pdf", bbox_inches="tight")
            plt.close(fig)

            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
            ax[0].scatter(frequencies, 20 * np.log10(np.abs(Y[:, 0, 1])), marker='o', facecolors='none', edgecolors='b',
                          linewidths=1.5, label=r'$Y_{dc-d}$')
            ax[0].scatter(frequencies, 20 * np.log10(np.abs(Y[:, 0, 2])), marker='x', c='r', linewidths=1.5,
                          label=r'$Y_{dc-q}$')
            ax[0].scatter(frequencies, 20 * np.log10(np.abs(Y[:, 1, 0])), marker='+', c='m', linewidths=1.5,
                          label=r'$Y_{d-dc}$')
            ax[0].scatter(frequencies, 20 * np.log10(np.abs(Y[:, 2, 0])), marker='.', c='g', linewidths=1.5,
                          label=r'$Y_{q-dc}$')
            # ax[0].set_yscale("log")
            ax[0].set_xscale("log")
            ax[0].set_xlim([frequencies[0], frequencies[-1]])
            ax[0].minorticks_on()
            ax[0].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
            ax[0].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
            ax[0].set_ylabel('Magnitude [dB]')
            ax[0].set_title('DUT admittance ― ' + str(len(frequencies)) + ' scanned frequencies')
            ax[0].legend(loc='upper right', ncol=2)

            ax[1].scatter(frequencies, np.angle(Y[:, 0, 1], deg=True), marker='o', facecolors='none', edgecolors='b',
                          linewidths=1.5, label=r'$Y_{dc-d}$')
            ax[1].scatter(frequencies, np.angle(Y[:, 0, 2], deg=True), marker='x', c='r', linewidths=1.5,
                          label=r'$Y_{dc-q}$')
            ax[1].scatter(frequencies, np.angle(Y[:, 1, 0], deg=True), marker='+', c='m', linewidths=1.5,
                          label=r'$Y_{d-dc}$')
            ax[1].scatter(frequencies, np.angle(Y[:, 2, 0], deg=True), marker='.', c='g', linewidths=1.5,
                          label=r'$Y_{q-dc}$')
            ax[1].set_xscale("log")
            ax[1].set_ylim([-200, 200])
            ax[1].set_yticks([-180, -90, 0, 90, 180])
            ax[1].set_xlim([frequencies[0], frequencies[-1]])
            ax[1].minorticks_on()
            ax[1].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
            ax[1].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
            ax[1].set_ylabel('Phase [°]')
            ax[1].set_xlabel('Frequency [Hz]')
            ax[1].legend(loc='upper right', ncol=2)
            filename_aux[1] = "Y_ACDC_couplings"
            fig.savefig(results_folder + '\\' + "#".join(filename_aux) + ".pdf", format="pdf", bbox_inches="tight")
            plt.close(fig)

    else:
        # NETWORK SCAN
        # Size of the admittance matrix
        if exploit_dq_sym and network.scan_type == "AC":
            N = 2*network.runs
        else:
            N = network.runs 
            
        # For each simulation (freq) a NxN matrix is computed where N = #buses for DC grids or 2 * #buses for AC grids
        deltaV = np.empty((len(frequencies), N, N), dtype='cdouble')  # Also dtype='csingle'
        deltaI = np.empty((len(frequencies), N, N), dtype='cdouble')
        Y = np.empty((len(frequencies), N, N), dtype='cdouble')

        data_ending = ["_"+str(pert) for pert in range(1,N+1)]  # Use the number of runs to define the file ending
        # The rows for each variable are based on the names of all_scans (topology), i.e. zblocks are sorted already
        # If removing the scan_type info is needed: "_".join(scan_name.split("_")[:-1])
        # row = {scan_name: idx for idx, scan_name in enumerate(network.all_scans)}

        # Build the small-signal complex voltage and current matrices at every frequency
        for sim, frequency in enumerate(frequencies):
            fft_idx = int(round(frequency * fft_periods * 1 / f_base))  # Index of the target FFT frequency
            for col, sim_type in enumerate(data_ending):
                current_row = 0
                for block_num, block in enumerate(zblocks):
                    if network.scan_type == "AC":
                        names = ['VDUTac:1','VDUTac:2','IDUTacA'+sides[block_num]+':1','IDUTacA'+sides[block_num]+':2']
                        if exploit_dq_sym and col < N/2:
                            # Exploiting dq-symmetry: # runs is half and only corresponds to d-axis perturbations
                            # The data for the q-axis perturbation is derived from that of the d-axis via symmetry
                            # Only half of the necessary perturbation data is used (exists): d-voltage perturbed
                            for name_pos, name in enumerate(names):
                                delta = block.perturbation_data[sim][name+sim_type][start_idx:] - block.snapshot_data[name][ start_idx:]
                                delta_FD = np.fft.rfft(delta, n=L, axis=0) * 2 / L
                                row = current_row + (name_pos % 2)  # _d variable followed by _q variable
                                # Retrieve the response at the target frequency
                                if "V" in name:
                                    deltaV[sim, row, 2*col] = delta_FD[fft_idx]  # Actual d-axis perturbation
                                else:
                                    deltaI[sim, row, 2*col] = delta_FD[fft_idx]  # Actual d-axis response data
                            
                            # The data for the q-axis is mirrored w.r.t. the d-axis data
                            deltaV[sim, current_row,   2*col+1] = - deltaV[sim, current_row+1, 2*col]
                            deltaV[sim, current_row+1, 2*col+1] = deltaV[sim, current_row,   2*col]
                            deltaI[sim, current_row,   2*col+1] = - deltaI[sim, current_row+1, 2*col]
                            deltaI[sim, current_row+1, 2*col+1] = deltaI[sim, current_row,   2*col]

                        elif not exploit_dq_sym:
                            # Computation using both d-axis and q-axis perturbations (no assumptions on the symmetry)
                            for name_pos, name in enumerate(names):
                                delta = block.perturbation_data[sim][name+sim_type][start_idx:] - block.snapshot_data[name][start_idx:]
                                delta_FD = np.fft.rfft(delta, n=L, axis=0) * 2 / L
                                row = current_row + (name_pos % 2)  # _d variable followed by _q variable
                                # Retrieve the response at the target frequency
                                if "V" in name:
                                    deltaV[sim, row, col] = delta_FD[fft_idx]
                                else:
                                    deltaI[sim, row, col] = delta_FD[fft_idx]
                    else:
                        names = ['VDUTdc', 'IDUTdcA' + sides[block_num]]
                        for name in names:
                            delta = block.perturbation_data[sim][name+sim_type][start_idx:] - block.snapshot_data[name][start_idx:]
                            delta_FD = np.fft.rfft(delta, n=L, axis=0) * 2 / L
                            # Retrieve the response at the target frequency
                            if "V" in name:
                                deltaV[sim, current_row, col] = delta_FD[fft_idx]
                            else:
                                deltaI[sim, current_row, col] = delta_FD[fft_idx]

                    current_row = current_row + int(len(names)/2)

            if network.enforce:
                # Enforce network connectivity
                Yextended = network.adj_matrix  # Matrix indicating the connectivity
                np.fill_diagonal(Yextended, 1)  # Fill the diagonal with ones (shunt or self admittance)
                # if sim == 0:
                #     print("Network connectivity matrix \n",Yextended)
                if network.scan_type == "AC":
                    Yextended = np.kron(Yextended,np.ones((2,2),dtype=int))  # Extend the matrix with dq-axes
                    if sim == 0: print("\n Network dq connectivity matrix \n", Yextended)
                # Yextended[m,n] = 1 <-> y[m,n] =/= 0, so we force the rest to zero
                Y[sim, ...] = np.multiply(Yextended,np.matmul(deltaI[sim, ...], np.linalg.inv(deltaV[sim, ...])))
            else:
                Y[sim, ...] = np.matmul(deltaI[sim, ...], np.linalg.inv(deltaV[sim, ...]))

        if results_folder is not None:
            filename = results_name+'#Y_'+network.scan_type+"#"+"#".join([zblocks[idx].name+"-"+sides[idx] for idx in range(len(sides))])
            results = [Y[:, row, col] for row in range(N) for col in range(N)]
            results.insert(0, frequencies)
            results = tuple(results)

            header = ["f"]  # Populate the header starting with the frequency and then the rest of the variables block by block
            for block_num, block in enumerate(zblocks):
                if block.type == "AC":
                    header.append(block.name+"-"+sides[block_num]+"_d")
                    header.append(block.name+"-"+sides[block_num]+"_q")
                else:
                    header.append(block.name+"-"+sides[block_num]+"_dc")
            header = "\t".join(header)

            np.savetxt(results_folder+r'\\'+filename+'#.txt',np.stack(results, axis=1),  delimiter='\t', header=header,comments='') #, comments="\t".join(network.all_scans)

            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
            for row in range(N):
                for col in range(N):
                    if network.enforce:
                        # Only plot the non-zero elements
                        if Yextended[row,col] == 1:
                            ax[0].scatter(frequencies, 20*np.log10(np.abs(Y[:,row,col])),linewidths=1.0,label=r'$Y_{'+str(row)+str(col)+'}$')
                    else:
                        ax[0].scatter(frequencies, 20 * np.log10(np.abs(Y[:,row,col])),linewidths=1.0,label=r'$Y_{'+str(row)+str(col)+'}$')
            # ax[0].set_yscale("log")
            ax[0].set_xscale("log")
            ax[0].set_xlim([frequencies[0], frequencies[-1]])
            ax[0].minorticks_on()
            ax[0].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
            ax[0].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
            ax[0].set_ylabel('Magnitude [dB]')
            ax[0].set_title('DUT admittance ― ' + str(len(frequencies)) + ' scanned frequencies')
            ax[0].legend(loc='upper right', ncol=4)
            # ,fancybox=True, shadow=True,
            for row in range(N):
                for col in range(N):
                    if network.enforce:
                        if Yextended[row,col] == 1:
                            ax[1].scatter(frequencies,np.angle(Y[:,row,col],deg=True),linewidths=1.0,label=r'$Y_{'+str(row)+str(col)+'}$')
                    else:
                        ax[1].scatter(frequencies,np.angle(Y[:,row,col],deg=True),linewidths=1.0,label=r'$Y_{'+str(row)+str(col)+'}$')
            ax[1].set_xscale("log")
            ax[1].set_ylim([-200, 200])
            ax[1].set_yticks([-180, -90, 0, 90, 180])
            ax[1].set_xlim([frequencies[0], frequencies[-1]])
            ax[1].minorticks_on()
            ax[1].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
            ax[1].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
            ax[1].set_ylabel('Phase [°]')
            ax[1].set_xlabel('Frequency [Hz]')
            fig.savefig(results_folder+'\\'+filename + ".pdf", format="pdf", bbox_inches="tight")
            plt.close(fig)

def admittance_multi_freq(f_base=None, freq_multi=None, fft_periods=1, sides=None, dt=None, exploit_dq_sym=False,
                          start_idx=None, zblocks=None, results_folder=None, results_name='Y', network=None):
    # Small-signal sinusoidal steady state computation and rFFT (no target frequency-based FFT distinction)
    L = int(fft_periods * 1 / f_base * 1.0 / dt)  # For the FFT computation

    # Size of the admittance matrix
    if exploit_dq_sym and network.scan_type == "AC":
        N = 2*network.runs
    else:
        N = network.runs

    frequencies = freq_multi.reshape(np.prod(freq_multi.shape))  # Frequencies
    # TODO: remove zero entries (pad) from frequencies and freq_multi!
    
    # For each frequency a NxN matrix is computed where N = #buses for DC grids or 2 * #buses for AC grids
    deltaV = np.empty((len(frequencies), N, N), dtype='cdouble')  # Also dtype='csingle'
    deltaI = np.empty((len(frequencies), N, N), dtype='cdouble')
    Y = np.empty((len(frequencies), N, N), dtype='cdouble')

    data_ending = ["_"+str(run) for run in network.runs_list]  # Use the number of runs to define the file ending
    # The rows for each variable are based on the names of all_scans (topology), i.e. zblocks are sorted already
    # If removing the scan_type info is needed: "_".join(scan_name.split("_")[:-1])
    # row = {scan_name: idx for idx, scan_name in enumerate(network.all_scans)}

    # Build the small-signal complex voltage and current matrices at every frequency
    for sim in range(freq_multi.shape[1]):
        # For every simulation, read the files and extract the FFT at the corresponding frequencies
        for col, sim_type in enumerate(data_ending):
            current_row = 0
            for block_num, block in enumerate(zblocks):
                if network.scan_type == "AC":
                    names = ['VDUTac:1','VDUTac:2','IDUTacA'+sides[block_num]+':1','IDUTacA'+sides[block_num]+':2']
                    if exploit_dq_sym and col < N/2:
                        # Exploiting dq-symmetry: # runs is half and only corresponds to d-axis perturbations
                        # The data for the q-axis perturbation is derived from that of the d-axis via symmetry
                        # Only half of the necessary perturbation data is used (exists): d-voltage perturbed
                        for name_pos, name in enumerate(names):
                            delta = block.perturbation_data[sim][name+sim_type][start_idx:] - block.snapshot_data[name][ start_idx:]
                            delta_FD = np.fft.rfft(delta, n=L, axis=0) * 2 / L
                            row = current_row + (name_pos % 2)  # _d variable followed by _q variable
                            for freq_file, freq in enumerate(freq_multi[..., sim]):
                                # Retrieve the response at the target frequencies
                                fft_idx = int(round(freq * fft_periods * 1 / f_base))  # Index of the target FFT frequency
                                freq_idx = sim + freq_file*freq_multi.shape[1]  # Index of the current frequency in the vector of all frequencies
                                if "V" in name:
                                    deltaV[freq_idx, row, 2*col] = delta_FD[fft_idx]  # Actual d-axis perturbation
                                else:
                                    deltaI[freq_idx, row, 2*col] = delta_FD[fft_idx]  # Actual d-axis response data

                        for freq_file, freq in enumerate(freq_multi[..., sim]):
                            freq_idx = sim + freq_file*freq_multi.shape[1]  # Index of the current frequency in the vector of all frequencies
                            # The data for the q-axis is mirrored w.r.t. the d-axis data
                            deltaV[freq_idx, current_row,   2*col+1] = - deltaV[freq_idx, current_row+1, 2*col]
                            deltaV[freq_idx, current_row+1, 2*col+1] = deltaV[freq_idx, current_row,   2*col]
                            deltaI[freq_idx, current_row,   2*col+1] = - deltaI[freq_idx, current_row+1, 2*col]
                            deltaI[freq_idx, current_row+1, 2*col+1] = deltaI[freq_idx, current_row,   2*col]

                    elif not exploit_dq_sym:
                        # Computation using both d-axis and q-axis perturbations (no assumptions on the symmetry)
                        for name_pos, name in enumerate(names):
                            delta = block.perturbation_data[sim][name+sim_type][start_idx:] - block.snapshot_data[name][start_idx:]
                            delta_FD = np.fft.rfft(delta, n=L, axis=0) * 2 / L
                            row = current_row + (name_pos % 2)  # _d variable followed by _q variable
                            for freq_file, freq in enumerate(freq_multi[..., sim]):
                                # Retrieve the response at the target frequencies
                                fft_idx = int(round(freq * fft_periods * 1 / f_base))  # Index of the target FFT frequency
                                freq_idx = sim + freq_file*freq_multi.shape[1]  # Index of the current frequency in the vector of all frequencies
                                if "V" in name:
                                    deltaV[freq_idx, row, col] = delta_FD[fft_idx]
                                else:
                                    deltaI[freq_idx, row, col] = delta_FD[fft_idx]
                else:
                    names = ['VDUTdc', 'IDUTdcA' + sides[block_num]]
                    for name_pos, name in enumerate(names):
                        delta = block.perturbation_data[sim][name+sim_type][start_idx:] - block.snapshot_data[name][start_idx:]
                        delta_FD = np.fft.rfft(delta, n=L, axis=0) * 2 / L
                        for freq_file, freq in enumerate(freq_multi[..., sim]):
                            # Retrieve the response at the target frequencies
                            fft_idx = int(round(freq * fft_periods * 1 / f_base))  # Index of the target FFT frequency
                            freq_idx = sim + freq_file*freq_multi.shape[1]  # Index of the current frequency in the vector of all frequencies
                            # Retrieve the response at the target frequency
                            if "V" in name:
                                deltaV[freq_idx, current_row, col] = delta_FD[fft_idx]
                            else:
                                deltaI[freq_idx, current_row, col] = delta_FD[fft_idx]
                
                current_row = current_row + int(len(names)/2)

    for freq_idx in range(len(frequencies)):
        if network.enforce:
            # Enforce network connectivity
            Yextended = network.adj_matrix  # Matrix indicating the connectivity
            np.fill_diagonal(Yextended, 1)  # Fill the diagonal with ones (shunt or self admittance)
            if network.scan_type == "AC":
                Yextended = np.kron(Yextended,np.ones((2,2),dtype=int))  # Extend the matrix with dq-axes
                if freq_idx == 0: print("\n Network dq connectivity matrix \n", Yextended)
            # Yextended[m,n] = 1 <-> y[m,n] =/= 0, so we force the rest to zero
            Y[freq_idx, ...] = np.multiply(Yextended,np.matmul(deltaI[freq_idx, ...], np.linalg.inv(deltaV[freq_idx, ...])))
        else:
            Y[freq_idx, ...] = np.matmul(deltaI[freq_idx, ...], np.linalg.inv(deltaV[freq_idx, ...]))

    if results_folder is not None:
        filename = results_name+'#Y_'+network.scan_type+"#"+"#".join([zblocks[idx].name+"-"+sides[idx] for idx in range(len(sides))])
        results = [Y[:, row, col] for row in range(N) for col in range(N)]
        results.insert(0, frequencies)
        results = tuple(results)

        header = ["f"]  # Populate the header starting with the frequency and then the rest of the variables block by block
        for block_num, block in enumerate(zblocks):
            if block.type == "AC":
                header.append(block.name+"-"+sides[block_num]+"_d")
                header.append(block.name+"-"+sides[block_num]+"_q")
            else:
                header.append(block.name+"-"+sides[block_num]+"_dc")
        header = "\t".join(header)
        
        np.savetxt(results_folder+r'\\'+filename+'#.txt',np.stack(results, axis=1),  delimiter='\t', header=header, comments='') #, comments="\t".join(network.all_scans)

        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
        for row in range(N):
            for col in range(N):
                if network.enforce:
                    # Only plot the non-zero elements
                    if Yextended[row,col] == 1:
                        ax[0].scatter(frequencies, 20*np.log10(np.abs(Y[:,row,col])),linewidths=1.0,label=r'$Y_{'+str(row)+str(col)+'}$')
                else:
                    ax[0].scatter(frequencies, 20 * np.log10(np.abs(Y[:,row,col])),linewidths=1.0,label=r'$Y_{'+str(row)+str(col)+'}$')
        # ax[0].set_yscale("log")
        ax[0].set_xscale("log")
        ax[0].set_xlim([frequencies[0], frequencies[-1]])
        ax[0].minorticks_on()
        ax[0].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
        ax[0].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
        ax[0].set_ylabel('Magnitude [dB]')
        ax[0].set_title('DUT admittance ― ' + str(len(frequencies)) + ' scanned frequencies')
        ax[0].legend(loc='upper right', ncol=4)
        # ,fancybox=True, shadow=True,
        for row in range(N):
            for col in range(N):
                if network.enforce:
                    if Yextended[row,col] == 1:
                        ax[1].scatter(frequencies,np.angle(Y[:,row,col],deg=True),linewidths=1.0,label=r'$Y_{'+str(row)+str(col)+'}$')
                else:
                    ax[1].scatter(frequencies,np.angle(Y[:,row,col],deg=True),linewidths=1.0,label=r'$Y_{'+str(row)+str(col)+'}$')
        ax[1].set_xscale("log")
        ax[1].set_ylim([-200, 200])
        ax[1].set_yticks([-180, -90, 0, 90, 180])
        ax[1].set_xlim([frequencies[0], frequencies[-1]])
        ax[1].minorticks_on()
        ax[1].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
        ax[1].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
        ax[1].set_ylabel('Phase [°]')
        ax[1].set_xlabel('Frequency [Hz]')
        fig.savefig(results_folder+'\\'+filename + ".pdf", format="pdf", bbox_inches="tight")
        plt.close(fig)

def admittance_generic(f_base=None, frequencies=None, fft_periods=1, sides=None, dt=None, exploit_dq_sym=False,
                       start_idx=None, zblocks=None, results_folder=None, results_name='Y', network=None):
    # Small-signal sinusoidal steady state computation and rFFT (no target frequency-based FFT distinction)
    L = int(fft_periods * 1 / f_base * 1.0 / dt)  # For the FFT computation
    
    if exploit_dq_sym and network.scan_type == "AC":
        N = 2*network.runs
    else:
        N = network.runs 
        
    # For each simulation (freq) a NxN matrix is computed where N = # DC buses + 2 * # AC buses
    deltaV = np.empty((len(frequencies), N, N), dtype='cdouble')  # Also dtype='csingle'
    deltaI = np.empty((len(frequencies), N, N), dtype='cdouble')
    Y = np.empty((len(frequencies), N, N), dtype='cdouble')

    data_ending = ["_"+str(run) for run in network.runs_list]  # Use the number of runs to define the file ending
    # The rows for each variable are based on the names of all_scans (topology), i.e. zblocks are sorted already
    # If removing the scan_type info is needed: "_".join(scan_name.split("_")[:-1])
    # row = {scan_name: idx for idx, scan_name in enumerate(network.all_scans)}

    # Build the small-signal complex voltage and current matrices at every frequency
    for sim, frequency in enumerate(frequencies):
        fft_idx = int(round(frequency * fft_periods * 1 / f_base))  # Index of the target FFT frequency
        for col, sim_type in enumerate(data_ending):
            current_row = 0
            for block_num, block in enumerate(zblocks):
                if block.type == "AC":
                    names = ['VDUTac:1','VDUTac:2','IDUTacA'+sides[block_num]+':1','IDUTacA'+sides[block_num]+':2']
                    if exploit_dq_sym and col < N/2:
                        # Exploiting dq-symmetry: # runs is half and only corresponds to d-axis perturbations
                        # The data for the q-axis perturbation is derived from that of the d-axis via symmetry
                        # Only half of the necessary perturbation data is used (exists): d-voltage perturbed
                        for name_pos, name in enumerate(names):
                            delta = block.perturbation_data[sim][name+sim_type][start_idx:] - block.snapshot_data[name][ start_idx:]
                            delta_FD = np.fft.rfft(delta, n=L, axis=0) * 2 / L
                            row = current_row + (name_pos % 2)  # _d variable followed by _q variable
                            # Retrieve the response at the target frequency
                            if "V" in name:
                                deltaV[sim, row, 2*col] = delta_FD[fft_idx]  # Actual d-axis perturbation
                            else:
                                deltaI[sim, row, 2*col] = delta_FD[fft_idx]  # Actual d-axis response data

                        # The data for the q-axis is mirrored w.r.t. the d-axis data
                        deltaV[sim, current_row,   2*col+1] = - deltaV[sim, current_row+1, 2*col]
                        deltaV[sim, current_row+1, 2*col+1] = deltaV[sim, current_row,   2*col]
                        deltaI[sim, current_row,   2*col+1] = - deltaI[sim, current_row+1, 2*col]
                        deltaI[sim, current_row+1, 2*col+1] = deltaI[sim, current_row,   2*col]

                    elif not exploit_dq_sym:
                        # Computation using both d-axis and q-axis perturbations (no assumptions on the symmetry)
                        for name_pos, name in enumerate(names):
                            delta = block.perturbation_data[sim][name+sim_type][start_idx:] - block.snapshot_data[name][start_idx:]
                            delta_FD = np.fft.rfft(delta, n=L, axis=0) * 2 / L
                            row = current_row + (name_pos % 2)  # _d variable followed by _q variable
                            # Retrieve the response at the target frequency
                            if "V" in name:
                                deltaV[sim, row, col] = delta_FD[fft_idx]
                            else:
                                deltaI[sim, row, col] = delta_FD[fft_idx]
                
                else:
                    names = ['VDUTdc', 'IDUTdcA' + sides[block_num]]
                    for name in names:
                        delta = block.perturbation_data[sim][name+sim_type][start_idx:] - block.snapshot_data[name][start_idx:]
                        delta_FD = np.fft.rfft(delta, n=L, axis=0) * 2 / L
                        # Retrieve the response at the target frequency
                        if "V" in name:
                            deltaV[sim, current_row, col] = delta_FD[fft_idx]
                        else:
                            deltaI[sim, current_row, col] = delta_FD[fft_idx]
                    
                current_row = current_row + int(len(names)/2)

        if network.enforce:
            # Enforce network connectivity
            Yextended = network.adj_matrix  # Matrix indicating the connectivity
            np.fill_diagonal(Yextended, 1)  # Fill the diagonal with ones (shunt or self admittance)
            # if sim == 0:
            #     print("Network connectivity matrix \n",Yextended)
            if network.scan_type == "AC":
                Yextended = np.kron(Yextended,np.ones((2,2),dtype=int))  # Extend the matrix with dq-axes
                if sim == 0: print("\n Network dq connectivity matrix \n", Yextended)
            # Yextended[m,n] = 1 <-> y[m,n] =/= 0, so we force the rest to zero
            Y[sim, ...] = np.multiply(Yextended,np.matmul(deltaI[sim, ...], np.linalg.inv(deltaV[sim, ...])))
        else:
            Y[sim, ...] = np.matmul(deltaI[sim, ...], np.linalg.inv(deltaV[sim, ...]))

    if results_folder is not None:
        filename = results_name+'#Y_'+network.scan_type+"#"+"#".join([zblocks[idx].name+"-"+sides[idx] for idx in range(len(sides))])
        results = [Y[:, row, col] for row in range(N) for col in range(N)]
        results.insert(0, frequencies)
        results = tuple(results)

        header = ["f"]  # Populate the header starting with the frequency and then the rest of the variables block by block
        for block_num, block in enumerate(zblocks):
            if block.type == "AC":
                header.append(block.name+"-"+sides[block_num]+"_d")
                header.append(block.name+"-"+sides[block_num]+"_q")
            else:
                header.append(block.name+"-"+sides[block_num]+"_dc")
        header = "\t".join(header)

        np.savetxt(results_folder+r'\\'+filename+'#.txt',np.stack(results, axis=1),  delimiter='\t', header=header,comments='') #, comments="\t".join(network.all_scans)

        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
        for row in range(N):
            for col in range(N):
                if network.enforce:
                    # Only plot the non-zero elements
                    if Yextended[row,col] == 1:
                        ax[0].scatter(frequencies, 20*np.log10(np.abs(Y[:,row,col])),linewidths=1.0,label=r'$Y_{'+str(row)+str(col)+'}$')
                else:
                    ax[0].scatter(frequencies, 20 * np.log10(np.abs(Y[:,row,col])),linewidths=1.0,label=r'$Y_{'+str(row)+str(col)+'}$')

        ax[0].set_xscale("log")
        ax[0].set_xlim([frequencies[0], frequencies[-1]])
        ax[0].minorticks_on()
        ax[0].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
        ax[0].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
        ax[0].set_ylabel('Magnitude [dB]')
        ax[0].set_title('DUT admittance ― ' + str(len(frequencies)) + ' scanned frequencies')
        ax[0].legend(loc='upper right', ncol=4)
        # ,fancybox=True, shadow=True,
        for row in range(N):
            for col in range(N):
                if network.enforce:
                    if Yextended[row,col] == 1:
                        ax[1].scatter(frequencies,np.angle(Y[:,row,col],deg=True),linewidths=1.0,label=r'$Y_{'+str(row)+str(col)+'}$')
                else:
                    ax[1].scatter(frequencies,np.angle(Y[:,row,col],deg=True),linewidths=1.0,label=r'$Y_{'+str(row)+str(col)+'}$')
        ax[1].set_xscale("log")
        ax[1].set_ylim([-200, 200])
        ax[1].set_yticks([-180, -90, 0, 90, 180])
        ax[1].set_xlim([frequencies[0], frequencies[-1]])
        ax[1].minorticks_on()
        ax[1].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
        ax[1].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
        ax[1].set_ylabel('Phase [°]')
        ax[1].set_xlabel('Frequency [Hz]')
        fig.savefig(results_folder+'\\'+filename + ".pdf", format="pdf", bbox_inches="tight")
        plt.close(fig)

def SISO_TF(f_base=None, frequencies=None, fft_periods=1, dt=None,
            start_idx=None, zblocks=None, results_folder=None, results_name='SISO_TF'):
    # Small-signal sinusoidal steady state computation and rFFT (no target frequency-based FFT distinction)
    L = int(fft_periods * 1 / f_base * 1.0 / dt)  # For the FFT computation

    # For each simulation (freq) a transfer function is computed
    delta_input = np.empty(len(frequencies), dtype='cdouble')  # Also dtype='csingle'
    delta_output= np.empty(len(frequencies), dtype='cdouble')
    Y = np.empty(len(frequencies), dtype='cdouble')

    # Build the small-signal complex voltage and current matrices at every frequency
    for sim, frequency in enumerate(frequencies):
        fft_idx = int(round(frequency * fft_periods * 1 / f_base))  # Index of the target FFT frequency
        for block in zblocks:
            for name in ['inputTF', 'outputTF']:
                delta = block.perturbation_data[sim][name+"_TF"][start_idx:] - block.snapshot_data[name][start_idx:]
                delta_FD = np.fft.rfft(delta, n=L, axis=0) * 2 / L
                # Retrieve the response at the target frequency
                if "input" in name:
                    delta_input[sim] = delta_FD[fft_idx]
                else:
                    delta_output[sim] = delta_FD[fft_idx]

            Y[sim] =delta_output[sim] / delta_input[sim]

    if results_folder is not None:
        filename = results_name+'#SISO_TF_'+"#"+zblocks[0].name
        results = [Y]
        results.insert(0, frequencies)
        results = tuple(results)
        header = "f\t"+"\t"+zblocks[0].name
        np.savetxt(results_folder+r'\\'+filename+'#.txt',np.stack(results, axis=1),  delimiter='\t', header=header,comments='')

        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
        ax[0].scatter(frequencies, 20 * np.log10(np.abs(Y)),linewidths=1.0,label=str(zblocks[0].name))
        ax[0].set_xscale("log")
        ax[0].set_xlim([frequencies[0], frequencies[-1]])
        ax[0].minorticks_on()
        ax[0].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
        ax[0].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
        ax[0].set_ylabel('Magnitude [dB]')
        ax[0].set_title('SISO Transfer function ― ' + str(len(frequencies)) + ' scanned frequencies')
        ax[0].legend(loc='upper right', ncol=4)
        
        ax[1].scatter(frequencies,np.angle(Y,deg=True),linewidths=1.0,label=str(zblocks[0].name))
        ax[1].set_xscale("log")
        ax[1].set_ylim([-200, 200])
        ax[1].set_yticks([-180, -90, 0, 90, 180])
        ax[1].set_xlim([frequencies[0], frequencies[-1]])
        ax[1].minorticks_on()
        ax[1].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
        ax[1].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
        ax[1].set_ylabel('Phase [°]')
        ax[1].set_xlabel('Frequency [Hz]')
        fig.savefig(results_folder+'\\'+filename + ".pdf", format="pdf", bbox_inches="tight")
        plt.close(fig)


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