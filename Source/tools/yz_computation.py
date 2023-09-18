__all__ = ['admittance']

import matplotlib.pyplot as plt
import numpy as np  # Numerical python functions
from matplotlib import rcParams  # Text's parameters for plots

rcParams['mathtext.fontset'] = 'cm'  # Font selection
rcParams['font.family'] = 'STIXGeneral'  # 'cmu serif'

def admittance(f_base=None, frequencies=None, fft_periods=1, scantype="AC", sides=None, dt=None,
               start_idx=None, zblocks=None, results_folder=None, results_name='Y', network=None):
    # Small-signal sinusoidal steady state computation and rFFT (no target frequency-based FFT distinction)
    L = int(fft_periods * 1 / f_base * 1.0 / dt)  # For the FFT computation
    if scantype is "AC":
        # Compute the small-signal sinusoidal steady-state waveforms
        # For each simulation (freq) it contains a matrix: [Vd_d Vd_q; Vq_d Vq_q]
        deltaV = np.empty((len(frequencies), 2, 2), dtype='cdouble')  # Also dtype='csingle'
        deltaI = np.empty((len(frequencies), 2, 2), dtype='cdouble')  # Also dtype='csingle'
        Y = np.empty((len(frequencies), 2, 2), dtype='cdouble')  # Also dtype='csingle'
        names = ['VDUTac:1', 'VDUTac:2', 'IDUTacA' + sides + ':1', 'IDUTacA' + sides + ':2']
        row = {names[0]: 0, names[1]: 1, names[2]: 0, names[3]: 1}
        for sim, frequency in enumerate(frequencies):
            idx = int(round(frequency * fft_periods * 1 / f_base))  # Index of the target frequency
            for col, sim_type in enumerate(["_d", "_q"]):
                # Debug
                # if sim == 21:
                #     for name in names:
                #         # Current and voltages
                #         fig, ax = plt.subplots()
                #         t = np.arange(len(zblocks.perturbation_data[sim][name + sim_type][start_idx:])) * dt
                #         ax.plot(t, zblocks.snapshot_data[name][start_idx:], 'b', linewidth=2, label="Snapshot")
                #         ax.plot(t, zblocks.perturbation_data[sim][name + sim_type][start_idx:], 'r', linewidth=2,
                #                 label="Perturbation")
                #         ax.set(xlabel='Time (s)', ylabel=name)
                #         ax.grid()
                #         ax.legend(loc='upper right', ncol=1)
                #         fig.savefig("test.png")
                #         plt.show()
                # End degub

                for name in names:
                    delta = zblocks.perturbation_data[sim][name+sim_type][start_idx:] - zblocks.snapshot_data[name][start_idx:]
                    delta_FD = np.fft.rfft(delta, n=L, axis=0) * 2 / L
                    # freq_FD = np.fft.rfftfreq(L, d=dt)  # rFFT frequency points
                    # Retrieve the response at the target frequency
                    if "V" in name:
                        deltaV[sim, row[name], col] = delta_FD[idx]
                    else:
                        deltaI[sim, row[name], col] = delta_FD[idx]
            Y[sim,...] = np.matmul(deltaI[sim,...], np.linalg.inv(deltaV[sim,...]))

        # Debugging
        # FFT
        # sim = 20
        # for sim_type in ["_d", "_q"]:
        #     fig, ax = plt.subplots(nrows=2, ncols=1)
        #     colors = {names[0]: 'b', names[1]: 'r', names[2]: 'm', names[3]: 'g'}
        #     for name in names:
        #         delta = zblocks.perturbation_data[sim][name + sim_type][start_idx:] - zblocks.snapshot_data[name][start_idx:]
        #         delta_FD = np.fft.rfft(delta, n=L, axis=0) * 2 / L
        #         ax[0].plot(freq_FD, np.abs(delta_FD), colors[name], linewidth=2, label=name)
        #         ax[1].plot(freq_FD, np.angle(delta_FD, deg=True), colors[name], linewidth=2, label=name)
        #     ax[0].legend(loc='upper right', ncol=1)
        #     ax[0].set_title(zblocks.name+' - Sim type: '+sim_type+' perturbation')
        #     ax[0].set(xlabel='Frequency [Hz]', ylabel="Magnitude")
        #     ax[0].set_xscale("log")
        #     ax[0].set_yscale("log")
        #     ax[0].grid()
        #     ax[0].minorticks_on()
        #     ax[0].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
        #     ax[0].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-',
        #                linewidth=0.5)
        #     ax[1].set(xlabel='Frequency [Hz]', ylabel="Phase [°]")
        #     ax[1].set_xscale("log")
        #     ax[1].grid()
        #     ax[1].minorticks_on()
        #     ax[1].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
        #     ax[1].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-',
        #                linewidth=0.5)
        #     fig.savefig("FFT.png")
            #plt.show()
        # End degub

        if results_folder is not None:
            filename = results_folder + '\\' + results_name + '_Y_' + zblocks.name + sides
            np.savetxt(filename+'.txt',
                       np.stack((frequencies, Y[:, 0, 0], Y[:, 0, 1], Y[:, 1, 0], Y[:, 1, 1]), axis=1), delimiter='\t',
                       header="f\td-d\td-q\tq-d\tq-q", comments='')

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
            ax[0].legend(loc='upper right',fancybox=True, shadow=True, ncol=2)

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
            fig.savefig(filename + ".pdf", format="pdf", bbox_inches="tight")
            fig.clear()

    elif scantype is "DC":
        # Use the actual steady-state waveforms (allows to remove periodicity)
        # For each simulation (freq) it contains a matrix: [Vd_d Vd_q; Vq_d Vq_q]
        deltaV = np.empty((len(frequencies),), dtype='cdouble')  # Also dtype='csingle'
        deltaI = np.empty((len(frequencies),), dtype='cdouble')  # Also dtype='csingle'
        names = ['VDUTdc', 'IDUTdcA' + sides]
        for sim, frequency in enumerate(frequencies):
            idx = int(round(frequency * fft_periods * 1 / f_base))  # Index of the target frequency
            for name in names:
                delta = zblocks.perturbation_data[sim][name + "_dc"][start_idx:] - zblocks.snapshot_data[name][start_idx:]
                delta_FD = np.fft.rfft(delta, n=L, axis=0) * 2 / L
                # Retrieve the response at the target frequency
                if "V" in name:
                    deltaV[sim] = delta_FD[idx]
                else:
                    deltaI[sim] = delta_FD[idx]
        Y = deltaI / deltaV

        if results_folder is not None:
            filename = results_folder + '\\' + results_name + '_Ydc_' + zblocks.name + sides
            np.savetxt(filename+'.txt', np.stack((frequencies, Y), axis=1), delimiter='\t', header="f\tdc", comments='')
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
            fig.savefig(filename + ".pdf", format="pdf", bbox_inches="tight")
            fig.clear()

    elif scantype is "ACDC":
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
            idx = int(round(frequency * fft_periods * 1 / f_base))  # Index of the target frequency
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
                            deltaV[sim, row[name], col] = delta_FD[idx]
                        else:
                            deltaI[sim, row[name], col] = delta_FD[idx]
            Y[sim, ...] = np.matmul(deltaI[sim, ...], np.linalg.inv(deltaV[sim, ...]))

        if results_folder is not None:
            filename = results_folder + '\\' + results_name + '_Yacdc_' + zblocks[0].name + sides[0] + '_' + zblocks[
                1].name + sides[1]
            np.savetxt(filename+'.txt',
                       np.stack((frequencies, Y[:, 0, 0], Y[:, 0, 1], Y[:, 0, 2], Y[:, 1, 0], Y[:, 1, 1], Y[:, 1, 2],
                                 Y[:, 2, 0], Y[:, 2, 1], Y[:, 2, 2]), axis=1),  delimiter='\t',
                       header="f\tdc-dc\tdc-d\tdc-q\td-dc\td-d\td-q\tq-dc\tq-d\tq-q", comments='')

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
            ax[0].legend(loc='upper right',fancybox=True, shadow=True, ncol=2)

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
            ax[1].legend(loc='upper right',fancybox=True, shadow=True, ncol=2)
            fig.savefig(
                results_folder + '\\' + results_name + '_Yac_' + zblocks[0].name + sides[0] + '_' + zblocks[1].name +
                sides[1] + ".pdf",
                format="pdf", bbox_inches="tight")
            fig.clear()

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
            ax[0].legend(loc='upper right',fancybox=True, shadow=True, ncol=2)

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
            ax[1].legend(loc='upper right',fancybox=True, shadow=True, ncol=2)
            fig.savefig(
                results_folder + '\\' +results_name+'_Ydc_' + zblocks[0].name + sides[0] + '_' + zblocks[1].name + sides[1] + ".pdf",
                format="pdf", bbox_inches="tight")
            fig.clear()

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
            ax[0].legend(loc='upper right',fancybox=True, shadow=True, ncol=2)

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
            ax[1].legend(loc='upper right',fancybox=True, shadow=True, ncol=2)
            fig.savefig(
                results_folder + '\\' +results_name+'_Ycoup_' + zblocks[0].name + sides[0] + '_' + zblocks[1].name + sides[
                    1] + ".pdf",
                format="pdf", bbox_inches="tight")

    else:
        # Passive network scan
        N = network.runs  # Size of the admittance matrix
        # For each simulation (freq) a NxN matrix is computed where N = #buses for DC grids or 2 * #buses for AC grids
        deltaV = np.empty((len(frequencies), N, N), dtype='cdouble')  # Also dtype='csingle'
        deltaI = np.empty((len(frequencies), N, N), dtype='cdouble')
        Y = np.empty((len(frequencies), N, N), dtype='cdouble')

        sim_type_ending = ["_"+str(pert) for pert in range(1,N+1)]  # Use the number of runs to define the file ending
        # The rows for each variable are based on the names of all_scans (topology), i.e. zblocks are sorted already
        # If removing the scan_type info is needed: "_".join(scan_name.split("_")[:-1])
        # row = {scan_name: idx for idx, scan_name in enumerate(network.all_scans)}

        for sim, frequency in enumerate(frequencies):
            idx = int(round(frequency * fft_periods * 1 / f_base))  # Index of the target frequency
            for col, sim_type in enumerate(sim_type_ending):
                for block_num, block in enumerate(zblocks):
                    if network.scan_type == "AC":
                        names = ['VDUTac:1','VDUTac:2','IDUTacA'+sides[block_num]+':1','IDUTacA'+sides[block_num]+':2']
                        for name in names:
                            delta = block.perturbation_data[sim][name+sim_type][start_idx:] - block.snapshot_data[name][start_idx:]
                            delta_FD = np.fft.rfft(delta, n=L, axis=0) * 2 / L
                            row = 2*block_num + int(name[-1]) - 1  # _d variable followed by _q variable
                            # Retrieve the response at the target frequency
                            if "V" in name:
                                deltaV[sim, row, col] = delta_FD[idx]
                            else:
                                deltaI[sim, row, col] = delta_FD[idx]
                    else:
                        names = ['VDUTdc', 'IDUTdcA' + sides[block_num]]
                        for name in names:
                            delta = block.perturbation_data[sim][name+sim_type][start_idx:] - block.snapshot_data[name][start_idx:]
                            delta_FD = np.fft.rfft(delta, n=L, axis=0) * 2 / L
                            # Retrieve the response at the target frequency
                            if "V" in name:
                                deltaV[sim, block_num, col] = delta_FD[idx]
                            else:
                                deltaI[sim, block_num, col] = delta_FD[idx]

            # if sim == 0:
            #     print("\ndeltaI\n", deltaI[sim, ...])
            #     print("deltaV\n", deltaV[sim, ...])
            #     print("Topology adjacent matrix\n", network.adj_matrix)

            if network.enforce:
                # Enforce network connectivity
                if network.scan_type == "AC":
                    Yextended = np.kron(network.adj_matrix,np.ones((2,2),dtype=int))  # Extend the matrix with dq-axes
                else:
                    Yextended = network.adj_matrix  # No extension needed
                # Yextended[m,n] = 1 <-> y[m,n] =/= 0, so we force the rest to zero
                Y[sim, ...] = np.matmul(Yextended * deltaI[sim, ...], np.linalg.inv(deltaV[sim, ...]))
            else:
                Y[sim, ...] = np.matmul(deltaI[sim, ...], np.linalg.inv(deltaV[sim, ...]))

        if results_folder is not None:
            filename = results_folder+'\\'+results_name+'_Y_'+network.scan_type+"_".join([zblocks[idx].name+sides[idx] for idx in range(len(sides))])
            results = [Y[:, row, col] for row in range(N) for col in range(N)]
            results.insert(0, frequencies)
            elements = [str(row)+"-"+str(col) for row in range(N) for col in range(N)]
            results = tuple(results)
            header = "f\t"+"\t".join(elements)
            np.savetxt(filename+'.txt',np.stack(results, axis=1),  delimiter='\t', header=header, comments="\n"+"\t".join(network.all_scans))

            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
            for row in range(N):
                for col in range(N):
                    if network.enforce:
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
            ax[0].legend(loc='upper right',fancybox=True, shadow=True, ncol=4)

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
            fig.savefig(filename + ".pdf", format="pdf", bbox_inches="tight")
            fig.clear()