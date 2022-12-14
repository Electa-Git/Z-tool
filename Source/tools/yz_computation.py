__all__ = ['admittance']

import matplotlib.pyplot as plt
import numpy as np  # Numerical python functions
from matplotlib import rcParams  # Text's parameters for plots
rcParams['mathtext.fontset'] = 'cm'  # Font selection
rcParams['font.family'] = 'STIXGeneral'  # 'cmu serif'

def admittance(f_base=None, frequencies=None, fft_periods=1, start_fft=None,
               ss=None, vi1_td=None, vi2_td=None, td=None, results_folder=None, results_name='Y'):
    dt = np.mean([td[i + 1] - td[i] for i in range(min(len(td), 100))])  # Sampling time [s]

    # Subtracts the steady state values from the signals after the system reached steady state
    start_idx = find_nearest(td, start_fft)
    if ss is None:
        # Steady-state as the average after the transient period (assuming a LTI sys i.e. no LTP = no steady-state osc.)
        # WARNING: THIS IS NOT VALIDATED YET
        deltavi1_td = vi1_td[start_idx:, :] - np.tile(np.mean(vi1_td[start_idx:, :], axis=0), (td[start_idx:].size, 1))
        deltavi2_td = vi2_td[start_idx:, :] - np.tile(np.mean(vi2_td[start_idx:, :], axis=0), (td[start_idx:].size, 1))
    else:
        # Use the steady-state waveforms
        if td[0] != 0.0 and ss[0, 0] != td[0]: ss = ss[1:, :]  # Shift the snapshot by dt
        ss_ext = np.tile(ss[start_idx:, 1:], (1, int(round(vi1_td.shape[1] / 4))))  # Extend the steady-state matrix
        deltavi1_td = vi1_td[start_idx:, :] - ss_ext  # Small-signal steady state computation
        deltavi2_td = vi2_td[start_idx:, :] - ss_ext

    Y_dd = np.empty((len(frequencies),), dtype='cdouble')  # Also dtype='csingle'
    Y_dq = np.empty((len(frequencies),), dtype='cdouble')
    Y_qd = np.empty((len(frequencies),), dtype='cdouble')
    Y_qq = np.empty((len(frequencies),), dtype='cdouble')

    # Option 1: No target frequency-based FFT distinction
    L = int(fft_periods * 1 / f_base * 1.0 / dt)
    deltavi1_fd = np.fft.rfft(deltavi1_td, n=L, axis=0) * 2 / L
    deltavi2_fd = np.fft.rfft(deltavi2_td, n=L, axis=0) * 2 / L
    # freqs = np.fft.rfftfreq(L, d=dt) # rFFT frequency points
    for sim, frequency in enumerate(frequencies):
        idx = int(round(frequency * fft_periods * 1 / f_base))
        Vdq = np.matrix([[deltavi1_fd[idx, 4 * sim + 0], deltavi2_fd[idx, 4 * sim + 0]],
                         [deltavi1_fd[idx, 4 * sim + 1], deltavi2_fd[idx, 4 * sim + 1]]])
        Idq = np.matrix([[deltavi1_fd[idx, 4 * sim + 2], deltavi2_fd[idx, 4 * sim + 2]],
                         [deltavi1_fd[idx, 4 * sim + 3], deltavi2_fd[idx, 4 * sim + 3]]])

        Y = Idq * np.linalg.inv(Vdq)
        Y_dd[sim] = Y[0, 0]
        Y_dq[sim] = Y[0, 1]
        Y_qd[sim] = Y[1, 0]
        Y_qq[sim] = Y[1, 1]

        # Z = Vdq * np.linalg.inv(Idq)
        # Z_dd[sim] = Z[0,0]
        # Z_dq[sim] = Z[0,1]
        # Z_qd[sim] = Z[1,0]
        # Z_qq[sim] = Z[1,1]

    # Option 2: Target frequency-based FFT distinction (more efficient) To be added

    if results_folder is not None:
        np.savetxt(results_folder + '\\' + results_name + '_Y.txt',
                   np.stack((frequencies, Y_dd, Y_dq, Y_qd, Y_qq), axis=1), header="f\tdd\tdq\tqd\tqq", delimiter='\t',
                   comments='')

        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
        ax[0].scatter(frequencies, 20*np.log10(np.abs(Y_dd)), marker='o', facecolors='none', edgecolors='b',
                      linewidths=1.5, label=r'$Y_{dd}$')
        ax[0].scatter(frequencies, 20*np.log10(np.abs(Y_dq)), marker='x', c='r', linewidths=1.5, label=r'$Y_{dq}$')
        ax[0].scatter(frequencies, 20*np.log10(np.abs(Y_qd)), marker='+', c='m', linewidths=1.5, label=r'$Y_{qd}$')
        ax[0].scatter(frequencies, 20*np.log10(np.abs(Y_qq)), marker='.', c='g', linewidths=1.5, label=r'$Y_{qq}$')
        # ax[0].set_yscale("log")
        ax[0].set_xscale("log")
        ax[0].set_xlim([frequencies[0], frequencies[-1]])
        ax[0].minorticks_on()
        ax[0].grid(visible=True, which='major', color='k', linestyle='-', linewidth=0.5)
        ax[0].grid(visible=True, which='minor', color='tab:gray', alpha=0.5, linestyle='-', linewidth=0.5)
        ax[0].set_ylabel('Magnitude [dB]')
        ax[0].set_title('DUT admittance ― ' + str(len(frequencies)) + ' scanned frequencies')
        ax[0].legend(loc='upper right', ncol=2)

        ax[1].scatter(frequencies, np.angle(Y_dd, deg=True), marker='o', facecolors='none', edgecolors='b',
                      linewidths=1.5, label=r'$Y_{dd}$')
        ax[1].scatter(frequencies, np.angle(Y_dq, deg=True), marker='x', c='r', linewidths=1.5, label=r'$Y_{dq}$')
        ax[1].scatter(frequencies, np.angle(Y_qd, deg=True), marker='+', c='m', linewidths=1.5, label=r'$Y_{qd}$')
        ax[1].scatter(frequencies, np.angle(Y_qq, deg=True), marker='.', c='g', linewidths=1.5, label=r'$Y_{qq}$')
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
        fig.savefig(results_folder + '\\' + results_name + "_Y.pdf", format="pdf", bbox_inches="tight")

def find_nearest(array, value):  # Efficient function to find the nearest value to a given one and its position
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or np.abs(value - array[idx - 1]) < np.abs(value - array[idx])):
        return idx - 1
    else:
        return idx