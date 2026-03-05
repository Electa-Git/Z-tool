"""
This program contains several useful functions for frequency-domain analysis such as spectrum and Bode plots.

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

__all__ = ['bode_plot', 'spectrum_plot']

import matplotlib.pyplot as plt
import numpy as np  # Numerical python functions
from matplotlib import rcParams  # Text's parameters for plots
import pickle

rcParams['mathtext.fontset'] = 'cm'  # Font selection
rcParams['font.family'] = 'STIXGeneral'  # 'cmu serif'

def spectrum_plot(signals=None, time=None, results_folder=None, file_name='spectrum', labels=None, style='line', show_plot=False,  save_data=False, show_dc=False, save_pickle=False):
    """This function produces a spectrum plot of the signals given as argument."""
    if results_folder is None and not show_plot:
        raise ValueError("The destination folder must be provided via the 'results_folder' argument")
    if signals is None or time is None:
        raise ValueError("The list of signals and time must be provided via the 'signals' and 'time' arguments")
    
    # Compute the FFT
    N = signals.shape[0] # Number of sample points
    T = time[1] - time[0]  # Sampling period [s]
    yf = np.fft.rfft(signals, axis=0) # FFT of the signal
    xf = np.fft.fftfreq(N, T)[:N//2]  # Frequencies

    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 4))
    for signal_col in range(signals.shape[1]):
        label = r'$'+str(signal_col)+'$' if labels is None else labels[signal_col]
        label = label + ', DC: ' + format(1.0/N*np.abs(yf[0, signal_col]), f".{4}e")
        ax[0].plot(time, signals[:, signal_col], linewidth=1.0, label=label)
        if style == 'line':
            ax[1].plot(xf, 2.0/N*np.abs(yf[0:N//2, signal_col]), linewidth=1.0, label=label)
        else:
            ax[1].scatter(xf, 2.0/N*np.abs(yf[0:N//2, signal_col]), linewidths=1, label=label)

    ax[0].set_xlim([time[0], time[-1]])
    ax[0].minorticks_on()
    ax[0].grid(visible=True, which='major', color='k', linestyle='-', alpha=0.5, linewidth=0.5)
    ax[0].set_ylabel('Signal')
    ax[0].set_xlabel('Time [s]')
    ax[0].set_title('Time domain and FFT plots')
    
    ax[1].set_yscale("log")
    if show_dc:
        ax[1].set_xlim([xf[0], xf[-1]])
    else:
        ax[1].set_xscale("log")
        ax[1].set_xlim([xf[1], xf[-1]])
    ax[1].minorticks_on()
    ax[1].grid(visible=True, which='major', color='k', linestyle='-', alpha=0.5, linewidth=0.5)
    ax[1].set_ylabel('Amplitude')
    ax[1].set_xlabel('Frequency [Hz]')
    ax[1].legend(loc='lower right', ncol=int(np.floor(np.sqrt(np.prod(signals.shape[1:])))), prop={'size': 6})    # ,fancybox=True, shadow=True,

    if results_folder is not None:
        fig.savefig(results_folder + '\\' + file_name + ".pdf", format="pdf", bbox_inches="tight")
        if save_pickle:
            with open(results_folder + '\\' + file_name + ".pickle", 'wb') as f: pickle.dump(fig, f)
        if save_data:
            results = [yf[0:N//2, col] for col in range(signals.shape[1])]
            results.insert(0, xf)
            results = tuple(results)
            np.savetxt(results_folder+r'\\'+file_name+'#.txt',np.stack(results, axis=1),  delimiter='\t', header='f [Hz]\t'+'\t'.join([str(signal_col) if labels is None else labels[signal_col] for signal_col in range(signals.shape[1])]),comments='')
    if show_plot: plt.show()  # Visualize the plot 
    plt.close(fig)

def bode_plot(Y=None, frequencies=None, results_folder=None, file_name='Bode_plot', style='scatter', show_plot=False, legend=[], save_data=False, return_plot=False, fig_handle=None, title=None, save_pickle=False, linear_mag=False):
    """This function produces a Bode plot of the complex matrix and frequency list given as argument."""
    if results_folder is None and (not show_plot and not return_plot):
        raise ValueError("The destination folder must be provided via the 'results_folder' argument")
    if frequencies is None:
        raise ValueError("The list of frequencies must be provided via the 'frequencies' argument")
    
    rows, cols = Y[0].shape if Y.ndim == 3 else (max(Y[0].shape if len(Y[0].shape) > 0 else [1]), 1)    
    if rows == cols and rows > 1:
        # Square matrix: generate pairwise labels' combinations
        if legend is not None and len(legend) == rows:
            labels = [ [legend[i]+", "+legend[j] for j in range(rows)] for i in range(rows) ]
        elif legend is not None:
            labels = [ [str(i)+", "+str(j) for j in range(rows)] for i in range(rows) ]
    else:
        # Otherwise, the labels correspond to the rows of the matrix
        if legend is not None and len(legend) == rows:
            labels = [legend[i] for i in range(rows)]
        elif legend is not None:
            labels = [str(i) for i in range(rows)]


    if fig_handle is not None:
        fig, ax = fig_handle # Retieve the figure and axes objects from the provided handle
    else:
        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 6)) # Create the figure and axes

    if Y.ndim == 3: # MIMO
        for row in range(Y.shape[1]):
            for col in range(Y.shape[2]):
                if np.abs(Y[0,row,col]) != 0: # Only plot if non-zero
                    if style == 'scatter':
                        ax[0].scatter(frequencies, 20*np.log10(np.abs(Y[:,row,col])) if not linear_mag else np.abs(Y[:,row,col]),linewidths=1.0,label=r'$'+labels[row][col]+r'$' if legend is not None else '_nolegend_')
                        ax[1].scatter(frequencies, np.angle(Y[:,row,col],deg=True),linewidths=1.0,label=r'$'+labels[row][col]+r'$' if legend is not None else '_nolegend_')
                    else:
                        ax[0].plot(frequencies, 20*np.log10(np.abs(Y[:,row,col])) if not linear_mag else np.abs(Y[:,row,col]),linewidth=2.0,label=r'$'+labels[row][col]+r'$' if legend is not None else '_nolegend_')
                        ax[1].plot(frequencies, np.angle(Y[:,row,col],deg=True),linewidth=2.0,label=r'$'+labels[row][col]+r'$' if legend is not None else '_nolegend_')
    elif Y.ndim == 2: # SIMO
        for row in range(Y.shape[1]):
            if np.abs(Y[0,row]) != 0: # Only plot if non-zero
                if style == 'scatter':
                    ax[0].scatter(frequencies, 20*np.log10(np.abs(Y[:,row])) if not linear_mag else np.abs(Y[:,row]),linewidths=1.0,label=r'$'+labels[row]+r'$' if legend is not None else '_nolegend_')
                    ax[1].scatter(frequencies, np.angle(Y[:,row],deg=True),linewidths=1.0,label=r'$'+labels[row]+r'$' if legend is not None else '_nolegend_')
                else:
                    ax[0].plot(frequencies, 20*np.log10(np.abs(Y[:,row])) if not linear_mag else np.abs(Y[:,row]),linewidth=2.0,label=r'$'+labels[row]+r'$' if legend is not None else '_nolegend_')
                    ax[1].plot(frequencies, np.angle(Y[:,row],deg=True),linewidth=2.0,label=r'$'+labels[row]+r'$' if legend is not None else '_nolegend_')
    else: # SISO
        if style == 'scatter':
            ax[0].scatter(frequencies, 20*np.log10(np.abs(Y)),linewidths=1.0,label=r'$'+labels[0]+r'$' if legend is not None else '_nolegend_')
            ax[1].scatter(frequencies, np.angle(Y,deg=True),linewidths=1.0,label=r'$'+labels[0]+r'$' if legend is not None else '_nolegend_')
        else:
            ax[0].plot(frequencies, 20*np.log10(np.abs(Y)),linewidth=2.0,label=r'$'+labels[0]+r'$' if legend is not None else '_nolegend_')
            ax[1].plot(frequencies, np.angle(Y,deg=True),linewidth=2.0,label=r'$'+labels[0]+r'$' if legend is not None else '_nolegend_')

    if frequencies[0]>0:
        ax[0].set_xscale("log")
        ax[1].set_xscale("log")
    ax[0].set_xlim([frequencies[0], frequencies[-1]])
    ax[0].minorticks_on()
    ax[0].grid(visible=True, which='major', color='k', linestyle='-', alpha=0.5, linewidth=0.5)
    ax[0].set_ylabel('Magnitude [dB]') if not linear_mag else ax[0].set_ylabel('Magnitude')
    if title is None:
        ax[0].set_title('Bode plot for ' + str(len(frequencies)) + ' frequencies')
    else:
        ax[0].set_title(title)
    if legend is not None and np.prod(Y.shape[1:]) < 6*6: # Limit the number of legend entries to avoid overcrowding the figure
        ax[0].legend(loc='lower left', ncol=int(np.ceil(np.sqrt(np.prod(Y.shape[1:])))) if Y.ndim >= 2 else 1, prop={'size': 6}) # 'upper right' ,fancybox=True, shadow=True,
    ax[1].set_ylim([-190, 190])
    ax[1].set_yticks([-180, -90, 0, 90, 180])
    ax[1].set_xlim([frequencies[0], frequencies[-1]])
    ax[1].minorticks_on()
    ax[1].grid(visible=True, which='major', color='k', linestyle='-', alpha=0.5, linewidth=0.5)
    ax[1].set_ylabel('Phase [°]')
    ax[1].set_xlabel('Frequency [Hz]')

    if results_folder is not None:
        fig.savefig(results_folder + '\\' + file_name + ".pdf", format="pdf", bbox_inches="tight")
        if save_pickle:
            with open(results_folder + '\\' + file_name + ".pickle", 'wb') as f: pickle.dump(fig, f)
    if show_plot: plt.show() 
    if not return_plot: plt.close(fig)

    if save_data and results_folder is not None:
        filename = results_folder + '\\' + file_name + '.txt'
        header = ["f"]
        for label in legend: header.append(label)
        np.savetxt(filename, np.column_stack((frequencies, Y if Y.ndim < 3 else Y.reshape(Y.shape[:-2] + (-1,), order='C'))), delimiter='\t', header='\t'.join(header), comments='')
    
    if return_plot: return fig, ax


bode_plot.__doc__ = """
Draws a Bode plot of the complex vector or array Y: to each frequency in the list 'frequencies',
the entries in 'Y' are plotted in terms of their magnitude (in dB) and phase (in degrees).

Required
    Y               Complex vector or array to be plotted.
    frequencies     List of frequencies at which Y is computed [Hz]
    results_folder  Absolute path where the plot is to be stored.

Optional
    file_name       Name of the plot files. Default = 'Bode_plot'.
    style           Plot style, either 'scatter' or 'line'. Default = 'scatter'.
    show_plot       If True, the plot is shown interactively. Default = False. 
    save_data       If True, the data is saved in a text file. Default = False.
    legend          List of labels corresponding to the variables in Y. If set to None, no legend is shown. The default = [] results in labels corresponding to the indices of the variables.
    fig_handle      Tuple with existing figure and axes objects to be used for plotting. Default = None.
    return_plot     If True, the figure and axes objects are returned; otherwise, the figure is closed. Default = False.
    title           Title of the plot. Default = None.
    save_pickle     If True, the figure is saved in a pickle file. Default = False
"""

spectrum_plot.__doc__ = """
Draws a spectrum plot of the real-valued signals provided in the 'signals' argument.

Required
    signals         2D array where each column is a signal to be plotted.
    time            1D array with the time instants corresponding to the signals.
    results_folder  Absolute path where the plot is to be stored.
Optional
    file_name       Name of the plot files. Default = 'spectrum'.
    labels          List of labels for each signal. Default = None.
    style           Plot style, either 'scatter' or 'line'. Default = 'line'.
    show_plot       If True, the plot is shown interactively. Default = False. 
    save_data       If True, the data is saved in a text file. Default = False.
    show_dc         If True, the DC component is plotted, otherwise a x-log scale is used. Default = False.
    save_pickle     If True, the figure is saved in a pickle file. Default = False.
 """