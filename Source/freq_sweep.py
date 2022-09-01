#!python3.7
""" Find the documentation for this function at the end of the script or by typing help(frequency_sweep) """

import time as t  # Relative time
from os import getcwd, chdir, makedirs, path, listdir  # Operating system features
import matplotlib.pyplot as plt
import numpy as np  # Numerical python functions
from matplotlib import rcParams  # Text's parameters for plots
from mhi.pscad import launch  # PSCAD automation functions

rcParams['mathtext.fontset'] = 'cm'  # Font selection
rcParams['font.family'] = 'STIXGeneral'  # 'cmu serif'

def frequency_sweep(t_snap=None, t_sim=None, t_step=None, sample_step=None, v_perturb_mag=None,
                    freq=None, f_points=None, f_base=None, f_min=None, f_max=None, working_dir=None,
                    delay_inj=0.001, freq_text_file='frequencies.txt', snapshot_file='Snapshot', take_snapshot=True,
                    project_name='DUT', workspace_name='DUTscan',
                    fortran_ext=r'.gf46', num_parallel_sim=8, component_parameters=None,
                    results_folder=None, output_files='Perturbation', compute_yz=False, save_td=False,
                    fft_periods=1, start_fft=1):
    """ --- Input data handling --- """
    if (t_snap or t_sim or t_step or v_perturb_mag or ((f_points or f_base or f_max or f_min) and freq)) is None:
        print(
            'One or more required arguments are missing!! \n Check the function documentation by typing help('
            'frequency_sweep) \n')
        return
    if sample_step is None: sample_step = t_step  # If the sample time is not specified, it is set equal to the simulation time step
    if working_dir is None:
        working_dir = getcwd() + '\\'  # Location of the PSCAD workspace
    else:
        chdir(working_dir.encode('unicode_escape'))  # Location of the PSCAD workspace
        working_dir = getcwd() + '\\'
    print('\nRunning from ' + working_dir + '\n')
    if freq is None:
        freq = np.logspace(np.log10(f_min), np.log10(f_max), num=int(
            f_points))  # Create the frequency perturbation vector, float by deafult, use dtype='int16' for integers
        for j in range(freq.size): freq[j] = freq[j] - (
                    freq[j] % f_base)  # Modify the list so the values are multiples of a base frequency
        freq = np.unique(
            freq)  # Since it is a brute force calculation, there are repeated values: let's get rid of them
        multiples = np.arange(f_min, f_max + f_base, f_base)
        if len(freq) < int(f_points) and len(multiples) > len(freq):  # Not enough values
            scope = np.setxor1d(freq, multiples)  # XOR operator
            to_be_added = min(int(f_points) - len(freq), len(multiples) - len(freq))
            idx = np.floor(len(scope) / to_be_added) * np.arange(to_be_added)  # Add the values indexed uniformly
            for i in idx: freq = np.append(freq, scope[int(i)])
            freq.sort()

    # Frequency vector to XY table
    f_points = len(freq)
    with open(freq_text_file, 'w') as f:  # Create the .txt file
        f.write('! This file stores the frequency swept values \n')  # Write header
        for j in range(f_points): f.write(str(j + 1) + '\t' + str(freq[j]) + '\n')  # Write values
        f.write('ENDFILE:')  # Write end of the file
    f.close()

    if (results_folder is not None) and (not path.exists(results_folder)): makedirs(
        results_folder)  # Create the folder if it does not exist
    if (results_folder is None) and (save_td or compute_yz): results_folder = working_dir

    """ --- Main program --- """
    t0 = t.time()  # Initial time
    print('Launching PSCAD')
    pscad = launch(minimize=True)
    wait4pscad(time=1, pscad=pscad)
    t.sleep(10)  # Wait a bit more just in case PSCAD is still loading some stuff
    pscad.load(working_dir + workspace_name + ".pswx")  # Load workspace
    project = pscad.project(project_name)  # Where the model is defined
    main = project.canvas('Main')  # Main of the project
    project.focus()  # Switch PSCAD’s focus to the project

    print(' Setting parameters and configuration')
    # Set global parameters for all simulations
    selector = main.find('master:consti', 'selector')
    selector.parameters(Value=0)
    XYlist_freq = main.find('XYlist_freq')
    XYlist_freq.parameters(File=freq_text_file, path='RELATIVE', npairs=f_points)
    Tdelay_inj = main.find('Tdelay_inj')
    Tdelay_inj.parameters(X=delay_inj + t_snap)
    V_p_mag = main.find('V_perturb_mag')
    V_p_mag.parameters(Value=v_perturb_mag)
    if component_parameters is not None:
        Parameters = [main.find('master:const', 'Param1'), main.find('master:const', 'Param2'),
                      main.find('master:const', 'Param3'), main.find('master:const', 'Param4'),
                      main.find('master:const', 'Param5'), main.find('master:const', 'Param6')]
        for i in range(min(len(component_parameters), len(Parameters))): Parameters[i].parameters(
            Value=component_parameters[i])

    # Set simulation-specific parameters
    simset = pscad.simulation_set('Perturbation')
    simset_task = simset.tasks()[0]  # The task is extracted

    if take_snapshot:  # It performs the snapshots
        print(' Running snapshot')
        t1 = t.time()
        selector.parameters(Value=0)  # No injection for the snapshot
        simset_task.parameters(volley=1, affinity_type='DISABLE_TRACING',
                               ammunition=1)  # affinity_type = 'DISABLE_TRACING' disables the plotting
        simset_task.overrides(duration=t_snap + t_sim, time_step=t_step, plot_step=sample_step, start_method=0,
                              timed_snapshots=1, snapshot_file=snapshot_file + '.snp', snap_time=t_snap,
                              save_channels_file=snapshot_file + '.out', save_channels=1)
        simset.run()
        print(' Snapshot done in', round((t.time() - t1), 2), 'seconds')

    if results_folder is not None and (save_td or compute_yz):
        wait4pscad(time=1, pscad=pscad)
        t1 = t.time()
        ss = read_and_save_ss(original_folder=working_dir + project_name + fortran_ext,
                              target_filename=simset_task.overrides()['save_channels_file'][:-4],
                              new_folder=results_folder, save=save_td,
                              output=compute_yz)
        print(' Snapshot results collected in', round((t.time() - t1), 2), 'seconds')

    # d-axis injection
    print('\n Running single frequency d-axis injection simulation')
    t1 = t.time()
    selector.parameters(Value=1)
    simset_task.parameters(volley=num_parallel_sim, affinity_type='DISABLE_TRACING', ammunition=f_points)
    simset_task.overrides(duration=t_sim, time_step=t_step, plot_step=sample_step, start_method=1,
                          timed_snapshots=0, startup_inputfile=snapshot_file + '.snp',
                          save_channels_file=output_files + '_d.out', save_channels=1)
    simset.run()
    print(' d-axis injection finished in', round((t.time() - t1), 2), 'seconds')
    if results_folder is not None and (save_td or compute_yz):
        wait4pscad(time=1, pscad=pscad)
        t2 = t.time()
        d_axis = read_and_save_ms(original_folder=working_dir + project_name + fortran_ext,
                                  target_filename=simset_task.overrides()['save_channels_file'][:-4],
                                  new_folder=results_folder, save=save_td,output=compute_yz, n_sim=f_points)
        print(' d-axis injection results collected in', round((t.time() - t2), 2), 'seconds')

    # q-axis injection
    print(' Running single frequency q-axis injection simulation')
    t1 = t.time()
    selector.parameters(Value=2)
    simset_task.parameters(volley=num_parallel_sim, affinity_type='DISABLE_TRACING', ammunition=f_points)
    simset_task.overrides(duration=t_sim, time_step=t_step, plot_step=sample_step, start_method=1,
                          timed_snapshots=0, startup_inputfile=snapshot_file + '.snp',
                          save_channels_file=output_files + '_q.out', save_channels=1)
    simset.run()
    print(' q-axis injection finished in', round((t.time() - t1), 2), 'seconds')
    if results_folder is not None and (save_td or compute_yz):
        wait4pscad(time=1, pscad=pscad)
        t2 = t.time()
        q_axis = read_and_save_ms(original_folder=working_dir + project_name + fortran_ext,
                                  target_filename=simset_task.overrides()['save_channels_file'][:-4],
                                  new_folder=results_folder, save=save_td, output=compute_yz, n_sim=f_points)
        print(' q-axis injection results collected in', round((t.time() - t2), 2), 'seconds')

    if compute_yz:
        t2 = t.time()
        yz_computation(f_base=f_base, frequencies=freq, fft_periods=fft_periods, start_fft=start_fft,
                       ss=ss[find_nearest(ss[:, 0], t_snap):, :],
                       vi1_td=d_axis[:, 1:], vi2_td=q_axis[:, 1:], td=d_axis[:, 0], results_folder=results_folder,
                       results_name=output_files)
        print(' Admittance computation finished in ', round((t.time() - t2), 2), 'seconds')

    project.save()  # Save the project changes
    print(' Quitting PSCAD')
    wait4pscad(time=1, pscad=pscad)
    pscad.quit()  # Quit PSCAD
    print('\nTotal execution time', round((t.time() - t0) / 60, 2), 'minutes\n')


def read_and_save_ss(original_folder=None, target_filename=None, new_folder=None, output=False, save=False):
    n = 4  # Number of simulation variables to be retreived (excluding time)
    values = np.loadtxt(original_folder + '\\' + target_filename + "_01.out", skiprows=1, usecols=np.arange(0, n + 1))
    if save: np.savetxt(new_folder + '\\' + target_filename + '.txt', values, delimiter='\t', comments='')
    if output: return values


def read_and_save_ms(n_sim=None, original_folder=None, target_filename=None, new_folder=None, output=False, save=False):
    # Filter the file names to identify the target multiple output files
    files = [file for file in listdir(original_folder) if
             (file.endswith(".out") and (file.count(target_filename) > 0))]  # end in .out and contain name
    # More file filtering: file name followed by an _ and another two _ near the end
    files_filtered = [file for file in files if (
                file[len(target_filename)] == '_' and file.count("_", len(target_filename)) == 2 and len(file) > len(target_filename))]
    if len(files_filtered) != n_sim: files_filtered = [file for file in files_filtered if file.endswith(
        "_01.out")]  # If multiple output vars, it only reads the 1st file
    n = 4  # Number of simulation variables to be retreived (excluding time)
    first_values = np.loadtxt(original_folder + '\\' + files_filtered[0], skiprows=1,
                              usecols=np.arange(0, n + 1))  # Only for the first file, it loads the first column (time)
    values = np.empty((first_values.shape[0], n * len(files_filtered) + 1),
                      dtype='d')  # Preallocation of all the memory space
    values[:, 0:n + 1] = first_values  # Save the first file
    del first_values, files
    j = n + 1  # File counter
    for file in files_filtered[1:n_sim]:
        values[:, j:j + n] = np.loadtxt(original_folder + '\\' + file, skiprows=1,
                                        usecols=np.arange(1, n + 1))  # Save the rest without the time / first column
        j += n
    if save: np.savetxt(new_folder + '\\' + target_filename + '.txt', values, delimiter='\t', comments='')
    if output: return values


def yz_computation(f_base=None, frequencies=None, fft_periods=1, start_fft=None,
                   ss=None, vi1_td=None, vi2_td=None, td=None, results_folder=None, results_name='Y'):
    dt = np.mean([td[i + 1] - td[i] for i in range(min(len(td), 100))])  # Sampling time [s]

    # Subtracts the steady state values from the signals after the system reached steady state
    start_idx = find_nearest(td, start_fft)
    if td[0] == dt: ss = ss[1:, :]  # If the simulations do not start at 0 but at dt, then shift the snapshot by dt
    ss_ext = np.tile(ss[start_idx:, 1:], (1, int(round(vi1_td.shape[1] / 4))))  # Extend the steady-state matrix by the number of simulations
    deltavi1_td = vi1_td[start_idx:, :] - ss_ext  # Small-signal steady state computation
    deltavi2_td = vi2_td[start_idx:, :] - ss_ext

    Y_dd = np.empty((len(frequencies),), dtype='cdouble')  # also dtype='csingle'
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
        if not path.exists(results_folder): makedirs(results_folder)  # Create the folder if it does not exist
        ##                np.savetxt(results_folder+'\\'+'frequencies.txt', frequencies, header="f", delimiter='\t',comments='')
        ##                np.savetxt(results_folder+'\\'+'Y.txt', np.stack((Y_dd,Y_dq,Y_qd,Y_qq),axis=1), header="dd\tdq\tqd\tqq", delimiter='\t',comments='')
        np.savetxt(results_folder + '\\' + results_name + '.txt',
                   np.stack((frequencies, Y_dd, Y_dq, Y_qd, Y_qq), axis=1), header="f\tdd\tdq\tqd\tqq", delimiter='\t',
                   comments='')
        ##                np.savetxt(results_folder+'\\'+'Y_dd.txt', Y_dd,  delimiter='\t')
        ##                np.savetxt(results_folder+'\\'+'Y_dq.txt', Y_dq,  delimiter='\t')
        ##                np.savetxt(results_folder+'\\'+'Y_qd.txt', Y_qd,  delimiter='\t')
        ##                np.savetxt(results_folder+'\\'+'Y_qq.txt', Y_qq,  delimiter='\t')

        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
        ax[0].scatter(frequencies, 20 * np.log10(np.abs(Y_dd)), marker='o', facecolors='none', edgecolors='b',
                      linewidths=1.5, label=r'$Y_{dd}$')
        ax[0].scatter(frequencies, 20 * np.log10(np.abs(Y_dq)), marker='x', c='r', linewidths=1.5, label=r'$Y_{dq}$')
        ax[0].scatter(frequencies, 20 * np.log10(np.abs(Y_qd)), marker='+', c='m', linewidths=1.5, label=r'$Y_{qd}$')
        ax[0].scatter(frequencies, 20 * np.log10(np.abs(Y_qq)), marker='.', c='g', linewidths=1.5, label=r'$Y_{qq}$')
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
        fig.savefig(results_folder + '\\' + "Admittance.pdf", format="pdf", bbox_inches="tight")


def wait4pscad(time=1, pscad=None):
    busy = pscad.is_busy()
    while busy:
        t.sleep(time)  # Wait a bit
        busy = pscad.is_busy()


def find_nearest(array, value):  # Very efficient function to find the nearest value to a given one and its position
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or np.abs(value - array[idx - 1]) < np.abs(value - array[idx])):
        return idx - 1
    else:
        return idx


frequency_sweep.__doc__ = """
Author: Francisco Javier Cifuentes García
V0.1 [03/08/2022]
PSCAD automation
The values of the frequency list are multiples of the base frequency.
The simulation configuration is set to perform one snapshot and then the frequency sweeps.
The function accepts several input arguments to customize the frequency sweep:

Required
        t_snap			Time when the snapshot is taken [s].
        take_snapshot		Bool: Does the user want to take a snapshot? Default = False. A previous snapshot can still be used for the Y computation.
                                The snapshot simulation runs for t_snap + t_sim so as to save the steady-state unperturbed waveforms.
        t_sim			Duration of each frequency injection simulation [s].
        t_step			Simulation time step [us].
        sample_step		Sample time of the output channels [us].
        v_perturb_mag		Voltage perturbation peak value [kV].
        freq			Frequencies to perform the injections [Hz]. Alternatively, the user can provide info to compute the list.
      
Optional
        f_points			Number of frequency perturbation points.
        f_base		 	Base frequency [Hz].
        f_min			Start frequency [Hz].
        f_max			End frequency [Hz].
        fft_periods 		Number of periods used to compute the FFT. Default = 1.
        start_fft		Time for the DUT to reach the new steady-state (injections) [s] . 
	delay_inj 	 	Waiting time until the injections are performed [s]. Default = 0.001 s.
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

read_and_save_ms.__doc__ = """
Function that reads Multiple Simulations (ms) results and saves them into a dedicated file and/or loads them into memory for further processing.
The function accepts several input arguments to customize the reading:
Required
        n_sim		 Number of simulation results to be read: total number of simulations.
        original_folder	 Absolute path of the folder where the results can be found.
        target_filename	 Common name of the files to be read.       

Optional
        output		 Boolean variable to control if the read results should be returned as an array. output = True returns an array with the results.
        save		 Boolean variable to control the saving of the read data. False means the data is not saved.
        new_folder	 Absolute path of the destination folder in case save = True.
        output_filename	 Output file name: no extension needed. By default it saves the data into .txt files where the first column is the time and the rest
                         are the concatenated extracted variables from each simulation in the order that they where performed.

Furthermore, there is an additional internal parameter that indicates how many variables are read from each file. For impedance computation
only two currents and voltages are needed, thus n = 4 variables.

"""
read_and_save_ss.__doc__ = """
Function that reads Single Simulation (ss) results and saves them into a dedicated file and/or loads them into memory for further processing.
The function accepts several input arguments to customize the reading:
Required
        original_folder	 Absolute path of the folder where the results can be found.
        target_filename	 Common name of the file to be read.       

Optional
        output		 Boolean variable to control if the read results should be returned as an array. output = True returns an array with the results.
        save		 Boolean variable to control the saving of the read data. False means the data is not saved.
        new_folder	 Absolute path of the destination folder in case save = True.
        output_filename	 Output file name: no extension needed. By default it saves the data into .txt files where the first column is the time and the rest
                         are the extracted variables.

Furthermore, there is an additional internal parameter that indicates how many variables are read from each file. For impedance computation
only two currents and voltages are needed, thus n = 4 variables.

"""
