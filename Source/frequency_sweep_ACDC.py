""" Find the documentation for this function at the end of the script or by typing help(frequency_sweep) """

import time as t  # Relative time
from os import listdir, chdir, getcwd, path, makedirs
import numpy as np  # Numerical python functions
from mhi.pscad import launch  # PSCAD automation functions
from Source.tools import *

def frequency_sweep_ACDC(t_snap=None, t_sim=None, t_step=None, sample_step=None, v_perturb_mag=None,
                         freq=None, f_points=None, f_base=None, f_min=None, f_max=None, working_dir=None,
                         delay_inj=0.0, freq_text_file='frequencies.txt', snapshot_file='Snapshot', take_snapshot=True,
                         project_name='DUT', workspace_name='DUTscan', scanid=None,
                         fortran_ext=r'.gf46', num_parallel_sim=8, component_parameters=None,
                         results_folder=None, output_files='Perturbation', compute_yz=False, save_td=False,
                         fft_periods=1, start_fft=1):
    """ --- Input data handling --- """
    if (t_snap or t_sim or t_step or v_perturb_mag or ((f_points or f_base or f_max or f_min) and freq)) is None:
        print(
            'One or more required arguments are missing!! \n Check the function documentation by typing help('
            'frequency_sweep) \n')
        return
    # If the sample time is not provided, it is set to half of the minimum required value (multiple of step_time)
    if sample_step is None:
        if f_max is None: f_max = max(freq)
        sample_step = round(t_step * np.floor((1e6 * 0.5 * 0.5 / f_max + t_step / 2) / t_step), 3)  # [us]
    if working_dir is None:
        working_dir = getcwd() + '\\'  # Location of the PSCAD workspace
    else:
        chdir(working_dir.encode('ascii', 'backslashreplace'))  # ('unicode_escape'))  Location of the PSCAD workspace
        working_dir = getcwd() + '\\'
    print('\nRunning from ' + working_dir + '\n')

    # The snapshot time must be a multiple of the sampling time
    t_snap = round(np.ceil(t_snap / (sample_step * 1e-6)) * sample_step * 1e-6, 6)

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
    # XYlist_freq.parameters(File=freq_text_file, path='RELATIVE', npairs=f_points)
    # Retreive scan block and type
    if scanid is not None:
        ScanBlock = main.find_first(Name=scanid)  # Try to find it by name
        if ScanBlock is None: ScanBlock = main.component(scanid)  # If it is not found, use the component id
    else:
        ScanBlock = main.find_first("Z_tool:ACscan")  # This assumes a single scan block in the main canvas
        if ScanBlock is None: ScanBlock = main.find_first("Z_tool:DCscan")
        if ScanBlock is None: ScanBlock = main.find_first("Z_tool:DCscanPM")
    if "DC" in ScanBlock.defn_name[1]:
        scantype = "DC"
    else:
        scantype = "AC"

    ScanBlock.parameters(V_perturb_mag=v_perturb_mag, Tdelay_inj=delay_inj + t_snap, selector=0)

    if component_parameters is not None:
        Parameters = [main.find('master:const', 'Param1'), main.find('master:const', 'Param2'),
                      main.find('master:const', 'Param3'), main.find('master:const', 'Param4'),
                      main.find('master:const', 'Param5'), main.find('master:const', 'Param6')]
        for i in range(min(len(component_parameters), len(Parameters))): Parameters[i].parameters(
            Value=component_parameters[i])

    all_pgb = main.find_all("master:pgb")  # Find all output channel in the main page
    for pgb in all_pgb: pgb.disable()  # Disable the outputs
    # selected_pgb = [main.find("master:pgb", "V_pcc_dq"), main.find("master:pgb", "I_pcc_dq")]
    # for pgb in selected_pgb:
    #     pgb.enable()  # Enable the selected outputs
    #     pgb.parameters(enab=0) # Disable PSCAD ploting

    # Set simulation-specific parameters
    if 'Perturbation' not in pscad.simulation_sets():
        simset = pscad.create_simulation_set('Perturbation')
        simset.add_tasks(project_name)
    else:
        simset = pscad.simulation_set('Perturbation')
    simset_task = simset.tasks()[0]  # The task is extracted

    if take_snapshot:  # It performs the snapshots
        print(' Running snapshot')
        t1 = t.time()
        ScanBlock.parameters(selector=0)  # No injection for the steady state
        simset_task.parameters(volley=1, affinity_type='DISABLE_TRACING',
                               ammunition=1)  # affinity_type = 'DISABLE_TRACING' disables the plotting
        simset_task.overrides(duration=t_snap + t_sim, time_step=t_step, plot_step=sample_step, start_method=0,
                              timed_snapshots=1, snapshot_file=snapshot_file + '.snp', snap_time=t_snap,
                              save_channels_file=snapshot_file + '.out', save_channels=1)
        simset.run()
        print(' Snapshot done in', round((t.time() - t1), 2), 'seconds')
    else:  # It performs the unperturbed simulation starting from the given snapshot
        print(' Running steady-state simulation')
        t1 = t.time()
        ScanBlock.parameters(selector=0)  # No injection for the steady state
        simset_task.parameters(volley=1, affinity_type='DISABLE_TRACING',
                               ammunition=1)  # affinity_type = 'DISABLE_TRACING' disables the plotting
        simset_task.overrides(duration=t_sim, time_step=t_step, plot_step=sample_step, start_method=1,
                              timed_snapshots=0, startup_inputfile=snapshot_file + '.snp',
                              save_channels_file=snapshot_file + '.out', save_channels=1)
        simset.run()
        print(' Steady-state simulation completed in', round((t.time() - t1), 2), 'seconds')

    if save_td or compute_yz:
        wait4pscad(time=1, pscad=pscad)
        t1 = t.time()
        ss = read_and_save.single_s(original_folder=working_dir + project_name + fortran_ext,
                                    target_filename=simset_task.overrides()['save_channels_file'][:-4],
                                    new_folder=results_folder, save=save_td, output=compute_yz, scantype=scantype)
        if take_snapshot: ss = ss[find_nearest(ss[:, 0], t_snap):, :]  # Get rid of the transient results, i.e. t<t_snap
        print(' Unperturbed simulation results collected in', round((t.time() - t1), 2), 'seconds')
    if scantype is "AC":
        # AC-type bus scan
        # d-axis injection
        print('\n Running single frequency d-axis injection simulations')
        t1 = t.time()
        ScanBlock.parameters(selector=1)  # Injection d
        simset_task.parameters(volley=num_parallel_sim, affinity_type='DISABLE_TRACING', ammunition=f_points)
        simset_task.overrides(duration=t_sim, time_step=t_step, plot_step=sample_step, start_method=1,
                              timed_snapshots=0, startup_inputfile=snapshot_file + '.snp',
                              save_channels_file=output_files + '_d.out', save_channels=1)
        simset.run()
        print(' d-axis injection finished in', round((t.time() - t1), 2), 'seconds')
        if save_td or compute_yz:
            wait4pscad(time=1, pscad=pscad)
            t2 = t.time()
            d_axis = read_and_save.multiple_s(original_folder=working_dir + project_name + fortran_ext,
                                              scantype=scantype,
                                              target_filename=simset_task.overrides()['save_channels_file'][:-4],
                                              new_folder=results_folder, save=save_td, output=compute_yz,
                                              n_sim=f_points)
            print(' d-axis injection results collected in', round((t.time() - t2), 2), 'seconds')

        # q-axis injection
        print(' Running single frequency q-axis injection simulations')
        t1 = t.time()
        ScanBlock.parameters(selector=2)  # Injection q
        simset_task.parameters(volley=num_parallel_sim, affinity_type='DISABLE_TRACING', ammunition=f_points)
        simset_task.overrides(duration=t_sim, time_step=t_step, plot_step=sample_step, start_method=1,
                              timed_snapshots=0, startup_inputfile=snapshot_file + '.snp',
                              save_channels_file=output_files + '_q.out', save_channels=1)
        simset.run()
        print(' q-axis injection finished in', round((t.time() - t1), 2), 'seconds')
        if save_td or compute_yz:
            wait4pscad(time=1, pscad=pscad)
            t2 = t.time()
            q_axis = read_and_save.multiple_s(original_folder=working_dir + project_name + fortran_ext,
                                              scantype=scantype,
                                              target_filename=simset_task.overrides()['save_channels_file'][:-4],
                                              new_folder=results_folder, save=save_td, output=compute_yz,
                                              n_sim=f_points)
            print(' q-axis injection results collected in', round((t.time() - t2), 2), 'seconds')

        if compute_yz:
            t2 = t.time()
            yz_computation.admittance(f_base=f_base, frequencies=freq, fft_periods=fft_periods, start_fft=start_fft,
                                      ss=ss,
                                      vi1_td=d_axis[:, 1:], vi2_td=q_axis[:, 1:], td=d_axis[:, 0], scantype=scantype,
                                      results_folder=results_folder, results_name=output_files)
            print(' Admittance computation finished in ', round((t.time() - t2), 2), 'seconds')
        # project.save()  # Save the project changes
        print(' Quitting PSCAD')
        wait4pscad(time=1, pscad=pscad)
        pscad.quit()  # Quit PSCAD
        print('\nTotal execution time', round((t.time() - t0) / 60, 2), 'minutes\n')
    else:
        # DC-type bus scan
        # Sinusoidal injection
        print('\n Running single frequency sinusoidal injection simulations')
        t1 = t.time()
        ScanBlock.parameters(selector=1)  # Single injection
        simset_task.parameters(volley=num_parallel_sim, affinity_type='DISABLE_TRACING', ammunition=f_points)
        simset_task.overrides(duration=t_sim, time_step=t_step, plot_step=sample_step, start_method=1,
                              timed_snapshots=0, startup_inputfile=snapshot_file + '.snp',
                              save_channels_file=output_files + '_dc.out', save_channels=1)
        simset.run()
        print(' Injections finished in', round((t.time() - t1), 2), 'seconds')
        if save_td or compute_yz:
            wait4pscad(time=1, pscad=pscad)
            t2 = t.time()
            dc_axis = read_and_save.multiple_s(original_folder=working_dir + project_name + fortran_ext,
                                               scantype=scantype,
                                               target_filename=simset_task.overrides()['save_channels_file'][:-4],
                                               new_folder=results_folder, save=save_td, output=compute_yz,
                                               n_sim=f_points)
            print(' Injection results collected in', round((t.time() - t2), 2), 'seconds')
        if compute_yz:
            t2 = t.time()
            yz_computation.admittance(f_base=f_base, frequencies=freq, fft_periods=fft_periods, start_fft=start_fft,
                                      ss=ss, vi1_td=dc_axis[:, 1:], td=dc_axis[:, 0], scantype=scantype,
                                      results_folder=results_folder, results_name=output_files)
            print(' Admittance computation finished in ', round((t.time() - t2), 2), 'seconds')
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


frequency_sweep_ACDC.__doc__ = """
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
