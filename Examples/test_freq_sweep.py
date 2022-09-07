#!python3.7
""" Simple script to test the impedance measurement tool """
from Source import frequency_sweep
from os import getcwd

pscad_folder = getcwd() + r'\\Scan'
results_folder = pscad_folder + r'\\Results'
snapshot_file = 'Snapshot'
output_files = 'RL'

perturbations = 8 * 6  # Number of frequencies
f_base = 1  # Hz
f_min = 1  # Hz
f_max = 2000  # Hz

start_fft = 0.5  # [s] Time for the DUT to reach the new steady-state (injections)
fft_periods = 1  # Number of periods used in the FFT for the lowest frequency

t_sim = start_fft + fft_periods / f_base  # [s]
t_step = 2  # [us]
t_snap = 1  # [s]
v_perturb_mag = 1  # [kV]

print('Results for a sample time of', str(sample_step), 'us')
frequency_sweep.frequency_sweep(t_snap=t_snap, t_sim=t_sim, t_step=t_step, v_perturb_mag=v_perturb_mag,
                                f_points=perturbations, f_base=f_base, f_max=f_max, f_min=f_min,
                                start_fft=start_fft, fft_periods=fft_periods,
                                compute_yz=True, results_folder=results_folder,
                                working_dir=pscad_folder, project_name='DUT', workspace_name='DUTscan',
                                num_parallel_sim=8, snapshot_file=snapshot_file, output_files=output_files)


