""" 
This script demonstrates how to use the Transfer Function (TF) scan function to perform a frequency sweep of a model in PSCAD.
The script sets up the parameters for the frequency sweep, including the time settings, perturbation magnitude, and frequency range.
It then calls the frequency_sweep_TF function to execute the analysis and store the results in a specified folder.
"""
from ztoolacdc.frequency_sweep import frequency_sweep_TF
from os import getcwd

""" -------------------- PSCAD PROJECT ---------------------- """
pscad_folder = getcwd()  # Absolute location of the PSCAD workspace
results_folder = getcwd() + r'\Results'  # Location of the folder to store the results (if it doesn't exit, it is created)
workspace_name = "TF_test"  # Name of the PPSCAD workspace
project_name = "Test_TF"  # Name of the project
fortran_ext='.gf81' # Fortran compiler extension. Also '.gf46' can be used depending on the compiler

""" -------------------- SCAN SETTINGS ---------------------- """
num_parallel_sim = 8 # Number of parallel EMTDC simulations: min(EMTDC allowance by PSCAD license, number of computer's cores)
# The scan time is reduced proportionally to num_parallel_sim
multi_freq_scan = True # True = 8-tone sinusoidal perturbation, False = single-tone perturbation
# Number of frequencies to be scanned: for efficiency it should be a multiple of num_parallel_sim or of 8*num_parallel_sim when multi_freq_scan=True
f_points = 8*8*2  # Number of frequencies to be scanned
f_base = 0.5  # Base frequency in Hz (determines the frequency resolution)
f_min = 1.0  # Minimum frequency in Hz
f_max = 100.0  # Maximum frequency in Hz

start_fft = 1.0  # [s] Time for the system to reach steady-state after every injection
fft_periods = 1  # Number of periods used in the FFT for the lowest frequency
dt_injections = 0.0 # [s] Time after the decoupling to reach steady-state: it can be set to zero

t_snap = 0.2  # Time for the cold-start (snapshot) [s]
t_sim = start_fft + fft_periods / f_base  # Simulation time during the sinusoidal perturbation [s]
t_step = 5.0  # Simulation time step [us]
v_perturb_mag = 0.01 # In per unit w.r.t. the steady-state input

output_name = 'TF'  # Desired name for the output files

frequency_sweep_TF(t_snap=t_snap, t_sim=t_sim, t_step=t_step, v_perturb_mag=v_perturb_mag, num_parallel_sim=num_parallel_sim, dt_injections=dt_injections,
                   f_points=f_points, f_base=f_base, f_max=f_max, f_min=f_min, multi_freq_scan=multi_freq_scan, delete_PSCAD_output_files=True,
                   start_fft=start_fft, fft_periods=fft_periods, working_dir=pscad_folder, fortran_ext=fortran_ext, workspace_name=workspace_name,
                   project_name=project_name, results_folder=results_folder, output_files=output_name, plot_snapshot=True, plot_perturbation=1, target_blocks=['CL', 'OL'])