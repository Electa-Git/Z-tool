#!python3.7
""" Simple script to test the impedance measurement tool """
from Source.freq_sweep import *

PSCAD_folder = r'C:\Users\fcifuent\Desktop\KU Leuven\Projects\Z based analysis\Freq measurement\PSCAD_MODEL'
results_folder = PSCAD_folder + r'\\Results'
snapshot_file = 'Snapshot'
output_files = 'Perturbation'

perturbations   = 8    # Number of frequencies
fbase           = 2     # Hz
fmin            = 2     # Hz
fmax            = 20    # Hz

delay_inj       = 0.001  # [s] Delay between the steady-state signal and the injections
start_FFT       = 4      # [s] Time for the DUT to reach the new steady-state (injections)
FFT_periods     = 1      # Number of periods used in the FFT for the lowest frequency

t_snap          = 3     # [s]
t_sim           = delay_inj + start_FFT + FFT_periods/fbase # [s]
t_step          = 1     # [us]
sample_step     = 25    # [us]
V_perturb_mag   = 0.5 # [kV]

#print('Results for a sample time of',str(sample_step),'us')
frequency_sweep(t_snap  = t_snap, t_sim = t_sim, t_step = t_step, sample_step = sample_step, V_perturb_mag = V_perturb_mag,
                fpoints = perturbations, fbase = fbase, fmax = fmax, fmin = fmin, delay_inj = delay_inj,
                working_dir = PSCAD_folder, snapshot_file = snapshot_file, output_files = output_files,
                compute_YZ = True, start_FFT = start_FFT, FFT_periods = 1, results_folder = results_folder) 
