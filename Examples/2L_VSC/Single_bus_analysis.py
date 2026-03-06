""" Simple script to illustrate the Z-tool for a single-bus analysis
The script performs the frequency scan of a converter (2L-VSC) connected to an RLC Thevenin equivalent in PSCAD.
It then retrieves the admittances of both the converter and the grid equivalent and performs a closed-loop 
stability analysis via the Generalized Nyquist Criterion (GNC).
Finally, it performs a parametric stability assessment for different series compensation levels
by re-using the scanned converter admittance and updated the grid-side impedance with additional compensation.
This Subsynchronous Control Interactions (SSCI) screening then can be validated against time-domain simulations.
"""
from ztoolacdc import *
from os import getcwd
import numpy as np
import time

""" -------------------- PSCAD PROJECT ---------------------- """
pscad_folder = getcwd()  # Absolute location of the PSCAD workspace
results_folder = getcwd() + r'\Results'  # Location of the folder to store the results (if it doesn't exit, it is created)
workspace_name = "Single_bus_example"  # Name of the PPSCAD workspace
project_name = "Simple_2L_VSC_RLC"  # Name of the project
fortran_ext='.gf81' # Fortran compiler extension. Also '.gf46' can be used depending on the compiler

""" -------------------- SCAN SETTINGS ---------------------- """
num_parallel_sim = 8 # Number of parallel EMTDC simulations: min(EMTDC allowance by PSCAD license, number of computer's cores)
# The scan time is reduced proportionally to num_parallel_sim
multi_freq_scan = True # True = 8-tone sinusoidal perturbation, False = single-tone perturbation
# Number of frequencies to be scanned: for efficiency it should be a multiple of num_parallel_sim or of 8*num_parallel_sim when multi_freq_scan=True
f_points = 8*8*6  # Number of frequencies to be scanned
f_base = 1.0  # Base frequency in Hz (determines the frequency resolution)
f_min = 1.0  # Minimum frequency in Hz
f_max = 500.0  # Maximum frequency in Hz

start_fft = 1.0  # [s] Time for the DUT to reach steady-state after every injection
fft_periods = 1  # Number of periods used in the FFT for the lowest frequency
dt_injections = 1  # [s] Time after the decoupling to reach steady-state

t_snap = 1.0  # Time for the cold-start (snapshot) [s]
t_sim = start_fft + fft_periods / f_base  # Simulation time during the sinusoidal perturbation [s]
t_step = 20.0  # Simulation time step [us]
v_perturb_mag = 0.01/2.83 if multi_freq_scan else 0.01 # In per unit w.r.t. the steady-state voltage at each bus

output_files = 'PSCAD_case'  # Desired name for the output files

""" -------------------- Frequency scan ---------------------- """
freq = create_freq.loglist(f_min=f_min, f_max=f_max, f_points=f_points, f_base=f_base)

# frequency_sweep.frequency_sweep(t_snap=t_snap, t_sim=t_sim, t_step=t_step, dt_injections=dt_injections, f_base=f_base,
#                                 freq=freq, start_fft=start_fft, fft_periods=fft_periods, v_perturb_mag=v_perturb_mag,
#                                 working_dir=pscad_folder, workspace_name=workspace_name, project_name=project_name,
#                                 results_folder=results_folder, output_files=output_files, show_powerflow=True, 
#                                 fortran_ext=fortran_ext, multi_freq_scan=multi_freq_scan)
# Note that you can change any Main canvas constants via the component_parameters argument for parametric changes
# E.g. the whole analysis/script can be run with ["pq_tau_meas",0.002] above for a different measurement filter time constant = different dynamics

# Retreive admittances
Y_VSC = read_admittance.read_admittance(path=results_folder, involved_blocks="PCC-1", file_root=output_files)  # Side 1 of the PCC block is connected to the VSC
Y_grid = read_admittance.read_admittance(path=results_folder, involved_blocks="PCC-2", file_root=output_files)  # Side 2 of the PCC block is connected to the grid equivalent

""" -------------------- Stability analysis ---------------------- """
print("\nAnalysis of the PSCAD case")
# stability.stability_analysis() evaluates the main small-signal analysis functions as commented below, but this case is simple enough to do it step-by-step
# stability.stability_analysis(results_folder=results_folder, file_root=output_files, indentations=[50.0], node_blocks=["PCC-1"]) # node_blocks specifies which blocks are considered as admittance in the stability analysis
L = np.matmul(np.linalg.inv(Y_grid.y), Y_VSC.y)  # Loop gain matrix
stability.nyquist(L, Y_VSC.f, results_folder=results_folder, filename=output_files)  # Application of the Generalized Nyquist Criterion
stability.EVD(G= Y_grid.y + Y_VSC.y, frequencies=Y_VSC.f, results_folder=results_folder, filename=output_files, Z_closedloop=False, bus_names=["d axis","q axis"]) # Oscillation mode identification
stability.passivity(G=Y_VSC.y, frequencies=Y_VSC.f, results_folder=results_folder, filename=output_files,Yedge=Y_grid.y) # Passivity index of all subsystems
stability.small_gain(G1=np.linalg.inv(Y_grid.y),G2=Y_VSC.y, frequencies=Y_VSC.f, results_folder=results_folder, filename=output_files) # Evaluate the gain of all subsystems

""" -------------------- Stability analysis with different series compensation values ---------------------- """
print("\nAnalysis of different series compensation values using previously scanned converter admittance")
results_folder_SSCI = results_folder+r"\SSCI_screening"  # Folder to store the SSCI screening results
filename_SSCI = "SSCI_case_"  # Prefix for the filenames of the SSCI screening results
t0_SSCI_screening = time.time() # Initial time for the SSCI screening

Wpu = np.array([[0,1],[-1,0]])  # Coupling matrix due to the abc-to-dq transformation T*dT^(-1)/dt
f0 = 50 # Fundamental frequency [Hz]
w0 = 2*np.pi*f0  # Fundamental angular frequency [rad/s]

Z_RL = np.linalg.inv(Y_grid.y)  # Grid-side impedance (RL in this case)
X_g = np.real(Z_RL[1,0,1])  # Extract the grid-side fundamental frequency reactance from the scanned data
print(" The grid-side inductance is",round(X_g/w0,5),"H")

comp_level = np.arange(0.05,0.7,0.01)  # Compensation level from 5% to 70% of the grid inductance
Y_C = np.empty((len(Y_grid.f), 2, 2), dtype='cdouble')  # Initialize the capacitance admittance matrix
stability_assessment = []  # Stability analysis results

# Iterate the compensation levels: (1) compute the capacitance, (2) compute the new grid-side admittance, (3) check the stability
for case in range(len(comp_level)):
    C_g = 1/(w0*comp_level[case]*X_g)  # Xc = 1/(w*C) = 5-70% X_g
    for f_point, f in enumerate(Y_grid.f):
        Y_C[f_point,...] = 1j*2*np.pi*f*C_g*np.identity(2) + w0*C_g*Wpu  # dq-frame admittance matrix of a capacitor in SI
    Z_RLC = np.linalg.inv(Y_C) + Z_RL  # Series-capacitor compensated line
    # Evaluate the stability via the GNC using eigenvalue decomposition: do plot all cases to speed up the screening but save the results as text files
    stable = stability.nyquist(np.matmul(Z_RLC, Y_VSC.y), Y_VSC.f, results_folder=results_folder_SSCI, filename=filename_SSCI+str(case), verbose=False, indentations =[f0], make_plot=False)
    stability_assessment.append(stable)
    if sum(stability_assessment) == case:
        # Show the participation factors and plot the Nyquist curve only for the first unstable case
        print(" First unstable case found for a series compensation level of",round(comp_level[case]*100,0),"%")
        stability.nyquist(np.matmul(Z_RLC, Y_VSC.y), Y_VSC.f, results_folder=results_folder_SSCI, filename=filename_SSCI+str(case), verbose=False, indentations =[f0], make_plot=True)
        stability.EVD(G=np.linalg.inv(Z_RLC) + Y_VSC.y, frequencies=Y_VSC.f, results_folder=results_folder_SSCI, filename=filename_SSCI+str(case), verbose=True, Z_closedloop=False, bus_names=["d axis","q axis"])

print("\nCompensation larger than",round(comp_level[sum(stability_assessment)]*100,0),"% might result in small-signal instability")
print("SSCI screening completed in", round(time.time()-t0_SSCI_screening,1),"seconds\n")

""" -------------------- Example of the admittance conversion functions ---------------------- """
# Convert from the dq frame to the stationary alpha-beta frame and plot the results
Y_ab_VSC, freq_ab = frame_conversion.dq2ab(Y_old_frame=Y_VSC.y, frequencies=Y_VSC.f, results_folder=results_folder, file_name="PCC-1", interpolate=False)
Y_ab_grid, freq_ab = frame_conversion.dq2ab(Y_old_frame=Y_grid.y, frequencies=Y_grid.f, results_folder=results_folder, file_name="PCC-2", interpolate=False)
plot_utils.bode_plot(Y=Y_ab_VSC,  frequencies=freq_ab, results_folder=results_folder, file_name="Y_ab_PCC-1", legend=["alpha","beta"])
plot_utils.bode_plot(Y=Y_ab_grid, frequencies=freq_ab, results_folder=results_folder, file_name="Y_ab_PCC-2", legend=["alpha","beta"])

# Convert from the alpha-beta frame to the positive-negative sequence frame and plot the results
Y_pn_VSC = frame_conversion.ab2pn(Y_old_frame=Y_ab_VSC, frequencies=freq_ab, results_folder=results_folder, file_name="PCC-1")
Y_pn_grid = frame_conversion.ab2pn(Y_old_frame=Y_ab_grid, frequencies=freq_ab, results_folder=results_folder, file_name="PCC-2")
plot_utils.bode_plot(Y=Y_pn_VSC[freq_ab>0],  frequencies=freq_ab[freq_ab>0], results_folder=results_folder, file_name="Y_pn_PCC-1", legend=["p","n"], style="solid")
plot_utils.bode_plot(Y=Y_pn_grid[freq_ab>0], frequencies=freq_ab[freq_ab>0], results_folder=results_folder, file_name="Y_pn_PCC-2", legend=["p","n"], style="solid")