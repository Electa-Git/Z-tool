""" Simple script to test the Z-tool for a single-bus analysis """
from ztoolacdc import *
from os import getcwd
import numpy as np

""" -------------------- PSCAD PROJECT ---------------------- """
pscad_folder = getcwd()  # Absolute location of the PSCAD workspace
results_folder = getcwd() + r'\Results'  # Location of the folder to store the results (if it doesn't exit, it is created)
workspace_name = "Single_bus_example"  # Name of the PPSCAD workspace
project_name = "Simple_2L_VSC_RLC"  # Name of the project

""" -------------------- SCAN SETTINGS ---------------------- """
f_points = 8 * 50  # Number of frequencies to be scanned
f_base = 0.5  # Base frequency in Hz (determines the frequency resolution)
f_min = 1.0  # Minimum frequency in Hz
f_max = 500.0  # Maximum frequency in Hz

start_fft = 1.0  # [s] Time for the DUT to reach steady-state after every injection
fft_periods = 1  # Number of periods used in the FFT for the lowest frequency
dt_injections = 1  # [s] Time after the decoupling to reach steady-state

t_snap = 10  # Time for the cold-start (snapshot) [s]
t_sim = start_fft + fft_periods / f_base  # Simulation time during the sinusoidal perturbation [s]
t_step = 20.0  # Simulation time step [us]
v_perturb_mag = 0.02  # In per unit w.r.t. the steady-state voltage at each bus

output_files = 'single_bus_example'  # Desired name for the output files

""" -------------------- Frequency scan ---------------------- """
freq = create_freq.loglist(f_min=f_min, f_max=f_max, f_points=f_points, f_base=f_base)

frequency_sweep.frequency_sweep(t_snap=t_snap, t_sim=t_sim, t_step=t_step, dt_injections=dt_injections, f_base=f_base,
                                freq=freq, start_fft=start_fft, fft_periods=fft_periods, v_perturb_mag=v_perturb_mag,
                                working_dir=pscad_folder, workspace_name=workspace_name, project_name=project_name,
                                results_folder=results_folder, output_files=output_files, show_powerflow=True,
                                component_parameters=[["pq_tau_meas",0.001]])
# Note that you can change any Main canvas constants via the component_parameters argument for parametric changes

# Retreive admittances
Y_VSC = read_admittance.read_admittance(path=results_folder, involved_blocks="PCC-1", file_root=output_files)  # Side 1 of the PCC block is connected to the VSC
Y_grid = read_admittance.read_admittance(path=results_folder, involved_blocks="PCC-2", file_root=output_files)  # Side 2 of the PCC block is connected to the grid equivalent

""" -------------------- Stability analysis ---------------------- """
print("\nAnalysis of the PSCAD case")
# You can call stability.stability_analysis but this case is simple enough so we can do it step-by-step
L = np.matmul(np.linalg.inv(Y_grid.y), Y_VSC.y)  # Loop gain matrix
stability.nyquist(L, Y_VSC.f, results_folder=results_folder, filename="PSCAD_case")  # Application of the Generalized Nyquist Criterion
stability.EVD(G=np.linalg.inv(Y_grid.y + Y_VSC.y), frequencies=Y_VSC.f, results_folder=results_folder, filename="PSCAD_case") # Oscillation mode identification
stability.passivity(G=Y_VSC.y, frequencies=Y_VSC.f, results_folder=results_folder, filename="PSCAD_case",Yedge=Y_grid.y) # Passivity index of all subsystems
stability.small_gain(G1=np.linalg.inv(Y_grid.y),G2=Y_VSC.y, frequencies=Y_VSC.f, results_folder=results_folder, filename="PSCAD_case") # Evaluate the gain of all subsystems

""" -------------------- Stability analysis with different series compensation values ---------------------- """
print("\nAnalysis of different series compensation values using previously scanned converter admittance")
Wpu = np.array([[0,1],[-1,0]])  # Coupling matrix due to the abc-to-dq transformation T*dT^(-1)/dt
w0 = 2*np.pi*50  # Fundamental angular frequency [Hz]
Z_RL = np.linalg.inv(Y_grid.y)  # Grid-side impedance (RL in this case)
X_g = np.real(Z_RL[1,0,1])  # Extract the grid-side fundamental frequency reactance from the scanned data
print(" The grid-side inductance is",X_g/(2*np.pi*50),"H")
comp_level = np.arange(0.05,0.7,0.01)  # Compensation level from 5% to 70% of the grid inductance
Y_C = np.empty((len(Y_grid.f), 2, 2), dtype='cdouble')  # Initialize the capacitance admittance matrix
stability_assessment = []  # Stability analysis results
# Iterate the compensation levels: (1) compute the capacitance, (2) compute the new grid-side admittance, (3) check the stability
for case in range(len(comp_level)):
    C_g = 1/(w0*comp_level[case]*X_g)  # Xc = 1/(w*C) = 5-70% X_g
    for f_point, f in enumerate(Y_grid.f):
        Y_C[f_point,...] = 1j*2*np.pi*f*C_g*np.identity(2) + w0*C_g*Wpu  # dq-frame admittance matrix of a capacitor in SI
    Z_RLC = np.linalg.inv(Y_C) + Z_RL  # Series-capacitor compensated line
    # Evaluate the stability via both the GNC using eigenvalue decomposition and also the determinant (equivalent)
    stable = stability.nyquist(np.matmul(Z_RLC, Y_VSC.y), Y_VSC.f, results_folder=results_folder+r"\Compensation", filename="PSCAD_case_"+str(case), verbose=False)
    stability.nyquist_det(L=np.matmul(Z_RLC, Y_VSC.y), frequencies=Y_VSC.f, results_folder=results_folder+r"\Compensation", filename="PSCAD_case_"+str(case))
    stability_assessment.append(stable)
    if not stable:
        print(round(comp_level[case]*100,1),"% series compensation is small-signal unstable")
        stability.EVD(G=np.linalg.inv(np.linalg.inv(Z_RLC) + Y_VSC.y), frequencies=Y_VSC.f, results_folder=results_folder+r"\Compensation", filename="PSCAD_case_"+str(case),verbose=False)
        stability.nyquist_det(L=np.matmul(Z_RLC, Y_VSC.y), frequencies=Y_VSC.f, results_folder=results_folder+r"\Compensation", filename="PSCAD_case_"+str(case))

print("\nCompensation larger than",round(comp_level[sum(stability_assessment)]*100,1),"% might result in instability")
