#!python3.7
""" Simple script to test the Z-tool for the energy hub """
from Source.frequency_sweep import frequency_sweep
from Source.tools.stability import stability_analysis
from os import getcwd, system

""" -------------------- PSCAD PROJECT ---------------------- """
pscad_folder = getcwd() + r"\PSCAD model"  # Absolute location of the PSCAD workspace
topology = getcwd() + r'\topology.txt'  # Topology file
workspace_name = "Test_case"  # Name of the PSCAD workspace
project_name="Energy_hub"  # Name of the project

""" -------------------- SCAN SETTINGS ---------------------- """
perturbations = 8 * 50   # Number of frequencies to be scanned
f_base = 1  # Base frequency in Hz (determines the frequency resolution)
f_min = 1  # Minimum frequency in Hz
f_max = 1000  # Maximum frequency in Hz

start_fft = 1.0  # [s] Time for the DUT to reach steady-state after every injection
fft_periods = 1  # Number of periods used in the FFT for the lowest frequency
dt_injections = 2  # [s] Time after the decoupling to reach steady-state

t_snap = 28   # Time for the cold-start (snapshot) [s]
t_sim = start_fft + fft_periods / f_base  # Simulation time during the sinusoidal perturbation [s]
t_step = 10.0  # Simulation time step [us]
v_perturb_mag = 0.02  # In per unit w.r.t. the steady-state voltage at each bus

output_files = 'ISGT_stable'  # Desired name for the output files
results_folder = getcwd() + r'\Results stable'  # Location of the folder to store the results (if it doesn't exit, it is created)

""" -------------------- RUN THE ANALYSIS ---------------------- """
print("\n Stable case \n")
frequency_sweep(t_snap=t_snap, t_sim=t_sim, t_step=t_step, dt_injections=dt_injections,
                f_points=perturbations, f_base=f_base, f_max=f_max, f_min=f_min, start_fft=start_fft, fft_periods=fft_periods,v_perturb_mag=v_perturb_mag,
                working_dir=pscad_folder, workspace_name=workspace_name, project_name=project_name, results_folder=results_folder, output_files=output_files,
                topology=topology, edge_dq_sym=True, scan_passives=True, show_powerflow=True, visualize_network=True)

stability_analysis(topology=topology, results_folder=results_folder, file_root=output_files)

print("\n Unstable case re-using previous snapshot \n")
output_files = 'ISGT_unstable'  # Desired name for the output files
results_folder = getcwd() + r'\Results unstable'  # Location of the folder to store the results (if it doesn't exit, it is created)
# component_parameters=[["Name1", Value1], ["Name1", Value1]] is a list of ["Name",Value] where "Name" is the name of a constant in PSCAD to be set to Value
frequency_sweep(t_snap=t_snap, t_sim=t_sim, t_step=t_step, dt_injections=dt_injections, snapshot_file="Snapshot_stable",
                f_points=perturbations, f_base=f_base, f_max=f_max, f_min=f_min, start_fft=start_fft, fft_periods=fft_periods,v_perturb_mag=v_perturb_mag,
                working_dir=pscad_folder, workspace_name=workspace_name, project_name=project_name, results_folder=results_folder, output_files=output_files,
                topology=topology, scan_passives=False, show_powerflow=True, visualize_network=True, component_parameters=[["Control_switch", 1]])

stability_analysis(topology=topology, results_folder=results_folder, file_root=output_files)


print("\n ALL CASES COMPLETED!! \n")
