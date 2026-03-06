"""
ISGT 2024 AC/DC Energy Hub system analysis example: Z-tool frequency sweep and stability analysis 
For more information on the system and the analysis check the following paper 10.1109/ISGTEUROPE62998.2024.10863484
"""
from ztoolacdc.frequency_sweep import frequency_sweep
from ztoolacdc.stability import stability_analysis
from os import getcwd, path, walk

""" -------------------- PSCAD PROJECT INFO ---------------------- """
pscad_folder = getcwd() + r"\PSCAD model"  # Absolute location of the PSCAD workspace
topology = getcwd() + r'\topology.txt'  # Topology file
workspace_name = "Test_case"  # Name of the PSCAD workspace
project_name="Energy_hub"  # Name of the project
fortran_ext= '.gf81' # Fortran compiler extension, e.g. '.gf46' or '.gf81'

""" -------------------- SCAN SETTINGS ---------------------- """
num_parallel_sim = 8 # Number of parallel EMTDC simulations: min(allowance by PSCAD license, number of computer's cores)
# The scan time is reduced proportionally to num_parallel_sim
multi_freq_scan = True # True = 8-tone sinusoidal perturbation, False = single-tone perturbations
perturbations = 36*8*2 # Number of frequencies: for efficiency set to a multiple of num_parallel_sim or 8*num_parallel_sim for multi-tone scans

f_base = 1  # Base frequency in Hz (determines the frequency resolution)
f_min = 1  # Minimum frequency in Hz
f_max = 1000  # Maximum frequency in Hz

t_step = 20.0  # Simulation time step [us]
t_snap = 5.0  # Time for the cold-start (snapshot) [s]
dt_injections = 0.0  # [s] Time after the decoupling (t_snap) to reach steady-state (usually zero)

start_fft = 1.0  # [s] Time for the subsystems to reach steady state after every perturbation
fft_periods = 1  # Number of periods used in the FFT for the lowest frequency
t_sim = start_fft + fft_periods / f_base  # Simulation time during the sinusoidal perturbation [s]

v_perturb_mag = 0.01/2.83 if multi_freq_scan else 0.01  # In per unit w.r.t. the steady-state voltage at each bus

output_files = 'ISGT2024_stable'  # Desired name for the output files
results_folder = getcwd() + r'\Results'  # Location of the folder to store the results: if it doesn't exist, it is created

""" -------------------- RUN THE SCAN ---------------------- """
print("\n Case 1: stable system with base parameters \n")
frequency_sweep(t_snap=t_snap, t_sim=t_sim, t_step=t_step, dt_injections=dt_injections, multi_freq_scan=multi_freq_scan, num_parallel_sim=num_parallel_sim,
                f_points=perturbations, f_base=f_base, f_max=f_max, f_min=f_min, f_exclude = [50.0], start_fft=start_fft, fft_periods=fft_periods,v_perturb_mag=v_perturb_mag,
                snapshot_file="Snapshot", take_snapshot=True, fortran_ext=fortran_ext, release_certificates=False, 
                working_dir=pscad_folder, workspace_name=workspace_name, project_name=project_name, results_folder=results_folder, output_files=output_files,
                topology=topology, edge_dq_sym=True, show_powerflow=True, visualize_network=True, component_parameters=[["Control_switch", 1]])

## component_parameters=[["Name1", Value1], ["Name1", Value1], ...] is a list of ["Name",Value] where "Name" is the name of a constant in PSCAD to be set to Value before the scan

""" -------------------- RUN THE STABILITY ANALYSIS ---------------------- """
stability_analysis(topology=topology, results_folder=results_folder, file_root=output_files, indentations = [50.0], run_nyquist_det=False) # Nyquist contour indentation around the fundamental frequency

""" -------------------- CASE 2: UPDATED PARAMETERS: NEW SCAN & ANALYSIS ---------------------- """
print("\n Case 2: updated control parameters at MMC 2\n")
output_files_case2 = 'ISGT2024_unstable'  # Desired name for the output files
# Only control parameters are changed by setting "Control_switch" to 0 in PSCAD: no need to scan the passive grids -> scan_multi_ports=False
# The snapshot can be re-used as well if the operating point is the same: take_snapshot = False and snapshot_file="Snapshot" so it re-uses the existing snapshot
frequency_sweep(t_snap=t_snap, t_sim=t_sim, t_step=t_step, dt_injections=dt_injections, multi_freq_scan=multi_freq_scan, num_parallel_sim=num_parallel_sim,
                f_points=perturbations, f_base=f_base, f_max=f_max, f_min=f_min, f_exclude = [50.0], start_fft=start_fft, fft_periods=fft_periods,v_perturb_mag=v_perturb_mag,
                snapshot_file="Snapshot", take_snapshot=False, fortran_ext=fortran_ext, release_certificates=True, delete_PSCAD_output_files=True,
                working_dir=pscad_folder, workspace_name=workspace_name, project_name=project_name, results_folder=results_folder, output_files=output_files_case2,
                topology=topology, scan_multi_ports=False, component_parameters=[["Control_switch", 0]])

""" -------------------- RUN THE STABILITY ANALYSIS ---------------------- """
# The passive network scan can be re-used from case 1: copy & paste the files from "results_folder" with the root "output_files"
print("\nCopying the passive network files from case 1 to case 2...")
for root, _, files in walk(results_folder):
    for name in files:
        new_file_path = path.join(results_folder,name.replace(output_files,output_files_case2)) # For all files in the folder: copy and paste them renamed
        if name.endswith(".txt") and output_files+"#" in name and not path.exists(new_file_path): # Only IF the new file name does not exist already
            with open(path.join(root,name),"r") as file_source, open(new_file_path,"w") as file_copy:
                file_copy.write(file_source.read())
                print(" "+str(name),"copied")
print(" All passive network files copied \n")

# Now run the stability analysis using the scans with the updated controls
stability_analysis(topology=topology, results_folder=results_folder, file_root=output_files_case2, indentations = [50.0], run_nyquist_det=False) # Nyquist contour indentation around the fundamental frequency

print("\n ALL CASES COMPLETED! \n")
