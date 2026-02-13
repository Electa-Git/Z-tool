"""
Device admittance scanning at different operating points.
This script performs a parametric frequency sweep of a PSCAD model at different operating points by changing the powers and AC voltage magnitude.
The parameters are given in the Main canvas constants and changed via the component_parameters argument of the frequency_sweep function.
The frequency scan is performed for all combinations of the specified parameter values.
The scanned admittances flattened are stored in a single text file together with the corresponding operating point parameters.
Finally, an interpolation of the admittance at a different operating point is performed by using the scanned data and compared with a direct scan at that point.
"""
from ztoolacdc.frequency_sweep import frequency_sweep
from ztoolacdc.read_admittance import read_admittance
from ztoolacdc.plot_utils import bode_plot
from os import getcwd, path
from shutil import copyfile
import time as t  # Relative time
import numpy as np
from scipy.interpolate import LinearNDInterpolator

""" -------------------- PSCAD PROJECT ---------------------- """
pscad_folder = getcwd()  # Absolute location of the PSCAD workspace
results_folder = getcwd() + r'\Results'  # Location of the folder to store the results (if it doesn't exit, it is created)
workspace_name = "Parametric_study"  # Name of the PPSCAD workspace
project_name = "Parametric_sweep"  # Name of the project
fortran_ext='.gf81' # Fortran compiler extension. Also '.gf46' can be used depending on the compiler

""" -------------------- SCAN SETTINGS ---------------------- """
num_parallel_sim = 8 # Number of parallel EMTDC simulations: min(EMTDC allowance by PSCAD license, number of computer's cores)
# The scan time is reduced proportionally to num_parallel_sim
multi_freq_scan = True # True = 8-tone sinusoidal perturbation, False = single-tone perturbation
# Number of frequencies to be scanned: for efficiency it should be a multiple of num_parallel_sim or of 8*num_parallel_sim when multi_freq_scan=True
f_points = 8*8*6  # Number of frequencies to be scanned
f_base = 1.0  # Base frequency in Hz (determines the frequency resolution)
f_min = 1.0  # Minimum frequency in Hz
f_max = 1000.0  # Maximum frequency in Hz

start_fft = 0.5  # [s] Time for the DUT to reach steady-state after every injection: it depends on the device under test
fft_periods = 1  # Number of periods used in the FFT for the lowest frequency

t_snap = 0.5  # Time for the cold-start (snapshot) [s]
t_sim = start_fft + fft_periods / f_base  # Simulation time during the sinusoidal perturbation [s]
t_step = 20.0  # Simulation time step [us]
v_perturb_mag = 0.01/2.83 if multi_freq_scan else 0.01 # In per unit w.r.t. the steady-state voltage at each bus

output_name = 'IBR' # Desired name for the output files

""" -------------------- OPERATING POINTS ------------------- """
# Few cases just for testing: the total number of cases grows very fast with the number of values
# V_pu = [0.95, 1.0] # Voltage magnitude in per unit
# p_ref = [0.8, 0.9] # Active power reference in per unit
# q_ref = [0, 0.1] # Reactive power reference in per unit

# Many cases: capability curve sampling
V_pu = [0.95, 1.0, 1.05] # Voltage magnitude in per unit
p_ref = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] # Active power reference in per unit
q_ref = [-0.333, -0.166, 0, 0.166, 0.333] # Reactive power reference in per unit

parameter_lists = [V_pu, p_ref, q_ref]
parameter_names = ["V_pu", "p_ref", "q_ref"]
parameter_num   = len(parameter_lists)

# Generate all combinations of the parameters
grid = np.meshgrid(*parameter_lists, indexing='ij')
cases = np.stack(grid, axis=-1).reshape(-1, parameter_num)
number_of_cases = cases.shape[0]

admittance_size = 2  # Size of the admittance matrix
results_file = path.join(results_folder,output_name+".txt")
print('There are '+str(number_of_cases)+" cases involving", ", ".join(parameter_names), "parameters")

""" -------------------- PARAMETRIC FREQUENCY SCAN ---------------------- """
# The Main canvas constants corresponding to the "parameter_names" are changed iteratively via the component_parameters argument and the frequency scan is performed for each case
# The resulting admittances are stored in a single text file together with the case-specific parameter values
t0 = t.time()  # Initial time
# Create the file storing the results: name of the parameters and sliced admittance elements as FrequencyIndex_MatrixIndex
header = " ".join([parameter_names[i] for i in range(parameter_num)] + [str(f+1)+"_"+str(y+1) for f in range(f_points) for y in range(admittance_size*admittance_size)])
with open(results_file, "w") as f:
    f.write(header + "\n")
    for case in range(number_of_cases):
        t1 = t.time()
        output_files = output_name + "_" + str(case)
        component_parameters = [ [parameter_names[i], cases[case,i]] for i in range(parameter_num) ] # Generate the list of parameter name-value pairs
        
        frequency_sweep(t_snap=t_snap, t_sim=t_sim, t_step=t_step, f_base=f_base, num_parallel_sim=num_parallel_sim,
                        f_min=f_min, f_max=f_max, f_points=f_points, start_fft=start_fft, fft_periods=fft_periods, v_perturb_mag=v_perturb_mag,
                        working_dir=pscad_folder, workspace_name=workspace_name, project_name=project_name, results_folder=results_folder,
                        output_files=output_files, fortran_ext=fortran_ext, component_parameters=component_parameters, snapshot_file="Snapshot",
                        take_snapshot=True, show_powerflow=False, multi_freq_scan=multi_freq_scan, make_plot=False, delete_PSCAD_output_files=False,
                        release_certificates=False, launch_and_load_PSCAD= True if case==0 else False, close_PSCAD=True if case==number_of_cases-1 else False) # Only launch PSCAD for the first case and keep it open until the last case is done
        
        Y = read_admittance(path=results_folder, involved_blocks=["PCC-1"], file_root=output_files) # "involved_blocks" specifies the side of the scan block to retrieve the admittance: IBR is at side 1 of the PCC block
        Y_flat = Y.y.ravel(order='C') # Flatten the admittance
        # Combine the parameter values and the admittance into a single vector
        combined = np.concatenate([np.asarray([comp[1] for comp in component_parameters],dtype='double'), Y_flat])
        np.savetxt(f, combined[None, :], fmt="%.10e") # Add the flattened vector as a row
        
        print('Case '+str(case+1)+"/"+str(number_of_cases)+" completed in", round((t.time() - t1), 2), 'seconds')

# Copy the frequency file separately: it is the same for all cases so a separate file reduces the total storage usage
copyfile(pscad_folder+r'\Scan_options\frequencies.txt', results_folder+r'\frequencies.txt')
with open(results_folder+r'\frequencies.txt', 'r') as freqs:
    data = freqs.read().splitlines(True)
with open(results_folder+r'\frequencies.txt', 'w') as freqs:
    freqs.write("Index\tf\n")
    freqs.writelines(data[1:-1]) # Remove first and last lines

print("All scans completed in", round((t.time() - t0)/60, 2), 'minutes \n')

""" -------------------- INTERPOLATION AT A NEW OPERATING POINT ---------------------- """
# Load the scanned data
scanned_data = np.loadtxt(results_file, dtype='cdouble', skiprows=1)
frequencies = np.loadtxt(results_folder+r'\frequencies.txt', skiprows=1)
frequencies = frequencies[:,1] # Retrieve the frequency vector

# Create the interpolating function
interpolating_function = LinearNDInterpolator(np.real(scanned_data[:,0:parameter_num]), scanned_data[:,parameter_num:])

# New operating point:   V_pu, p_ref, q_ref
new_op_point = np.array([0.98, 0.85, 0.05]) # It can be any point in the convex hull of the scanned points
component_parameters = [ [parameter_names[i], new_op_point[i]] for i in range(parameter_num) ]

# Interpolate at the new operating point- and save the results in a text file
Y_interpolated_flat = interpolating_function(new_op_point)
N = int(np.sqrt(Y_interpolated_flat.size/len(frequencies)))  # The size of the admittance matrix, i.e. admittance_size
Y_interpolated = Y_interpolated_flat.reshape((len(frequencies), N, N), order="C")

# Plot and save the interpolated admittance
print("\nThe interpolated admittance at the new operating point has been computed and plotted.")

# Compute the actual admittance at the new op. point for comparison
output_files = output_name + "_direct_scan"
frequency_sweep(t_snap=t_snap, t_sim=t_sim, t_step=t_step, f_base=f_base, num_parallel_sim=num_parallel_sim,
                f_min=f_min, f_max=f_max, f_points=f_points, start_fft=start_fft, fft_periods=fft_periods, v_perturb_mag=v_perturb_mag,
                working_dir=pscad_folder, workspace_name=workspace_name, project_name=project_name, results_folder=results_folder,
                output_files=output_files, fortran_ext=fortran_ext, component_parameters=component_parameters, snapshot_file="Snapshot",
                take_snapshot=True, show_powerflow=False, multi_freq_scan=multi_freq_scan, delete_PSCAD_output_files=True)

"""-------------- Compare the actual and interpolated admittances -----------------"""
legend = ["d", "q"]
Y_scan = read_admittance(path=results_folder, involved_blocks="PCC-1", file_root=output_files)  # Side 1 of the PCC block is connected to the VSC
# Plot the scanned admittance at the new operating point
fig_comparison = bode_plot(Y_scan.y, Y_scan.f, return_plot=True, style="line", legend=False) 
# Add the interpolated admittance to the same plot for comparison
fig_comparison, ax_comparison = bode_plot(Y_interpolated, frequencies, return_plot=True, legend=legend, fig_handle=fig_comparison) 
ax_comparison[0].set_title('Comparison of scanned (solid lines) and interpolated (dots) admittances for ' + str(len(frequencies)) + ' frequencies')
fig_comparison.savefig(results_folder + '\\Scan_vs_interpolated.pdf', format="pdf", bbox_inches="tight") # Save the figure
