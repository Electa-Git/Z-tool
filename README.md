# Z-tool
Z-tool is a Python-based implementation for the stability analysis of modern power systems.
The core functionalities are impedance/admittance scan and stability assessment.
The analysis relies on an existing system model in the EMT simulation software [PSCAD]([url](https://www.pscad.com/)).

The following features are currently implemented:
- [x] Voltage perturbation-based admittance scan at several buses, including HVDC converters and black-box components
- [x] Stability assessment based on the Generalized Nyquist Criteria applicable to MIMO systems
- [x] Computation of stability margins: phase, gain and vector margins
- [x] Oscillation mode identification, bus participation factors, controllability, and observability

## Installation
To use the tool, the following pre-requisites are needed.
1. Python 3.7 or higher together with
   * Numpy (already included in common python packages such as Anaconda)
   * Matplotlib (idem)
   * [PSCAD automation library]([url](https://www.pscad.com/webhelp-v5-al/index.html))
   
   Check the example [here](Examples) for more information on how to install the previous.

2. PSCAD v5 or higher

3. Download or install the Z-tool GitHub repository.

4. Add the location of the downloaded files, and specially the _Source_ folder containing the source code,
to the system path so Python can find the necessary modules:<br />
Environment Variables... -> System variables -> PYTHONPATH -> Directory where _Source_ is located
If PYTHONPATH does not exit, it needs to be created. This should be done automatically if installing the package.


## Usage
Follow the minimal example described [here](Examples). A complete documentation is currently under development.

## Contributors
* Fransciso Javier Cifuentes Garcia (KU Leuven / EnergyVille): Main developer
* Thomas Roose (KU Leuven / EnergyVille): Initial stability analysis functions
* Eros Avdiaj and Özgür Can Sakinci (KU Leuven / EnergyVille): Validation and development support

## Contact Details
For queries about the package or related work please feel free to reach out a 
## Future work
- [ ] Minimum simulation time before starting FFT (does it need to be at least as long as the period of the perturbation or could it be smaller?)
- [ ] Adapt the main function so a previous snapshot can be re-used
- [ ] Switch between current and voltage perturbation, e.g. for the scan of voltage-controlling devices
- [ ] Option to clear the temporary PSCAD files
- [ ] Allow for different computation of the PFs, e.g. admittance PFs
- [ ] Transformation to positive and negative sequence representation
- [ ] Frequency scan and stability analysis optimization based on the passivity properties of the converters

## 