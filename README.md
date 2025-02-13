# Z-tool
Z-tool is a Python-based implementation for the frequendy-domain analysis of modern power systems.
The core functionalities are admittance characterization and small-signal stability assessment.
The analysis relies on an existing system model in the EMT simulation software [PSCAD]([url](https://www.pscad.com/)).

The following features are currently implemented and validated:
- [x] Voltage perturbation-based admittance scan at several nodes, including MMC-based systems and black-box components
- [x] Stability assessment via Generalized Nyquist Criteria applicable to standalone-stable MIMO systems
- [x] Oscillation mode identification via eigenvalue decomposition (EVD) and bus participation factors
- [x] Passivity assessment and Singular Value Decomposition functions

The flowchart below summarizes a common usage of the tool for stability studies, including frequency-domain system identification ([frequency_sweep](Source/ztoolacdc/frequency_sweep.py)) and several stability analysis functions ([stability_analysis](Source/ztoolacdc/stability.py)):

![Tool flowchart](Doc/flowchart.png)
![Tool summary](Doc/Ztool_summary.png)

## Installation
To use the toolbox, the following pre-requisites are needed.
1. Python 3.7 or higher together with
   * [Numpy](https://numpy.org/), [Scipy](https://scipy.org/), and [Matplotlib](https://matplotlib.org/) (included in common python installations such as Anaconda)
   * [PSCAD automation library]([url](https://www.pscad.com/webhelp-v5-al/index.html))
2. PSCAD v5 or higher is recommended.
3. Install the Z-tool via cmd `py -m pip install ztoolacdc` or using the package files.

## Usage
A generic usage of the package can be summarized in the following steps:
1. Add the Z-tool PSCAD library to your PSCAD project
2. Place the tool's analysis blocks at the target buses and name them uniquely
3. Define the resulting connectivity of the scan blocks
4. Specify the basic simulation settings and frequency range for your study
5. Run the frequency scan and small-signal stability analysis functions
Follow the example(s) described [here](Examples/README.md) for more details. The GUI is currently under development.

## Citing Z-tool
If you find the Z-tool useful in your work, we kindly request that you cite the following publications which you can access [here](https://lirias.kuleuven.be/4201452&lang=en):

```bibtex
@INPROCEEDINGS{Cifuentes2024,
  author={Cifuentes Garcia, Francisco Javier and Roose, Thomas and Sakinci, Özgür Can and Lee, Dongyeong and Dewangan, Lokesh and Avdiaj, Eros and Beerten, Jef},
  booktitle={2024 IEEE PES Innovative Smart Grid Technologies Europe (ISGT EUROPE)}, 
  title={Automated Frequency-Domain Small-Signal Stability Analysis of Electrical Energy Hubs}, 
  year={2024},
  pages={1-6},
  doi={10.1109/ISGTEUROPE62998.2024.10863484}}

```

## Contact Details
For queries about the package or related work please feel free to reach out to [Fransciso Javier Cifuentes Garcia](https://www.kuleuven.be/wieiswie/en/person/00144512). You can find more open-source tools for power systems analysis in the [etch website](https://etch.be/en/research/open-source-tools).

## License
This is a free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. Z-tool is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

## Contributors
* Fransciso Javier Cifuentes Garcia: Main developer
* Thomas Roose: Initial stability analysis functions
* Eros Avdiaj and Özgür Can Sakinci: Validation and support

## Future work
- [x] Exploit the symmetric properties of the admittance matrix for AC (and DC) systems to reduce the scans (less simulation time)
- [x] Allow a previous snapshot to be re-used
- [ ] Snapshot simulation plot
- [ ] Support for non-topology specification: inefficient but easier to use
- [ ] Option to clear the temporary PSCAD files
- [ ] Allow for different computation of the PFs, e.g. admittance PFs
- [ ] Switch between current and voltage perturbation
- [ ] Computation of stability margins: phase, gain and vector margins
<!--- - [ ] Minimum simulation time before starting FFT (does it need to be at least as long as the period of the perturbation or could it be smaller?) --->
<!--- - [ ] Transformation to positive and negative sequence representation 
- [ ] Frequency scan and stability analysis optimization based on the passivity properties of the converters --->