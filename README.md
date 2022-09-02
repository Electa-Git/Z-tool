# Z-tool
Z-tool is a Python-based implementation for the analysis of modern power systems.
The basic functionalities are impedance/admittance scan and stability assessment.
The scan is done based on an already existing project/model in PSCAD according to the preferences of the user.
The following features are currently implemented:
- [x] Voltage perturbation

## Installation
To use the tool, the following pre-requisites are needed
1. Python 3.7 or higher together with
   * Numpy (already included in most python packages such as Anaconda)
   * Matplotlib (idem)
   * [PSCAD automation library]([url](https://www.pscad.com/webhelp-v5-al/index.html))

2. PSCAD v5 or higher

3. Dowload the Z-tool GitHub repository to your PC in a stable location.

4. Add the location of the downloaded files, and specially the _Source_ folder containing the source code,
to the system path so Python can find the necessary modules:<br />
Enviroment Variables... -> System variables -> PYTHONPATH -> Directory where _Source_ is located
If PYTHONPATH does not exit, it needs to be created.

5. If you opt for the MATLAB wrapper....

## Usage
### Python-based example
Follow the example described [here](Examples/README.md).

### MATLAB-based example

## Future work
- [ ] Current/voltage perturbation
- [ ] PSCAD library
- [ ] Stability assessment

## Contributors
