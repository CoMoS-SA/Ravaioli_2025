Repository for the version of the DSK (Dystopian Schumpeter meeting Keynes) model described and used in Ravaioli, G., Lamperti, F., Roventini, A. and Domingos, T., 2025. Tackling Emissions and Inequality: Policy Insights from an Agent-Based Model. Available at SSRN: https://ssrn.com/abstract=5223096 or http://dx.doi.org/10.2139/ssrn.5223096

# Licence
Copyright (C) 2025 Giacomo Ravaioli, Francesco Lamperti, Andrea Roventini, Tiago Domingos

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.

# Main files
The folder contains all files necessary to compile the model and simulate the scenarios shown in the paper.

The model is structured as follows:
- dsk_sfc_main.cpp is the main model script which contains the main simulation loop as well as a majority of model functions
- Some functions pertaining to climate, the energy sector, finance, and macroeconomic aspects are contained in separate scripts contained in the "modules" folder. Each .cpp file comes with an associated header file
- dsk_sfc_functions.h declares the model functions
- dsk_sfc_parameters.h declares the model parameters
- dsk_sfc_inits.h declares the model initial values
- dsk_sfc_flags.h declares the model flags (indicator variables used to specify different simulation settings)
- dsk_sfc_include.h declares all libraries and other files needed
- The folder "auxiliary" contains a range of functions used in the code
- ./input_files/dsk_sfc_inputs_baseline.json is the input file which supplies parameters, flags and initial values for the baseline scenario
- runPar.R is an example R script which can be used to perform parallel runs of the model.

# Compilation instructions

## Compilation and debug runs on Windows machines using WSL and VS Code

Compilation and debug runs on Windows machines using WSL and VS Code:
- Make sure that WSL is enabled in Windows features
- Install Ubuntu from the Microsoft store
- Open the Ubuntu console and execute the following commands:
```
sudo apt-get update
sudo apt-get install build-essential gdb cmake
```
- Install VS code
- In VS code, install the Remote-WSL extension
- Once the WSL extension is installed, open a remote VS Code window
- Having done so, install the following VS Code extensions in the remote 
    - C/C++
    - CMake
    - CMake Tools
- Inside the remote VS Code window, open the folder containing the model
- Change the path given in the CMakeSettings.json file, "wslPath" (line 25) to the appropriate path on your machine
- In the file CMakeLists.txt, ensure that the option CMAKE_BUILD_TYPE is set as set(CMAKE_BUILD_TYPE Debug)
- Compile the model by executing:
```
mkdir build
cd build
cmake .
make
```
- Using the launch configuration contained in the .vscode folder, you should be able to perform debug runs in VSCode by pressing F5 once the model has been compiled.

## Compilation on Linux
Compilation in Linux:

- Open the console and type and execute:
```
sudo apt-get update
sudo apt-get install build-essential gdb cmake
```
- frome the source directory, type and execute:
```
mkdir build
cd build
cmake ..
make
```
- The executable file dsk_SFC_G will be created in the build folder

## Compliation in macOS:
- Install Homebrew; open the terminal and execute:
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
- Install cmake; execute:
```
brew install cmake
```
- cd into the DSK-SFC source directory, then type and execute:
```
cmake .
make
```
# Running the model

The compiled executable dsk_SFC_G takes several arguments to be entered in the console when the executable is called:

1) The path to a json file containing parameter values, initial values and flags 
2) -r, A name for the run (without spaces), which will be appended, along with the seed, to the name of every output file generated (default="test")
3) -s, The seed, which should be a positive integer (default=1)
4) -f, A dummy, indicating whether full output (1) or reduced output (0) should be written to files (default=0)
5) -c, A dummy indicating whether error messages should be printed to the console (1) or only saved to the error log file (0) (default=0)
6) -v, A dummy indicating whether simulation progress updates should be printed to the console (1=yes, 0=no) (default=0)

Arguments 2-6 are optional; default values will be used if they are not supplied.
To see instructions, in the same folder where the executable dsk_SFC_G is located execute 
```
./dsk_SFC_G --help 
```
in Linux/macOS or

```
.\dsk_SFC_G --help 
```
in Windows.

Example Linux/macOS (in the same folder where the executable dsk_SFC_G is located): 
```
./dsk_SFC_G ./dsk_SFC_G ../input_files/dsk_sfc_inputs_baseline.json -r test -s 1 -f 0 -c 0 -v 0
```
Example Windows (in the same folder where the executable dsk_SFC_G is located): 
```
.\dsk_SFC_G .\dsk_SFC_G ..\input_files\dsk_sfc_inputs_baseline.json -r test -s 1 -f 0 -c 0 -v 0
```

All simulated data will be saved in a folder named "output" in the current working directory; this folder will be created if it does not exist yet.
All error log files will be saved in a folder named "output/errors" in the current working directory; this folder will be created if it does not exist yet.

# Scenarios

To simulate the scenarios shown in the paper, the relevant flags in the input file provided to the executable must be set to the values given below. All instructions are provided starting from the unchanged dsk_sfc_inputs_baseline.json file found in the input_files folder. To run a different scenario, it is possible to simply change the values inside the dsk_sfc_inputs_baseline.json file, or copy the file with a different name to specify when running the executable.

All scenarios shown in the paper are simulated for 300 periods, of which the first 200 are discarded as a transient during analysis of output data.
The model is simulated 300 times with seeds 1-300 for each Single policies and Additional scenarios, 50 with seeds 1-50 for Policy mixes.

## Baseline run
No changes needed w.r.t. configuration provided.

## Single policies
Single policies, as presented in section 3.1 and 4.1.

- **Progressive income tax**  
  set `flag_change_wage_tax_rate=4`  
  set `aliqw_progressivity_regime_shift=0.48`  

- **Shift taxes to capital**  
  set `flag_change_dividends_tax_rate=2`  
  set `aliqdiv_regime_shift=0.14`  

- **Higher tax Top 10%**  
  set `flag_highly_progressive_taxation=1`  
  set `aliq_households_increase=0.05`  

- **Lower tax Bottom 60%**  
  set `flag_highly_progressive_taxation=3`  
  set `aliq_households_increase=-0.10`  

- **Green capital subsidies**  
  set `flag_environmental_subsidies_C_firms=1`  
  set `env_subsidy_per_machine=0.04`  

- **Dirty capital taxation**  
  set `flag_environmental_subsidies_C_firms=2`  
  set `env_subsidy_per_machine=-0.13`  

- **Carbon tax**  
  set `flag_tax_CO2=6`  
  set `t_CO2_0=0.00045`  
  set `tc2=1.007`

## Policies combination
Polixy mixes, as presented in sections 3.2 and 4.2. 
To construct the Policy Mixes, we activated each policy alternatively and in combination with the others by changing the relative flag value and testing multiple values of the policy parameter.

- **Progressive income tax**  
  set `flag_change_wage_tax_rate=4`  
  set `aliqw_progressivity_regime_shift=0.34` or `aliqw_progressivity_regime_shift=0.48`  

- **Shift taxes to capital**  
  set `flag_change_dividends_tax_rate=2`  
  set `aliqdiv_regime_shift=0.07` or `aliqdiv_regime_shift=0.14`  

- **Higher tax Top 10%**  
  set `flag_highly_progressive_taxation=1`  
  set `aliq_households_increase=0.05`  

- **Lower tax Bottom 60%**  
  set `flag_highly_progressive_taxation=3`  
  set `aliq_households_increase=-0.05`  

- **Green capital subsidies**  
  set `flag_environmental_subsidies_C_firms=1`  
  set `env_subsidy_per_machine=0.025` or `env_subsidy_per_machine=0.045` or `env_subsidy_per_machine=0.075`    

- **Dirty capital taxation**  
  set `flag_environmental_subsidies_C_firms=2`  
  set `env_subsidy_per_machine=-0.06` or `env_subsidy_per_machine=-0.13` or `env_subsidy_per_machine=-0.4`  

- **Carbon tax**  
  set `flag_tax_CO2=6`  
  set `t_CO2_0=0.00025` or `t_CO2_0=0.00045` or `t_CO2_0=0.00075`  
  set `tc2=1.007` 


## Additional scenarios
Selected policy mix and Additional scenarios, as presented in sections 3.3 and 4.3. 
The scenarios presented were obtained as all possible combinations of the following 3.

- **Selected Policy Mix**  
  set `flag_change_wage_tax_rate=4`  
  set `aliqw_progressivity_regime_shift=0.48`   
  set `flag_change_dividends_tax_rate=2`  
  set `aliqdiv_regime_shift=0.14`   
  set `flag_highly_progressive_taxation=3`  
  set `aliq_households_increase=-0.05`  
  set `flag_environmental_subsidies_C_firms=1`  
  set `env_subsidy_per_machine=0.075`  
  set `flag_tax_CO2=6`  
  set `t_CO2_0=0.00045`  
  set `tc2=1.007`

- **Energy transition scenario**  
   set `flag_energy_exp=3`  
   set `K_ge_END_perc=0.7`  
   set `t_length_energy_transition=100`     
   set `renew_impact_on_p_e=0`  

- **Productive green technology scenario**  
   set `flag_correlate_prod_and_green_tech=4`   
   set `rs_uu_lp=0.97`  
   set `rs_uu_ee=0.58`  
   set `rs_uu_ef=0.12`  

Additional scenario presented in Appendix D.3  

- **Managers' bonuses**  
   set `bonuses_share=0.05` or `bonuses_share=0.1` or `bonuses_share=0.2`      