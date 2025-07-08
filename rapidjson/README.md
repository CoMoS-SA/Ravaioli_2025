# DSK - Dystopian Schumpeter meeting Keynes model

First attempt at an SFC version of the DSK model.
This version differs quite substantially from previous ones, having been modified to ensure consistency of the accounting system.
This repository also includes an input files for an empirical calibrations on EU27 data. During the process of empirical calibration, several of the behavioural assumptions were tweaked slightly to enable the model to better reproduce (also quantitatively) certain empirical regularities.

## Compilation and debug runs on Windows machines using WSL and VS Code

This version of the model can be compiled with cmake in WSL (using Ubuntu distro).
To compile it on a Windows machine and perform debug runs in VS Code, try following these steps:
- Make sure that WSL is enabled in Windows features; if your machine is still running WSL1, update to WSL2
- Install Ubuntu from the Microsoft store
- Open the Ubuntu console and execute the following commands:
```
sudo apt-get update
sudo apt-get install build-essential gdb cmake
```
- Install VS code (I am using version 1.64.2)
- In VS code, install the Remote-WSL extension
- Once the WSL extension is installed, open a remote VS Code window using the Ubuntu distro
- Having done so, install the following VS Code extensions in the remote 
    - C/C++ (I am using version 1.7.1; there seems to be a documented problem with doing debugging runs in the most recent version)
    - CMake
    - CMake Tools
- Inside the remote VS Code window, open the folder containing the DSK-SFC model
- Change the path given in the CMakeSettings.json file, "wslPath": "C:\\Users\\severindavid.reissl\\AppData\\Local\\Microsoft\\WindowsApps\\ubuntu.exe" to the appropriate path on your machine

Hopefully you should now be able to compile the model by executing

```
cmake .
make 
```

In the terminal in VS Code. In addition, using the launch configuration contained in the .vscode folder, you should be able to perform debug runs in VSCode by pressing F5 once the model has been compiled.

## Compilation on Windows using MinGW
- [Install Chocolatey](see https://chocolatey.org/install) using an Admin PowerShell (e.g. open via Win + X), then type
```
choco install mingw cmake
```
- Make sure the `Path` environmental variable includes `C:\ProgramData\chocolatey\bin` and `C:\Program Files\cmake\bin` (e.g. you should be able to invoke `gcc` and `cmake` in a command prompt)
- Open a command prompt and cd into DSK source dir, then type
```
del CMakeCache.txt
cmake -G "MinGW Makefiles" -D CMAKE_BUILD_TYPE=Debug .
cmake --build .
```

## Compilation on Linux

```
mkdir build
cd build
cmake ..
make
```

## Command line arguments

The compiled executable dsk_SFC takes several arguments to be entered in the console when the executable is called:

1) The path to a json file containing parameter values, initial values and flags 
2) -r, A name for the run (without spaces), which will be appended, along with the seed, to the name of every output file generated (default="test")
3) -s, The seed, which should be a positive integer (default=1)
4) -f, A dummy, indicating whether full output (1) or reduced output (0) should be written to files (default=0)
5) -c, A dummy indicating whether error messages should be printed to the console (1) or only saved to the error log file (0) (default=0)
6) -v, A dummy indicating whether simulation progress updates should be printed to the console (1=yes, 0=no) (default=0)

Arguments 2-6 are optional; default values will be used if they are not supplied.
To see instructions execute 
```
./dsk_SFC --help 
```
in Linux/macOS or

```
.\dsk_SFC --help 
```
in Windows.

Example Linux/macOS: 
```
./dsk_SFC dsk_sfc_inputs.json -r test -s 1 -f 0 -c 0 -v 0
```
Example Windows: 
```
.\dsk_SFC dsk_sfc_inputs.json -r test -s 1 -f 0 -c 0 -v 0
```

All simulated data will be saved in a folder named "output" in the same directory as the executable dsk_SFC; this folder will be created if it does not exist yet.
All error log files will be saved in a folder named "output/errors" in the same directory as the executable dsk_SFC; this folder will be created if it does not exist yet.

