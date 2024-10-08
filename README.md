# SPM_Project
Final project of the SPM course 2023/2024


# README

To compile and execute the files on my machine, I had to slightly modify a makefile definition in:
  `CXX = g++-13 -std=c++20`
But instead for operating on the shared cluster this definition caused errors due to the explicit version of g++

The `Makefile` provides several commands to perform common tasks. Below is a list of the available commands and their descriptions:

- **make all**  
  Compiles all the executables inside the `src` folder. Creates a folder `obj` to store the object files and avoids compiling each executable named `no_*`.

- **make allmpi**  
  Compiles all the executables inside the `src` folder named `mpi_*`. Used to compile all MPI-related executables.

- **make allnmpi**  
  Compiles all the executables that do not fall under the category compiled by `make allmpi`.

- **make `namefile.cpp`**  
  Compiles only `namefile.cpp` inside the `src` folder.

- **make cleanall**  
  Deletes all object files, compiled executables, and log files.

- **make cleanlog**  
  Deletes all log files inside the `log` folder.

- **make clean**  
  Deletes all object files and compiled executables.

You can add `DEBUG=0` to each compilation command via `MAKE` to change the value of the program's DEBUG flag. Use `DEBUG=0` for efficiency testing and `DEBUG=2` for terminal-readable results. See the code comments for a detailed explanation.

**Recommended compilation command:** 
    `make all DEBUG=2`

Once the file is compiled, run it normally. All executables create a `log` folder and a file inside with a default name where execution details are stored. If the file already exists, the results are appended to it. All executables have a `-h` flag that returns a description of all available options.

The `.sh` files in this folder are some tests for the programs created in this project. `test_general.sh` is a simple test that runs all the executables with default arguments.

**Recommended test:** 
    `./test_general.sh`

Python notebooks may not work depending on the execution environment, they are only used for plotting data.
