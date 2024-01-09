# thresholding
R package for Langston Lab's C++ thresholding techniques
Package still in development... 

For R installation (also still figuring out the right command to build):
1. devtools::install_github("snapplelevi/thresholding")  // use install.packages(devtools) in the R terminal if devtools isn't installed already)
2. First, clone the git repo onto your machine.

   Afterwards, tar the clone git repo directory (**tar -cvfz thresholding.tar.gz ./thresholding**).

   then, run the install.packages function with the newly created tar file:
   
   **install.packages("./thresholding.tar.gz", type="source", repos=NULL)**  
   Note: this command will install to a default library (first entry in R's .libPaths()).
   
   To move this to a desired location, run the command with
   
   **install.packages("./thresholding.tar.gz", type="source", repos=NULL, lib="/your/path/here")**
   
## Required tools for installation
These thresholding codes depend on the external **igraph** C library for graph creation and manipulation. 
Some external functions from **alglib** are also used. 

The following tools were used in development of this package in order to build and link the external 
libraries and functions to the existing lab code. The Rcpp library was used heaviliy in creation of this
R package.

### Windows users
- Rtools (version #)
    - Purpose of Rools
- CMake (minimum version 3.18)
    - CMake is required to create the shared object/.dll file from the igraph external dependency.
    - When R builds the package, it will run the Makevars.win in **src**, which runs the igraph CMake calls.

### *Nix users
- CMake (minimum version 3.18)
    - CMake is required to create the shared object/.dll file from the igraph external dependency.
    - When R builds the package, it will run the Makevars file in **src** which runs the igraph CMake calls.
- Current Issue: building igraph during package installation
   - Error during linking process of shared object that asks to recompiled certain files with -fPIC (position independent code) flag enabled.
     We are still looking for how to resolve this. This problem has happened during the install.packages() testing, but devtools::install_github()
     hasn't been tested yet. 
