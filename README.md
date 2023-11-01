# thresholding
R package for Langston Lab's C++ thresholding techniques
Package still in development... 

For R installation (also still figuring out the right command to build):
1. devtools::install_github("snapplelevi/threhsolding")  // use install.packages(devtools) if devtools does not exist already)
2. First, tar the clone git repo directory (**tar -cvfz thresholding.tar.gz ./thresholding**).

   then, run the install.packages function with the newly created tar file:
   
   **install.packages("./thresholding.tar.gz", type="source", repos=NULL)**  
   Note: this command will install to a default library (first entry in R's .libPaths()).
   To move this to a desired location, run the command with
   **install.packages("./thresholding.tar.gz", type="source", repos=NULL, lib="/your/path/here)**
## Required tools for installation
These thresholding codes depend on the external **igraph** C library for graph creation and manipulation. 
Some external functions from **alglib** are also used. 

The following tools were used in development of this package in order to build and link the external 
libraries and functions to the existing lab code. The Rcpp library was used heaviliy in creation of this
R package.

### Windows users
- Rtools (version #)
    - Purpose of Rools
- CMake (version #)
    - Purpose of CMake

### *Nix users
- List tools needed
