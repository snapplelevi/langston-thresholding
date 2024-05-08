# thresholding
R package for Langston Lab's C++ thresholding techniques
\
\
Navigate to installation section if system prerequisites are in place [**Installing the thresholding package**](#installing-the-thresholding-package) 


   
## Required tools for before package installation
These thresholding codes depend on the external **igraph** C library for graph creation and manipulation. 
Some external functions from **alglib** are also used and are included in the package already (they are compiled and linked during the package installation process as well)

The following tools were used in development of this package in order to build and link the external 
libraries and functions to the existing lab code. The [Rcpp package](https://www.rcpp.org/) library was used heaviliy in creation of this package.

### Windows users
- Rtools
    - Rtools is a toolchain used by R to build source packages (i.e., packages that need code to be compiled within the package).
       - This package wraps an existing C++ command utility for these functionalities and thus needs compiled to be used within R. The [Rcpp](https://www.rcpp.org/) package was used to incorporate a smoother integration of the C++ code into an R package.
       - The C++ code depends on the **[igraph C library]([https://igraph.org/c/)** for network creation and manipulation. The **igraph** library uses the CMake build system to build and install **igraph**. This packages needs to build a static library of **igraph** so that the C++ code can link to the **igraph** functions and data structures.
       - Fortunately, CMake is one of the tools that comes with Rtools (as well as normal '**make**')! In theory, Rtools should be the only thing Windows users need to install and use this thresholding package. We are still debugging the installation process on other machines and would love feedback on if the installation works. 
     - **If Rtools is not installed on your system yet, [the link to install any particular version of R is given here.](https://cran.r-project.org/bin/windows/Rtools/)**


### *Nix users
- CMake (recommended version 3.0 or higher - 3.27 was used in package development)
   - This package wraps an existing C++ command utility for these functionalities and thus needs compiled to be used within R. 
   - The C++ code depends on the **[igraph C library]([https://igraph.org/c/)** for network creation and manipulation. The **igraph** library uses the CMake build system to build and install **igraph**. This packages needs to build a static library of **igraph** so that the C++ code can link to the **igraph** functions and data structures.
   - If your machine does not come with CMake by default, then you can **[choose a version and install CMake with this link.](https://cmake.org/cmake/help/latest/release/index.html)**
       - This package was developed with CMake version = 3.27 if you aren't sure which version to install.
     
## Installing the thresholding package
1. **Install from devtools()**
   Use the following command in an R window if you have the `devtools` package installed already:
      `devtools::install_github("snapplelevi/thresholding")`

   You can use `install.packages(devtools)` in the R window if `devtools` isn't installed already)
4. **Install from GitHub source)**

    1. First, clone the git repo onto your machine using:
       - `git clone https://github.com/snapplelevi/thresholding.git`    (clone with URL over HTTPS)\
           or
       - `git clone git@github.com:snapplelevi/thresholding.git`        (clone using SSH if keys already exist)\
           or
       - `gh repo clone snapplelevi/thresholding`                       (clone using GitHub CLI)

    2. Afterwards, tar the cloned git repo directory using `tar -cvzf thresholding.tar.gz ./thresholding`.\
       \
      This creates the source tar.gz file for R's `install.packages()` to install from.

    4. After the tar.gz file of the package is createrun the install.packages function with the newly created tar file:
   
       `install.packages("./thresholding.tar.gz", type="source", repos=NULL)`\
       \
       Note: this command will install to a default library (first entry in R's `.libPaths()`).\
      \
       To move this to a desired location, run the command with:
      \
       `install.packages("./thresholding.tar.gz", type="source", repos=NULL, lib="/your/path/here")`


\ 
\ 
**Please contact *lhochste@vols.utk.edu* if you have any trouble with installing the package or if there are any bugs in using the package. Thank you for trying out the Langston Lab's thresholding code!**
 
