# thresholding - an R package for graph thresholding
## (the name is not finalized yet :-) )
An R package for the [Langston Lab's](https://web.eecs.utk.edu/~mlangsto/) C++ graph thresholding techniques


Navigate to [installation section](#installing-the-thresholding-package)  if [the required system prerequisites](#system-prerequisites) are already installed.


   
## System prerequisites
We assume the user has a version of R installed. While not necessary, having an IDE like [RStudio](https://posit.co/products/open-source/rstudio/) will streamline package installation and utilization.

These thresholding codes depend on the external **igraph** C library for graph creation and manipulation. 
Some external functions from **alglib** are also used and are included in the package already (they are compiled and linked during the package installation process as well)

The following tools were used in development of this package in order to build and link the external 
libraries and functions to the existing lab code. The [Rcpp package](https://www.rcpp.org/) library was used heaviliy in creation of this package.

### Windows users
- #### Rtools
    - Rtools is a toolchain used by R to build source packages (i.e., packages that need code to be compiled within the package).
       - This package wraps an existing C++ command utility for these functionalities and thus needs compiled to be used within R. The [Rcpp](https://www.rcpp.org/) package was used to incorporate a smoother integration of the C++ code into an R package.
       - The C++ code depends on the **[igraph C library]([https://igraph.org/c/)** for network creation and manipulation. The **igraph** library uses the CMake build system to build a static library of **igraph**. This package needs to this static **igraph** library so that the C++ code can link to  **igraph** functions and data structures.
       - Fortunately, CMake is one of the tools that comes with Rtools (as well as normal '**make**')! In theory, Rtools should be the only thing Windows users need to install and use this thresholding package. We are still debugging the installation process on other machines and would love feedback on if the installation works. 
     - **If Rtools is not installed on your system yet, [the link to install any particular version of R is given here.](https://cran.r-project.org/bin/windows/Rtools/)**


### *Nix users
- #### CMake (recommended version 3.0 or higher - 3.27 was used in package development)
   - This package wraps an existing C++ command utility for these functionalities and thus needs compiled to be used within R. 
   - The C++ code depends on the **[igraph C library]([https://igraph.org/c/)** for network creation and manipulation. The **igraph** library uses the CMake build system to build and install **igraph**. This packages needs to build a static library of **igraph** so that the C++ code can link to the **igraph** functions and data structures.
   - If your machine does not come with CMake by default, then you can **[choose a version and install CMake with this link.](https://cmake.org/cmake/help/latest/release/index.html)**
       - This package was developed with CMake version = 3.27 if you aren't sure which version to install.
     
## Installing the thresholding package
1. **Install from devtools()**
   Use the following command in an R window if you have the `devtools` package installed already:
   
      `devtools::install_github("snapplelevi/langston-thresholding")`

   You can use `install.packages(devtools)` in the R window if `devtools` isn't installed already).
   
   There are a few R packages that the **thresholding** package depends on and imports from. After running the `install_github` call, you will likely be prompted if you would like to *update* or *install* some of these packages. **Please install any packages if necessary.** For package updates, you can choose to take any action from the prompt including updating all packages, some packages, or no packages. You can skip the updating step by entering a blank input line (pressing the 'Enter' key once).
   
   `devtools::install_github()`  takes in the same arguments that you would pass to `install.packages()` which allows you to further specify how and where the package should be installed. 
3. **Clone from GitHub and Install as a Source Package**

    1. First, clone the git repo onto your machine using:
       - `git clone https://github.com/snapplelevi/langston-thresholding.git`    (clone with URL over HTTPS)
         
           or
       
       - `git clone git@github.com:snapplelevi/langston-thresholding.git`        (clone using SSH if keys already exist)
         
           or
   
       - `gh repo clone snapplelevi/langston-thresholding`                       (clone using GitHub CLI
       
      


    2. Afterwards, tar the cloned git repo directory using this command:

       `tar -cvzf thresholding.tar.gz ./langston-thresholding`
       
      This creates the tar.gz file of this package for R's `install.packages()` to install from source. This command assumes you are in the same directory of the cloned repo.


    3. After the tar.gz file for the package is created, run the `install.packages()` function with the newly created tar file:
   
       `install.packages("./thresholding.tar.gz", type="source", repos=NULL)`
  
       Note: this command will install to a default library (first entry in R's `.libPaths()`). To move the package installation to a specified location, run the command with:
       
       `install.packages("./thresholding.tar.gz", type="source", repos=NULL, lib="/your/path/here")`
   

   4. After R finishes installing the package, the `thresholding.tar.gz` file can be deleted.

       
**Please contact *lhochste@vols.utk.edu* if you have any trouble with installing the package or if there are any problems experienced while using the package. Thank you for trying out the Langston Lab's thresholding code!**
 
