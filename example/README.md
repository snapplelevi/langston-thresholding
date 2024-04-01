Tenative things to come back to:

- SEE COMMENTS IN src/R_main_analysis_driver.cpp

- UPDATE THE USER FUNCTIONALITY OF THE INDIVIDUAL get_iter_t_vals(), etc. FUNCTIONS
TO RETURN BOTH THE DF AND THE CORRESPONDING D LIST AS ONE RETURN LIST

- UPDATE THE PLOTTING FUNCTIONALITY TO WORK FROM A FILE NAME INSTEAD OF MULTIPLE PIECES

- Make sure the package builds with igraph as a .xz file (per CRAN docs)

- Figure out solution for unique names when multiple methods are passed to the analysis function (currently using the PID string for now)

- Figure out a more unique name than "thresholding"
  iterative file, which would make the results unreadable and not helpful to the user.
  
  grthresh, graphthresh, lgthresh (Langston graph thresholding),
  
  
- Add additional parameters for file naming inputs (like outfile_prefix)
  - Want to give users the versatility to rename their files while still
    being able to use the package. Right now the code is set up to run only
    if the user is following the .iterative.txt name. Maybe choose an optional 
    argument that controls whether or not multiple files are grabbed just based     on the prefix (or just do that by default)
    
- Documentation / error message output for igraph_ext_alg and igraph_ext_io

- Documentation and other helper functions with analysis
  - Do this once the thresholding paper comes in to have more concrete knowledge

- Documentation/examples for utility functions


- Updating or modifying thresholding::help() to get more specific instructions for thresholdAnalysis or expanding to multiple help() functions depending on the function

- Test installation for Mac users? Does this work on *Nix systems? Other Windows environments?

- Create installation and set-up instructions for user convenience in README.md

- Update README.md to have proper outline of the package

- **Structure package so that it can be submitted to CRAN potentially (get rid of tar files, exacutables, etc; shorten directory names)

- Error check the string argument parsing for invalid input in the str_methods
  - Confirm lists of system requirements
  - (i.e. text instead of numbers, expressions instead of numbers)
  - eventually could add the feature to type a certain method instead of use numbers, 
    but for now this will be the way the function will be written
    
- Add note in main README.md for troubleshooting such as:
  - If there is an issue with CMakeCache (navigate to and delete from igraph/build dir)
  - Attempt a different method of installation if devtools::install_github() doesn't work
    (.tar file extraction with devtools::install(), install.packages(), etc.)
  - Reach out for any issues with the build process, installation of tools, or any general questions
  
