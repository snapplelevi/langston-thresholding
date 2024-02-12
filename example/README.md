Tenative things to come back to:
- Figure out solution for unique names when multiple methods are passed to the analysis function
  - Specifying two iterative methods in the package's current state will lead to writing to the same
  iterative file, which would make the results unreadable and not helpful to the user.
- Documentation / error message output for igraph_ext_alg and igraph_ext_io
- Documentation and other helper functions with thresholdAnalysis
- Updating or modifying thresholding::help() to get more specific instructions for thresholdAnalysis
or expanding to multiple help() functions depending on the function
- Test installation for Mac users? Does this work on *Nix systems? Other Windows environments?
- Create installation and set-up instructions for user convience in README.md
- Update README.md to have proper outline of the package
- **Structure package so that it can be submitted to CRAN potentially (get rid of tar files, executables, etc; shorten directory names)
- Error check the string argument parsing for invalid input in the str_methods
  - (i.e. text instead of numbers, expressions instead of numbers)
  - eventually could add the feature to type a certain method instead of use numbers, 
    but for now this will be the way the function will be written
- Add note in main README.md for troubleshooting such as:
  - If there is an issue with CMakeCache (naivgate to and delete from igraph/build dir)
  - Attempt a different method of installation if devtools::install_github() doesn't work
    (.tar file extraction with devtools::install(), install.packages(), etc.)
  - Reach out for any issues with the build process, installation of tools, or any general questions
  
