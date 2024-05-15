# External Dependencies

## IGraph C library tar file
**igraph** version 0.9.9 will be contained in the file `igraph-0.9.9.tar.xz`. This file is untarred so the igraph library can be compiled during the package build process.

If compiling on Windows, R will read `src/Makefile.win` first for preliminary building instructions. Otherwise, it will use the `src/Makevars` file.

This package also uses the **alglib** library, which is provided in the `./alglib` subdirectory.