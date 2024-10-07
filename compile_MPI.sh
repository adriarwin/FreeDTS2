#!/bin/bash

# Weria Pezeshkian
# Niels Bohr International Academy
# Niels Bohr Institute
# University of Copenhagen

cd dts_src

# Compile all .cpp files with MPI support and define the MPI_DETECTED flag
mpic++ -c -O3 -DMPI_DETECTED *.cpp

# Link the object files to create the executable with MPI support
mpic++ -o DTS *.o

# Move the compiled executable to the parent directory
mv DTS ../

cd ..

