#!/bin/bash

# Weria Pezeshkian
# Niels Bohr International Academy
# Niels Bohr Institute
# University of Copenhagen

cd dts_src
g++ -c -O3 *.cpp
g++ -o DTS *.o
mv DTS ../
cd ..
