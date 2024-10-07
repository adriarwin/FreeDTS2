#!/bin/bash

# Weria Pezeshkian
# Niels Bohr International Academy
# Niels Bohr Institute
# University of Copenhagen

cd dts_src
mpic++ -c -O3 *.cpp
mpic++ -o DTS *.o
mv DTS ../
