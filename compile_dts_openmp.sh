#!/bin/bash

cd dts_src
g++ -c -O3 -fopenmp *.cpp
g++ -o DTS -fopenmp *.o
mv DTS ../
cd ..
