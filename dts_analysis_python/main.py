import argparse
import os
import numpy as np
import frame
import matplotlib.pyplot as plt
import time
import re
import universe


def main():
    parser = argparse.ArgumentParser(description="Process input parameters")
    parser.add_argument("-i", "--input_file", help="Path to input file", required=True)
    parser.add_argument("-d", "--directory_path", help="Path to data directory")
    parser.add_argument("-o", "--name_output_files", help="Name for output files")
    parser.add_argument("-f", "--name_output_folder", help="Name for output folder")
    parser.add_argument("-p", "--path_output_folder", help="Path to output folder")
    parser.add_argument("-e", "--initial_steps", help="Equilibration steps")
    parser.add_argument("-s", "--final_steps", help="Number of final steps")
    """
    parser.add_argument("-a", "--area", help="Area calculation", choices=["on", "off"])
    parser.add_argument("-n", "--inclusion_average_neighbours", help="Inclusion average neighbours calculation", choices=["on", "off"])
    parser.add_argument("-c", "--inclusion_cluster_statistics", help="Inclusion cluster statistics calculation", choices=["on", "off"])
    parser.add_argument("-m", "--membrane_thickness", help="Membrane thickness calculation", choices=["on", "off"])
    parser.add_argument("-g", "--energy", help="Energy calculation", choices=["on", "off"])"""
    
    import time

    start_time = time.time()
    # Your code or task here
    

    args = parser.parse_args()

    print(args)

    read_only_mode=False
    calculation = universe.Universe(args,read_only_mode)


    end_time = time.time()

    elapsed_time = end_time - start_time
    print("Elapsed time:", elapsed_time, "seconds")



if __name__ == "__main__":
    main()
