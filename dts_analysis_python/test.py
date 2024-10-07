import numpy as np
import universe

class Namespace:
    def __init__(self, input_file='input.txt', directory_path=None, name_output_files=None,
                 name_output_folder=None, path_output_folder=None, initial_steps=None,
                 final_steps=None):
        self.input_file = input_file
        self.directory_path = directory_path
        self.name_output_files = name_output_files
        self.name_output_folder = name_output_folder
        self.path_output_folder = path_output_folder
        self.initial_steps = initial_steps
        self.final_steps = final_steps

args2=Namespace()



read_only_mode=True
data=universe.Universe(args2,read_only_mode)

print(data.projected_area_array_energy_file[0])

