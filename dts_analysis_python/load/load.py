import os
import numpy as np

class Load:
    def __init__(self, universe_instance):
        self.universe = universe_instance

    def load_data(self):
        # Helper function to load npy file with extension if needed
        def np_load_with_extension(file_path):
            if not file_path.endswith('.npy'):
                file_path += '.npy'
            return np.load(file_path)

        def check_file_existence(file_path):
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"File '{file_path}' not found.")
        
        # Load frame-related data
        if self.universe.frame_iteration:
            file_path_frame_num = os.path.join(self.universe.output_folder_path, self.universe.frame_steps_filename + self.universe.name_output_files)
            check_file_existence(file_path_frame_num+'.npy')
            self.universe.frame_num_list = np_load_with_extension(file_path_frame_num)

            if self.universe.area_calculation == "on":
                file_path_area = os.path.join(self.universe.output_folder_path, self.universe.area_filename + self.universe.name_output_files)
                file_path_projected_area = os.path.join(self.universe.output_folder_path, self.universe.projected_area_filename + self.universe.name_output_files)
                check_file_existence(file_path_area+'.npy')
                check_file_existence(file_path_projected_area+'.npy')
                self.universe.area_array = np_load_with_extension(file_path_area)
                self.universe.projected_area_array = np_load_with_extension(file_path_projected_area)

            if self.universe.inclusion_average_neighbours_calculation == "on":
                file_path_inclusion_avg = os.path.join(self.universe.output_folder_path, self.universe.inclusion_average_neighbours_filename + self.universe.name_output_files)
                check_file_existence(file_path_inclusion_avg+'.npy')
                self.universe.inclusion_average_neighbours_array = np_load_with_extension(file_path_inclusion_avg)
            
            if self.universe.inclusion_cluster_statistics_calculation == "on":
                file_path_inclusion_cluster = os.path.join(self.universe.output_folder_path, self.universe.inclusion_cluster_filename + self.universe.name_output_files)
                check_file_existence(file_path_inclusion_cluster+'.npy')
                self.universe.inclusion_cluster_statistics_array = np_load_with_extension(file_path_inclusion_cluster)

            if self.universe.fluctuation_spectrum_planar_calculation == "on":
                file_path_fluctuation_spectrum = os.path.join(self.universe.output_folder_path, self.universe.fluctuation_spectrum_filename + self.universe.name_output_files)
                file_path_qvec = os.path.join(self.universe.output_folder_path, self.universe.qvec_filename + self.universe.name_output_files)
                check_file_existence(file_path_fluctuation_spectrum+'.npy')
                check_file_existence(file_path_qvec+'.npy')
                self.universe.fluctuation_spectrum_planar_array = np_load_with_extension(file_path_fluctuation_spectrum)
                self.universe.q2vec_array = np_load_with_extension(file_path_qvec)

            if self.universe.membrane_thickness_calculation == "on":
                file_path_membrane_thickness = os.path.join(self.universe.output_folder_path, self.universe.membrane_thickness_filename + self.universe.name_output_files)
                check_file_existence(file_path_membrane_thickness+'.npy')
                self.universe.membrane_thickness_array = np_load_with_extension(file_path_membrane_thickness)

        # Load non-frame-related data
        if self.universe.non_frame_iteration:
            file_path_energy = os.path.join(self.universe.output_folder_path, self.universe.energy_filename + self.universe.name_output_files)
            file_path_energy_steps = os.path.join(self.universe.output_folder_path, self.universe.energy_steps_filename + self.universe.name_output_files)
            file_path_projected_area_energy = os.path.join(self.universe.output_folder_path, self.universe.projected_area_energy_filename + self.universe.name_output_files)
            check_file_existence(file_path_energy +'.npy')
            check_file_existence(file_path_energy_steps+'.npy')
            check_file_existence(file_path_projected_area_energy+'.npy')
            self.universe.energy_array = np_load_with_extension(file_path_energy)
            self.universe.energy_MCsteps_array = np_load_with_extension(file_path_energy_steps)
            self.universe.projected_area_array_energy_file = np_load_with_extension(file_path_projected_area_energy)

    def load_data_pt(self):
        # Helper function to load npy file with extension if needed
        def np_load_with_extension(file_path):
            if not file_path.endswith('.npy'):
                file_path += '.npy'
            return np.load(file_path)

        def check_file_existence(file_path):
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"File '{file_path}' not found.")
            
        if self.universe.non_frame_iteration and self.universe.parallel_tempering_on:
            # Load tempering moves
            file_path_tempering_moves = os.path.join(self.universe.output_folder_path, self.universe.tempering_moves_filename + self.universe.name_output_files)
            check_file_existence(file_path_tempering_moves + '.npy')
            self.universe.tempering_moves_array = np_load_with_extension(file_path_tempering_moves)

            # Load beta list
            file_path_beta_list = os.path.join(self.universe.output_folder_path, self.universe.beta_list_filename + self.universe.name_output_files)
            check_file_existence(file_path_beta_list + '.npy')
            self.universe.beta_list = np_load_with_extension(file_path_beta_list)

            # Load each energy array individually
            self.universe.energy_array = []
            i = 0
            while True:
                file_path_energy = os.path.join(self.universe.output_folder_path, f"{self.universe.energy_filename}_{i}") + self.universe.name_output_files + '.npy'
                if os.path.exists(file_path_energy):
                    self.universe.energy_array.append(np_load_with_extension(file_path_energy))
                    i += 1
                else:
                    break

            # Load each energy steps array individually
            self.universe.energy_MCsteps_array = []
            i = 0
            while True:
                file_path_energy_steps = os.path.join(self.universe.output_folder_path, f"{self.universe.energy_steps_filename}_{i}") + self.universe.name_output_files + '.npy'
                if os.path.exists(file_path_energy_steps):
                    self.universe.energy_MCsteps_array.append(np_load_with_extension(file_path_energy_steps))
                    i += 1
                else:
                    break

            # Load each projected area energy file array individually
            self.universe.projected_area_array_energy_file = []
            i = 0
            while True:
                file_path_projected_area_ef = os.path.join(self.universe.output_folder_path, f"{self.universe.projected_area_energy_filename}_{i}") + self.universe.name_output_files + '.npy'
                if os.path.exists(file_path_projected_area_ef):
                    self.universe.projected_area_array_energy_file.append(np_load_with_extension(file_path_projected_area_ef))
                    i += 1
                else:
                    break


