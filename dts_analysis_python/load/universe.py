import argparse
import os
import numpy as np
import frame
import load
import matplotlib.pyplot as plt
import time
import re


class Universe:
    def __init__(self, input_arguments, read_only_mode=False):

        #predefined variables
        self.area_filename='area'
        self.projected_area_filename='projected_area'
        self.membrane_thickness_filename='thickness'
        self.fluctuation_spectrum_filename='fluctuation_spectrum'
        self.energy_filename='energy'
        self.energy_steps_filename='energy_MCsteps'
        self.projected_area_energy_filename='projected_area_ef'
        self.inclusion_average_neighbours_filename='inclusion_average_neighbours'
        self.inclusion_cluster_filename='inclusion_cluster_statistics'
        self.qvec_filename='q2vec'
        self.frame_steps_filename='frame_steps'
        self.output_analysis_filename='analysis_output.txt'
        self.tempering_moves_filename='tempering_moves'
        self.beta_list_filename='beta_list'
        

        #Paths and names of files constructed from others
        self.path_input_dts=None
        self.output_folder_path=None


        ############ BEGINNING OF THE INPUTS from input.txt ###

        #Initial and final frames subject to study
        self.initial_step = None
        self.final_step = None

        #names of directory where analysis is done
        self.directory_path = None

        #name of the output files and the output folder
        self.name_output_files = ''
        self.name_output_folder = None

        #name of TrjTSI folder and name of the output energy file
        self.name_TrjTSI_folder='TrjTSI'
        self.name_energy_filename='output-en.xvg'
        self.name_TrjTSI_files='dts'

        #temperature of your system
        self.beta=1

        #Frame data
        self.area_calculation = 'off'
        self.inclusion_average_neighbours_calculation = 'off'
        self.inclusion_cluster_statistics_calculation = 'off'
        self.fluctuation_spectrum_planar_calculation = 'off'
        self.membrane_thickness_calculation = 'off'
        #Non-Frame data
        self.projected_area_calculation = 'off'
        self.energy_calculation = 'off'
        #Parallel tempering method on or off
        self.parallel_tempering = 'off'
        self.version=1
        

        #######END of the INPUTS from input.txt#######

        #Frame-related-varaibles
        self.frame_path_list=None
        self.nframes=None
        self.frame_num_list=None
        self.ninclusion=0

        #Are we performing frame calculations or not?
        self.frame_iteration=False
        self.non_frame_iteration=False

        #Variables for calculation of fluctuation spectrum
        self.bx=None
        self.by=None

        #IDK what this does!

        self.fluctuation_spectrum_no_inclusions = "off"

        #Different parameters which are inputed with the input.dts file.
        #membrane parameters
        self.Spont_C=0.0
        self.kappa=20.0
        self.kappag=0.0
        #inclusion parameters
        #A number that tells how many types of inclusions (N) we have got
        self.number_of_inclusion_type=None
        #A list of type [[1,'S'],[2,'A'],...,[N,'S']] which tells the assigned
        #number to each inclusion and if it is symmetric ('S') or asymmetric ('A')
        self.inclusion_type=[]
        #Inclusion definition. A list of type [[K,Kg,C0],[Kp,Kl,Cp,Cl],...,[Kp,Kl,Cp,Cl]] depending on
        #the data of each inclusion
        self.inclusion_definition=[]
        #Inclusion density. A list of type [d1,d2,d3], where d1 is the inclusion density of inclusion type 1..
        self.inclusion_density=[]
        self.inclusion_density_in_input=False

        #Variable that tell us if there are inclusion
        self.inclusion=False

        #Parallel tempering related variables
        self.parallel_tempering_on=False
        self.name_TrjTSI_folder_pt='TrjTSI_temp_'
        self.name_energy_file_pt='output_temp_'
        self.name_tempering_moves_file='tempering_moves.txt'
        self.TrjTSI_path_list = []
        self.energy_files_list = []
        self.beta_list = None
        self.tempering_moves_array = None
        self.number_of_temperatures=None
        self.parallel_tempering_tau = None
        

        
        

        #Arrays that store data from frames
        self.area_array= None
        self.projected_area_array= None
        self.inclusion_average_neighbours_array = None
        self.inclusion_cluster_statistics_array = None
        self.fluctuation_spectrum_planar_array = None
        self.membrane_thickness_array=None

        #Arrays that store data from output energy file
        self.projected_area_array_energy_file=None
        self.energy_array=None
        self.energy_MCsteps_array=None

        #I do not know what this counter does!
        self.stupid_counter=0
        
        
        

        #assigning the values of the input file to our variables.
        self.parse_input_file(input_arguments.input_file)

        #overwriting with terminal inputs
        for attr, value in vars(input_arguments).items():
            if value is not None:
                setattr(self, attr, value)
        
        #Defining self.frame_variables and self.non_frame_variables, for the case
        #in which only non_frame_variables want to be worked out.

        self.frame_variables=[self.area_calculation,self.inclusion_average_neighbours_calculation, \
                              self.inclusion_cluster_statistics_calculation,self.fluctuation_spectrum_planar_calculation,
                              self.membrane_thickness_calculation]
        self.non_frame_variables=[self.projected_area_calculation,self.energy_calculation]

        if "on" in self.frame_variables:
            self.frame_iteration=True

        if "on" in self.non_frame_variables:
            self.non_frame_iteration=True

        self.path_input_dts=os.path.join(self.directory_path,"input.dts")

        if self.parallel_tempering=="on":
            self.parallel_tempering_on=True
            self.extract_pt_parameters_dts()

        #Reading input.dts. Should be improved as to include all the cases.
        
        self.extract_membrane_parameters_dts()
        self.extract_inclusion_parameters_dts()
        

        self.output_folder_path = os.path.join(self.directory_path, self.name_output_folder)


        if read_only_mode==False:
            #Creation of folder where I will keep my data (always to be done)
            # Create the folder if it doesn't exist
            if not os.path.exists(self.output_folder_path):
                os.makedirs(self.output_folder_path)


            if self.parallel_tempering_on==False:

                if self.frame_iteration==True:
                    self.generate_frame_path_list()
                    self.array_initialization_frames()
                    self.iteration_frames()
                
                if self.non_frame_iteration==True:
                    path=os.path.join(self.directory_path,self.name_energy_filename)
                    self.energy_MCsteps_array, self.energy_array,self.projected_area_array_energy_file=self.read_energy_file(path)

                self.save_data()

            elif self.parallel_tempering_on==True:

                temp_folders,energy_files,temperature_list=self.sort_and_order_lists(self.name_TrjTSI_folder_pt,self.name_energy_file_pt)
                self.number_of_temperatures=len(temperature_list)

                self.beta_list=np.array(temperature_list)
                self.energy_files_list = [os.path.join(self.directory_path, file) for file in energy_files]
                self.TrjTSI_path_list = [os.path.join(self.directory_path, folder) for folder in temp_folders]

                print(self.beta_list)
                print(self.energy_files_list)

                self.tempering_moves_array=self.read_tempering_moves(self.name_tempering_moves_file)

                if self.non_frame_iteration==True:
                    self.read_energy_file_pt(self.energy_files_list)

                self.save_data_pt()
                    
                


            

        elif read_only_mode==True:

            if self.parallel_tempering_on==False:
                loader = load.Load(self)
                loader.load_data()

            elif self.parallel_tempering_on==True:
                loader = load.Load(self)
                loader.load_data_pt()



    def generate_frame_path_list(self):

        self.frame_path_list=np.array([],dtype=str)
        #Checking existance of frames 
        for i in range(self.initial_step,self.final_step+1):
            
            #Checking existence of frame i
            file_path_existence,file_path=self.check_frame_existence(i)

            if file_path_existence:
                # Check if the file is empty
                if os.path.getsize(file_path) == 0:
                    print(f"File {file_path} exists but is empty.")
                else:
                    self.frame_path_list = np.append(self.frame_path_list, file_path)
                
        #Terminating program if there are no frames

        if self.frame_path_list.size==0:
            raise ValueError("There are no tsi files in the indicated path and tsi files range")   
        
        else:
            self.nframes=self.frame_path_list.size
            self.frame_num_list=range(0,self.nframes)
        
        self.frame_num_list=np.array(self.frame_num_list,dtype=int)

        #Checking number of inclusions is different than zero if one has calculations.
        frame_object=frame.frame(self.frame_path_list[0])
        self.ninclusion=frame_object.ninclusion
        if self.ninclusion==0:
            self.inclusion==False
            self.inclusion_average_neighbours_calculation = None
            self.inclusion_cluster_statistics_calculation = None

            if self.fluctuation_spectrum_planar_calculation=="on":
                self.fluctuation_spectrum_no_inclusions="on"
        
        del frame_object
        

    def extract_pt_parameters_dts(self):

        with open(self.path_input_dts, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith("Parallel_Tempering"):
                    _, value = line.split('=')
                    self.parallel_tempering_tau=int(value.split()[1])


    def extract_membrane_parameters_dts(self):

        with open(self.path_input_dts, 'r') as file:
            for line in file:
                line = line.strip()

                    # Extract kappa
                if line.lstrip().startswith("Kappa"):
                    _, value = line.split('=')
                    floats = [float(x) for x in value.split()]
                    print(floats)
                    self.kappa = floats[0]
                    self.kappag = floats[1]

                    # Extract Spont_C
                elif line.lstrip().startswith("Spont_C"):
                    _, value = line.split('=')
                    self.Spont_C = float(value)

    def extract_inclusion_parameters_dts(self):
        inclusion_section = False
        inclusion_interactions_section = False
        

        with open(self.path_input_dts, 'r') as file:
            for line in file:
                line = line.strip()

                # Detect the start of the INCLUSION section
                if line == "INCLUSION":
                    inclusion_section = True
                    inclusion_interactions_section = False
                    self.inclusion=True
                    continue

                # Detect the start of the Inclusion-Inclusion-Int section
                if line == "Inclusion-Inclusion-Int":
                    inclusion_interactions_section = True
                    inclusion_section = False
                    continue

                

                # Process inclusion parameters
                if inclusion_section:
                    # Detect number of inclusions
                    if line.startswith("Define"):
                        self.number_of_inclusion_type = int(line.split()[1])
                        continue

                    # Skip the header line
                    if line.startswith("SRotation"):
                        continue

                    # Detect inclusion types and their parameters
                    if re.match(r"^\d", line):
                        params = line.split()
                        inclusion_id = int(params[1][3]) # Extract the type (e.g., Pro1, Pro2)
                        
                        print(inclusion_id)
                        # Determine if the inclusion is symmetric or asymmetric
                        k_values = [float(params[2])-self.kappa, float(params[3]), float(params[6])]
                        kp_values = [float(params[4]), float(params[5]), float(params[7]), float(params[8])]

                        print(k_values)
                        print(kp_values)

                        if any(k_values):
                            inclusion_symmetry = 'S'
                            inclusion_data = k_values
                        elif any(kp_values):
                            inclusion_symmetry = 'A'
                            inclusion_data = kp_values
                        else:
                            inclusion_symmetry = 'U'  # Default if neither condition is met
                            inclusion_data = k_values

                        self.inclusion_type.append([inclusion_id, inclusion_symmetry])

                        # Store inclusion data
                        self.inclusion_definition.append(inclusion_data)

                    # Extract density
                    elif line.startswith("Density"):
                        self.inclusion_density_in_input=True
                        densities = line.split()[1:]
                        self.inclusion_density = [float(d) for d in densities]                   


        

    
    def write_output_file(self,output_path):
        content = f"Initial step={self.initial_step}\n" \
              f"Final step={self.final_step}\n" \
              f"kappa={self.kappa}\n" \
              f"kappa_g={self.kappag}\n" \
              f"Spont_C={self.Spont_C}\n" \
              
        
              
        if self.inclusion==True:
            content+=f"INCLUSIONS\n"
            for i in range(0,self.number_of_inclusion_type):
                string0=f"Inclusion type {i+1}\n"
                string1=f"Type={self.inclusion_type[i][0]},{self.inclusion_type[i][1]}\n"
                if (self.inclusion_type[i][1]=="S" or self.inclusion_type[i][1]=="U"):
                    string2=f"Data={self.inclusion_definition[i][0]},{self.inclusion_definition[i][1]},{self.inclusion_definition[i][2]}\n" 
                elif self.inclusion_type[i][1]=="A":
                    string2=f"Data={self.inclusion_definition[i][0]},{self.inclusion_definition[i][1]},{self.inclusion_definition[i][2]},{self.inclusion_definition[i][3]}\n"
                
                content+=string0+string1+string2

                if self.inclusion_density_in_input==True:
                    string3=f"Density={self.inclusion_density[i]}\n"
                    content+=string3

                

        # Write the content to the file
        with open(output_path, 'w') as file:
            file.write(content)


    def read_energy_file(self, path):
        print('obtaining energy')
        try:
            with open(path, 'r') as f:
                lines = f.readlines()
        except FileNotFoundError:
            print("File not found:", path)
            return

        # Initialize lists to store data
        mcstep = []
        energy = []
        projected_area = []

        # Parse each line and extract data
        for line in lines:
            # Skip lines starting with '#' or empty lines
            if line.startswith('#') or not line.strip():
                continue
            # Split line by whitespace and extract the first two columns
            columns = line.split()
            if len(columns) < 2:
                print("Skipping line with unexpected format:", line.strip())
                continue
            try:
                mcstep.append(int(columns[0]))
                energy.append(float(columns[1]))
                projected_area.append(float(columns[2])*float(columns[3]))
            except (IndexError, ValueError):
                print("Skipping line with unexpected format:", line.strip())

        # Convert lists to numpy arrays

        
        return np.array(mcstep), np.array(energy), np.array(projected_area)
    
    def read_energy_file_pt(self, energy_files_list):

        self.energy_MCsteps_array=[]
        self.energy_array=[]
        self.projected_area_array_energy_file=[]
        for i in range(0,len(energy_files_list)):

            mcstep1,energy1,projected_area1=self.read_energy_file(energy_files_list[i])
            self.energy_MCsteps_array.append(mcstep1)
            self.energy_array.append(energy1)
            self.projected_area_array_energy_file.append(projected_area1)




    def parse_input_file(self, input_file):
        with open(input_file, 'r') as file:
            for line in file:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                variable, value = line.split('=')
                variable = variable.strip()
                value = value.strip()

                if variable == 'fluctuation_spectrum_planar_calculation':
                    auxiliar=value.split()
                    value=auxiliar[0]
                    self.bx=int(auxiliar[1])
                    self.by=int(auxiliar[2])
                    
                
                if variable == 'initial_step':
                    value=int(value)

                if variable == 'final_step':
                    value=int(value)
                

                setattr(self, variable, value)

    def save_data(self):

        #Save data
        file_path=os.path.join(self.output_folder_path,self.output_analysis_filename) + self.name_output_files

        self.write_output_file(file_path)
        
    
        

        if self.frame_iteration==True:
            file_path=os.path.join(self.output_folder_path,self.frame_steps_filename)+self.name_output_files
            np.save(file_path,self.frame_num_list)
            if self.area_calculation=="on":
                file_path_area=os.path.join(self.output_folder_path,self.area_filename)+self.name_output_files
                file_path_projected_area=os.path.join(self.output_folder_path,self.projected_area_filename)+self.name_output_files
                np.save(file_path_area,self.area_array)
                np.save(file_path_projected_area,self.projected_area_array)

            
            if self.inclusion_average_neighbours_calculation=="on":
                file_path=os.path.join(self.output_folder_path,self.inclusion_average_neighbours_filename)+self.name_output_files
                np.save(file_path,self.inclusion_average_neighbours_array)
            
            if self.inclusion_cluster_statistics_calculation=="on":
                #Since maximum cluster size is the number of inclusions:
                file_path=os.path.join(self.output_folder_path,self.inclusion_cluster_filename)+self.name_output_files
                np.save(file_path,self.inclusion_cluster_statistics_array)

            if self.fluctuation_spectrum_planar_calculation=="on":
                file_path_fs=os.path.join(self.output_folder_path,self.fluctuation_spectrum_filename)+self.name_output_files
                file_path_qvec=os.path.join(self.output_folder_path,self.qvec_filename)+self.name_output_files
                np.save(file_path_fs,self.fluctuation_spectrum_planar_array)
                np.save(file_path_qvec,self.q2vec_array)

            if self.membrane_thickness_calculation=="on":
                file_path=os.path.join(self.output_folder_path,self.membrane_thickness_filename)+self.name_output_files
                np.save(file_path,self.membrane_thickness_array)


        if self.non_frame_iteration==True:
            file_path_energy=os.path.join(self.output_folder_path,self.energy_filename)+self.name_output_files
            file_path_energy_steps=os.path.join(self.output_folder_path,self.energy_steps_filename)+self.name_output_files
            np.save(file_path_energy,self.energy_array)
            np.save(file_path_energy_steps,self.energy_MCsteps_array)

            file_path_projected_area_ef=os.path.join(self.output_folder_path,self.projected_area_energy_filename)+self.name_output_files
            np.save(file_path_projected_area_ef,self.projected_area_array_energy_file)

    def save_data_pt(self):
        #Writing the output file
        file_path = os.path.join(self.output_folder_path, self.output_analysis_filename) + self.name_output_files
        self.write_output_file(file_path)
        #Saving the tempering moves
        file_path_tempering_moves=os.path.join(self.output_folder_path, self.tempering_moves_filename) + self.name_output_files
        np.save(file_path_tempering_moves,self.tempering_moves_array)
        #Saving the beta list
        file_path_beta_list=os.path.join(self.output_folder_path, self.beta_list_filename) + self.name_output_files
        np.save(file_path_beta_list,self.beta_list)

        if self.non_frame_iteration:
            # Save each energy array individually
            for i, energy in enumerate(self.energy_array):
                file_path_energy = os.path.join(self.output_folder_path, f"{self.energy_filename}_{i}") + self.name_output_files
                np.save(file_path_energy, energy)

            # Save each energy steps array individually
            for i, energy_steps in enumerate(self.energy_MCsteps_array):
                file_path_energy_steps = os.path.join(self.output_folder_path, f"{self.energy_steps_filename}_{i}") + self.name_output_files
                np.save(file_path_energy_steps, energy_steps)

            # Save each projected area energy file array individually
            for i, projected_area in enumerate(self.projected_area_array_energy_file):
                file_path_projected_area_ef = os.path.join(self.output_folder_path, f"{self.projected_area_energy_filename}_{i}") + self.name_output_files
                np.save(file_path_projected_area_ef, projected_area)


            



    def array_initialization_frames(self):

        if self.area_calculation=="on":
            self.area_array=np.zeros((self.nframes),float)
            self.projected_area_array=np.zeros((self.nframes),float)

        if self.inclusion_average_neighbours_calculation=="on":
            self.inclusion_average_neighbours_array=np.zeros((self.nframes),float)

        if self.inclusion_cluster_statistics_calculation=="on":
            #Since maximum cluster size is the number of inclusions:
            self.inclusion_cluster_statistics_array=np.zeros((self.nframes,self.ninclusion),float)

        if self.fluctuation_spectrum_planar_calculation=="on":
            self.q2vec_array=np.zeros((self.nframes,self.by+1,2*self.bx+1),dtype=float)
            if self.fluctuation_spectrum_no_inclusions=="off":
                self.fluctuation_spectrum_planar_array=np.zeros((self.nframes,self.by+1,2*self.bx+1,2,2),dtype=complex)
            else:
                self.fluctuation_spectrum_planar_array=np.zeros((self.nframes,self.by+1,2*self.bx+1),dtype=complex)

        if self.membrane_thickness_calculation=="on":
            self.membrane_thickness_array=np.zeros((self.nframes),float)
    

    def check_frame_existence(self,index):
        """Checks the existance of frame with a given index, and gives the path to the frame with index
        index."""

        # Construct the filename based on the index
        filename =self.name_TrjTSI_files+str(index)+".tsi"
        file_path = os.path.join(self.directory_path,self.name_TrjTSI_folder,filename)

        return os.path.exists(file_path),file_path

    


    def perform_calculation(self,file_path,i):
        
        frame_object=frame.frame(file_path)
        i=i+self.stupid_counter
        """
        self.inclusion_average_neighbours_calculation = None
        self.inclusion_cluster_statistics_calculation = None
        self.fluctuation_spectrum_planar_calculation = None
        self.bx=None
        self.by=None
        self.membrane_thickness_calculation = None"""

        try:
            if self.area_calculation=="on":
                frame_object.area_calculation()
                self.area_array[i]=frame_object.area
                self.projected_area_array[i]=frame_object.projected_area

                pass

            if self.inclusion_average_neighbours_calculation=="on":
                frame_object.inclusion_connectivity()
                self.inclusion_average_neighbours_array[i]=frame_object.inclusion_average_neighbours
                pass

            if self.inclusion_average_neighbours_calculation=="on" and self.inclusion_cluster_statistics_calculation=="on":
                frame_object.inclusion_cluster()
                self.inclusion_cluster_statistics_array[i,frame_object.inclusion_cluster_sizes]=frame_object.inclusion_cluster_frequency

            if self.inclusion_average_neighbours_calculation=="off" and self.inclusion_cluster_statistics_calculation=="on":
                frame_object.inclusion_connectivity()
                frame_object.inclusion_cluster()
                self.inclusion_cluster_statistics_array[i,frame_object.inclusion_cluster_sizes]=frame_object.inclusion_cluster_frequency

            if self.fluctuation_spectrum_planar_calculation=="on":
                if self.fluctuation_spectrum_no_inclusions=="off":
                    self.fluctuation_spectrum_planar_array[i,:,:,:,:],self.q2vec_array[i,:,:]=frame_object.ft_matrix(self.bx,self.by)
                else:
                    self.fluctuation_spectrum_planar_array[i,:,:],self.q2vec_array[i,:,:]=frame_object.ft_matrix_no_inc(self.bx,self.by)

            if self.membrane_thickness_calculation=="on":
                self.membrane_thickness_array[i]=frame_object.thickness

        except Exception as e:
            print("An error occurred while performing calculation for frame {}: {}".format(i, e))
            """Now I want this code to be modified. I would like to eliminate to pop the elemnt i from each
            of the lists after the conditionals. Also, I would like to modify the self.frame_num_list and 
            also eliminate its element i. """

            if self.area_calculation == "on":
                self.area_array = np.delete(self.area_array, i, axis=0)
                self.projected_area_array = np.delete(self.projected_area_array, i, axis=0)

            if self.inclusion_average_neighbours_calculation == "on":
                self.inclusion_average_neighbours_array = np.delete(self.inclusion_average_neighbours_array, i, axis=0)

            if self.inclusion_average_neighbours_calculation == "on" and self.inclusion_cluster_statistics_calculation == "on":
                self.inclusion_cluster_statistics_array = np.delete(self.inclusion_cluster_statistics_array, i, axis=0)

            if self.inclusion_average_neighbours_calculation == "off" and self.inclusion_cluster_statistics_calculation == "on":
                self.inclusion_cluster_statistics_array = np.delete(self.inclusion_cluster_statistics_array, i, axis=0)

            if self.fluctuation_spectrum_planar_calculation == "on":
                self.fluctuation_spectrum_planar_array = np.delete(self.fluctuation_spectrum_planar_array, i, axis=0)
                self.q2vec_array = np.delete(self.q2vec_array, i, axis=0)

            if self.membrane_thickness_calculation == "on":
                self.membrane_thickness_array = np.delete(self.membrane_thickness_array, i, axis=0)

            self.frame_num_list=np.delete(self.frame_num_list,i)

            self.stupid_counter-=1
        
            



        #perform calculations
    
    def iteration_frames(self):

        """Iterates over the desired frames and obtains quantities specified
        in the input file."""

        for i in range(0,len(self.frame_path_list)):

            #start_time = time.time()

            self.perform_calculation(self.frame_path_list[i],i)


            #end_time = time.time()

            #elapsed_time = end_time - start_time
            #print("Percentatge",i/(self.frame_num_list[-1]),"Elapsed time:", elapsed_time, "seconds")

    def read_temperature_folders(self,name):
        temp_folders = []
        for item in os.listdir(self.directory_path):
            if os.path.isdir(os.path.join(self.directory_path, item)) and item.startswith(name):
                temp_folders.append(item)
        return temp_folders
    
    def read_energy_files_pt(self, name):
        energy_files = []
        for item in os.listdir(self.directory_path):
            if os.path.isfile(os.path.join(self.directory_path, item)) and item.startswith(name) and item.endswith('.xvg'):
                energy_files.append(item)
        return energy_files
    
    
    def sort_and_order_lists(self,name_TrjTSI_folder,name_energy_files):
        # Step 1: Read temperature folders and energy files
        temp_folders = self.read_temperature_folders(name_TrjTSI_folder)
        energy_files = self.read_energy_files_pt(name_energy_files)

        print(temp_folders,energy_files)
        # Step 2: Extract temperatures from temp_folders
        temperatures = []
        for folder in temp_folders:
            try:
                temperature = float(folder.split('_')[-1])
                temperatures.append(temperature)
            except ValueError:
                pass
        print(temperatures)
        # Step 3: Convert temperatures to a numpy array and sort it
        temperatures_array = np.array(temperatures)
        sorted_temperatures = np.sort(temperatures_array)

        # Step 4: Order temp_folders and energy_files based on sorted_temperatures
        ordered_temp_folders = []
        ordered_energy_files = []

        for temp in sorted_temperatures:
            temp_str = f"TrjTSI_temp_{temp}"
            for folder in temp_folders:
                if temp_str in folder:
                    ordered_temp_folders.append(folder)
                    break

        for temp in sorted_temperatures:
            temp_str = f"output_temp_{temp}-en.xvg"
            for file in energy_files:
                if temp_str in file:
                    ordered_energy_files.append(file)
                    break

        return ordered_temp_folders, ordered_energy_files, sorted_temperatures
    
    def read_tempering_moves(self,name_tempering_moves_file):
        file_path = os.path.join(self.directory_path, name_tempering_moves_file)
        
        tempering_moves_list = []
        with open(file_path, 'r') as file:
            for line in file:
                # Split the line by spaces and convert each element to an integer
                tempering_moves_list.append([int(x) for x in line.split()])
        
        # Convert the list of lists to a numpy array
        return np.array(tempering_moves_list)

