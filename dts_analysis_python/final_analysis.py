import numpy as np
import os
import pickle

def count_lines(filename):
    with open(filename, 'r') as file:
        line_count = 0
        for line in file:
            line_count += 1
    return line_count

# Function to check if file is readable
def is_file_readable(file_path):
    return os.access(file_path, os.R_OK)

#Setting the grounds

filename = 'directory.txt'  # Change this to directory file path
nsimulations = count_lines(filename)
nparameters= int(nsimulations/4)


inclusion_average_neighbours=[]
thickness=[]
frame_steps=[]
area=[]
q2vec=[]
energy=[]
energy_MCsteps=[]
projected_area=[]
fluctuation_spectrum=[]
inclusion_cluster_statistics=[]

inclusion_parameters=[]

"""
inclusion_C0=np.array([])
inclusion_density=np.array([])
inclusion_kappa=np.array([])
inclusion_kappa_g=np.array([])"""

area_filename='area.npy'
projected_area_filename='projected_area.npy'
membrane_thickness_filename='thickness.npy'
fluctuation_spectrum_filename='fluctuation_spectrum.npy'
energy_filename='energy.npy'
energy_steps_filename='energy_MCsteps.npy'
inclusion_average_neighbours_filename='inclusion_average_neighbours.npy'
inclusion_cluster_filename='inclusion_cluster_statistics.npy'
qvec_filename='q2vec.npy'
frame_steps_filename='frame_steps.npy'
output_analysis_filename='analysis_output.txt'

not_working_arrays=np.array([],dtype=int)


# Function to read .npy files from a directory
def read_npy_files(directory):
    npy_files = []
    for file in os.listdir(directory):
        if file.endswith('.npy'):
            npy_files.append(os.path.join(directory, file))
    return npy_files

# Read directories from directory.txt
with open('directory.txt', 'r') as file:
    directories = file.readlines()

# Iterate over each group of 4 directories
for i in range(nparameters):
    # Determine the start and end indices for the current group of 4 directories
    start_index = i * 4
    end_index = min((i + 1) * 4, len(directories))
    counter=0
    check=0
    check2=0
    k=0
    thickness_s=[]
    frame_steps_s=[]
    area_s=[]
    q2vec_s=[]
    energy_MCsteps_s=[]
    energy_s=[]
    projected_area_s=[]
    fluctuation_spectrum_s=[]
    inclusion_cluster_statistics_s=[]
    inclusion_average_neighbours_s=[]
    print("Parameter percentage:", i/nparameters)
    # Process each directory within the current group
    for j in range(start_index, end_index):
        directory = directories[j].strip()  # Remove leading/trailing whitespace and newline characters
        analysis_output_dir = os.path.join(directory, 'analysis_output')



        # Check if the analysis output directory exists
        if os.path.exists(analysis_output_dir):


            analysis_output_file = os.path.join(analysis_output_dir, output_analysis_filename)

            # Check if the analysis output file exists
            if os.path.exists(analysis_output_file):
                #check=1

                if counter==0:
                    with open(analysis_output_file, 'r') as f:
                        lines = f.readlines()
                        # Extract and append values to corresponding variables
                        ik = float(lines[4].split('=')[1].strip())
                        ikg = float(lines[5].split('=')[1].strip())
                        iC0 = float(lines[6].split('=')[1].strip())
                        id = float(lines[7].split('=')[1].strip())

                        it=[id,iC0,ik,ikg]

                        inclusion_parameters.append(it)
                        counter+=1
                

            else:
                print(i,f"Analysis output file {analysis_output_file} does not exist.")
                not_working_arrays=np.append(not_working_arrays,j)



            """Next thing I would like to do here is to read, from analysis_output_dir, all the numpy files,
            which are identified with the filenames. To do so, first I construct the file to the path to each numpy file,
            which is given by os.path.join(analysis_output_dir, something_filename), and check its existence. 
            Then, I want to append, to each something_s, the corresponding information stored in something_filename using 
             numpy load(). """
            
             # Read .npy files in analysis_output_dir
            npy_files = read_npy_files(analysis_output_dir)
            if len(npy_files)!=0:
                check+=1
                """
                k=0
                for npy_file in npy_files:
                    filename = os.path.basename(npy_file)
                    if filename==inclusion_average_neighbours_filename:
                        k+=1
                        continue
                
                if k==0:
                    print("inclusion_average_neighbours does not exist in:", i,analysis_output_dir)"""
                    


                for npy_file in npy_files:
                    try:
                    # Load the data from npy_file
                        data = np.load(npy_file)

                        # Append the data to corresponding variables based on filename
                        filename = os.path.basename(npy_file)
                        if filename == area_filename:
                            area_s.append(data)
                        elif filename == projected_area_filename:
                            #data_row = data.reshape(1, -1)
                            #projected_area_s = np.append(projected_area_s, data_row,axis=0)
                            projected_area_s.append(data)
                            
                        elif filename == membrane_thickness_filename:
                            
                            thickness_s.append(data)
                        elif filename == fluctuation_spectrum_filename:
                            fluctuation_spectrum_s.append(data)
                        elif filename == energy_filename:
                            energy_s.append(data)
                        elif filename == energy_steps_filename:
                            energy_MCsteps_s.append(data)
                        elif filename == inclusion_average_neighbours_filename:
                            inclusion_average_neighbours_s.append(data)
                        elif filename == inclusion_cluster_filename:
                            inclusion_cluster_statistics_s.append(data)
                        elif filename == qvec_filename:
                            q2vec_s.append(data)
                        elif filename == frame_steps_filename:
                            frame_steps_s.append(data)
                        else:
                            print(f"Unknown filename: {filename}")
                    
                    except OSError as e:
                        print("it happeneded")
                        data = 0
                        # Append the data to corresponding variables based on filename
                        filename = os.path.basename(npy_file)

                        if filename == area_filename:
                            area_s.append(data)
                        elif filename == projected_area_filename:
                            #data_row = data.reshape(1, -1)
                            #projected_area_s = np.append(projected_area_s, data_row,axis=0)
                            projected_area_s.append(data)
                            
                        elif filename == membrane_thickness_filename:
                            
                            thickness_s.append(data)
                        elif filename == fluctuation_spectrum_filename:
                            fluctuation_spectrum_s.append(data)
                        elif filename == energy_filename:
                            energy_s.append(data)
                        elif filename == energy_steps_filename:
                            energy_MCsteps_s.append(data)
                        elif filename == inclusion_average_neighbours_filename:
                            inclusion_average_neighbours_s.append(data)
                        elif filename == inclusion_cluster_filename:
                            inclusion_cluster_statistics_s.append(data)
                        elif filename == qvec_filename:
                            q2vec_s.append(data)
                        elif filename == frame_steps_filename:
                            frame_steps_s.append(data)
                        else:
                            print(f"Unknown filename: {filename}")
                        
                        """And then, you just do not append :)"""
                    
                    
            else:
                print(i,f"Numpy files at {analysis_output_dir} do not exist.")



        else:
            print(f"Analysis output directory {analysis_output_dir} does not exist.")
    
    if check>0:
        inclusion_average_neighbours.append(inclusion_average_neighbours_s)
        thickness.append(thickness_s)
        frame_steps.append(frame_steps_s)
        area.append(area_s)
        q2vec.append(q2vec_s)
        energy_MCsteps.append(energy_MCsteps_s)
        energy.append(energy_s)
        projected_area.append(projected_area_s)
        fluctuation_spectrum.append(fluctuation_spectrum_s)
        inclusion_cluster_statistics.append(inclusion_cluster_statistics_s)
        

print(len(projected_area))

print(len(fluctuation_spectrum[2]))

print(fluctuation_spectrum[2][0].shape)

inclusion_parameters=np.array(inclusion_parameters)
print(inclusion_parameters)

np.save('inclusion_parameters',inclusion_parameters)
np.save('not_working_arrays',not_working_arrays)
print(not_working_arrays)

# Assuming 'data_list' is your Python list containing NumPy arrays
with open('inclusion_average_neighbours.pkl', 'wb') as f:
    pickle.dump(inclusion_average_neighbours, f)

with open('thickness.pkl', 'wb') as f:
    pickle.dump(thickness, f)

with open('fluctuation_spectrum.pkl', 'wb') as f:
    pickle.dump(fluctuation_spectrum, f)

with open('projected_area.pkl', 'wb') as f:
    pickle.dump(projected_area, f)



with open('area.pkl', 'wb') as f:
    pickle.dump(area, f)

with open('inclusion_cluster_statistics.pkl', 'wb') as f:
    pickle.dump(inclusion_cluster_statistics, f)

with open('energy_MCsteps.pkl', 'wb') as f:
    pickle.dump(energy_MCsteps, f)

with open('energy.pkl', 'wb') as f:
    pickle.dump(energy, f)

with open('q2vec.pkl', 'wb') as f:
    pickle.dump(q2vec, f)

with open('frame_steps.pkl', 'wb') as f:
    pickle.dump(frame_steps, f)


