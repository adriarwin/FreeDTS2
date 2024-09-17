
import numpy as np
import os
from scipy.sparse.csgraph import connected_components


class frame:
    def __init__(self, file_path):
        

        #INITIALIZATION VALUES

        #box size
        self.Lx=None
        self.Ly=None
        self.Lz=None

        #vertex,triangle and inclusion numbers
        self.nvertex=None
        self.ntriangle=None
        self.ninclusion=None

        #vertex,triangle lists
        self.vertex=None
        self.inclusion_type=None

        #List of inclusions and their vertex
        self.inclusion_vertex=None

        #List of vertexes per inclusion
        self.vertex_occupancy=None
        self.vertex_neighbours=None
        self.triangle=None

        #average thickness

        self.thickness=None

        #METHODS VARIABLES

        #Inclusion connectivity matrix
        self.inclusion_connectivity_matrix=None

        #Inclusion-inclusion average neighbours
        self.inclusion_average_neighbours=None

        #Inclusion sizes of different clusters
        self.inclusion_cluster_sizes=None

        #Frequency of each inclusion cluster
        self.inclusion_cluster_frequency=None

        #Area of the surface
        self.area=None
        
        
        



        self.load_data(file_path)


    def load_data(self,file_path):
        try:
            with open(file_path, 'r') as file:
                # Read the entire file
                file_contents = file.read().splitlines()
                #version variable
                version=file_contents[0]
                #determining the box size
                box=file_contents[1].split()

                self.Lx=float(box[1])
                self.Ly=float(box[2])
                self.Lz=float(box[3])

                self.projected_area=self.Lx*self.Ly

                #Determining the number of vertexes,triangles, and inclusions
                self.nvertex=int(file_contents[2].split()[1])
                self.ntriangle=int(file_contents[3+self.nvertex].split()[1])
                self.ninclusion=int(file_contents[4+self.nvertex+self.ntriangle].split()[1])
                #Vertex array
                self.vertex=np.zeros((self.nvertex,3),dtype=float)
                #triangle array
                self.triangle=np.zeros((self.ntriangle,3),dtype=int)
                #neighbours array
                # Create the empty array
                self.vertex_neighbours = np.zeros((self.nvertex,self.nvertex), dtype=int)
                #what vertex is the inclusion occupaying
                self.inclusion_vertex=np.zeros((self.ninclusion),dtype=int)
                #type of inclusion
                self.inclusion_type=np.zeros((self.ninclusion),dtype=int)
                #Vertex occupancy (what inclusion is on the vertex. -1 if the vertex is empty)
                self.vertex_occupancy=-1*np.ones((self.nvertex),dtype=int)

                #thickness
                

                #Filling data in for the vertex
                for i in range(0,self.nvertex):
                    coordinate=file_contents[i+3].split()
                    self.vertex[i]=[float(coordinate[1]),float(coordinate[2]),float(coordinate[3])]
                
                #eliminating PBC if they were crossed
                self.thickness=np.max(self.vertex[:,2])-np.min(self.vertex[:,2])
                #This part needs to be rewritten. What about if you cross PBC in the 
                #X and Y direction? For now we can leave it like that, sine all the
                #analysis is to be done in other places. But for a vesicle, this needs
                #to be changed!

                if self.thickness>0.9*self.Lz:
                    # Find elements smaller than 50% Lz
                    mask = self.vertex[:,2] < 0.5*self.Lz
                    # Add 100 to elements smaller than 50
                    self.vertex[:,2][mask] += self.Lz

                self.thickness=np.max(self.vertex[:,2])-np.min(self.vertex[:,2])
                #Filling data in for inclusions
                for i in range(0,self.ninclusion):
                    coordinate=file_contents[i+self.nvertex+self.ntriangle+5].split()
                    self.inclusion_vertex[i]=int(coordinate[2])
                    self.inclusion_type[i]=int(coordinate[1])
                    self.vertex_occupancy[self.inclusion_vertex[i]]=i
                #filling data in for triangles!
                for i in range(0,self.ntriangle):
                    coordinate=file_contents[i+self.nvertex+4].split()
                    self.triangle[i]=[int(coordinate[1]),int(coordinate[2]),int(coordinate[3])]
                    #determining the neighbours of each vertex
                    self.vertex_neighbours[self.triangle[i][0]][self.triangle[i][1]]=1
                    self.vertex_neighbours[self.triangle[i][0]][self.triangle[i][2]]=1
                    self.vertex_neighbours[self.triangle[i][1]][self.triangle[i][2]]=1
                    self.vertex_neighbours[self.triangle[i][1]][self.triangle[i][0]]=1
                    self.vertex_neighbours[self.triangle[i][2]][self.triangle[i][0]]=1
                    self.vertex_neighbours[self.triangle[i][2]][self.triangle[i][1]]=1

        except IOError:
            print("The file '{}' was not found.".format(file_path))
            return None
        except Exception as e:
            print("An error occurred: {}".format(e))
            return None
    

    

    def inclusion_connectivity(self):
        """This function calculates inclusion connectivity matrix, of dimensions ninclusion X ninclusion
        where its elements i and j are different than zero when inclusion i and j are nearest neighbours."""

        #defining connectivity matrix
        self.inclusion_connectivity_matrix=np.zeros((self.ninclusion,self.ninclusion),dtype=int)
        #defining number of neighbours
        self.inclusion_average_neighbours=0.0

        for i in range(0,self.ninclusion):
            #First, a vector with inclusion neighbours when it is not zero
            inclusion_indices=np.where(self.vertex_neighbours[self.inclusion_vertex[i]][self.inclusion_vertex]*self.vertex_occupancy[self.inclusion_vertex]!=0)[0]
            size_inclusion_indices=inclusion_indices.size
            if size_inclusion_indices!=0:
                self.inclusion_connectivity_matrix[i,inclusion_indices]=1
                self.inclusion_average_neighbours+=size_inclusion_indices

        self.inclusion_average_neighbours=self.inclusion_average_neighbours/(2*self.ninclusion)
    
    def inclusion_cluster(self):
        """This function calculates the number of connected components of a graph. Connectivity matrix
        must be calculated first. It returns inclusion cluster sizes, and the number of times each cluster
        size appears. """

        #Obtaining list of labels for each inclusion (to what connected graph each inclusion belongs )
        result=connected_components(self.inclusion_connectivity_matrix, directed=False, return_labels=True)
        #Using unique function, we find size of each connected graph
        unique_elements, counts = np.unique(np.array(result[1]), return_counts=True)
        #counts is a 1D array that contains sizes of all connected graphs. By applaying another unique,
        #to counts, we obtain number of clusters of a given size
        self.inclusion_cluster_sizes, self.inclusion_cluster_frequency = np.unique(np.array(counts), return_counts=True)

    def area_calculation(self):
        """This function calculates the area and proejcted area of our surface. It defines,
        for each triangle, two vectors P0P1 and P0P2 and eliminates PBC effects before calculating
        its area. The sum of the areas of all the triangles gives the toal area."""

        #First, we calculate P1P2 and P1P3 vectors
        P0P1=self.vertex[self.triangle[:,1]]-self.vertex[self.triangle[:,0]]
        P0P2=self.vertex[self.triangle[:,2]]-self.vertex[self.triangle[:,0]]

        #Removal of PBC in the X and Y directions for the triangles that cross the
        #edges of the box

        #X direction, positive vectors
        maskxP0P1 = P0P1[:,0] > 0.5*self.Lx
        maskxP0P2 = P0P2[:,0] > 0.5*self.Lx

        P0P1[:,0]-=maskxP0P1*self.Lx
        P0P2[:,0]-=maskxP0P2*self.Lx

        #negative vectors
        maskxP0P1 = P0P1[:,0] < -0.5*self.Lx
        maskxP0P2 = P0P2[:,0] < -0.5*self.Lx

        P0P1[:,0]+=maskxP0P1*self.Lx
        P0P2[:,0]+=maskxP0P2*self.Lx

        #Y direction, positive vectors
        maskxP0P1 = P0P1[:,1] > 0.5*self.Ly
        maskxP0P2 = P0P2[:,1] > 0.5*self.Ly

        P0P1[:,1]-=maskxP0P1*self.Ly
        P0P2[:,1]-=maskxP0P2*self.Ly

        #negative vectors
        maskxP0P1 = P0P1[:,1] < -0.5*self.Ly
        maskxP0P2 = P0P2[:,1] < -0.5*self.Ly

        P0P1[:,1]+=maskxP0P1*self.Ly
        P0P2[:,1]+=maskxP0P2*self.Ly

        uvec=np.cross(P0P1,P0P2)

        self.area=np.sum(np.sqrt(uvec[:,0]**2 + uvec[:,1]**2 + uvec[:,2]**2))/2


        self.projected_area=self.Lx*self.Ly

    def ft_height(self,bx,by):
        """Fourier transfrom of array A, where A[i,:]=r_i and r_i=(x,y,z). 
        The fourier transfrom is applied on the height, defined as h(x,y)=z. Lx and
        Ly define the size of the simulation box and bx and by are the boundaries of 
        q, defined as Lx/bx and Ly/by. """
        z_vec=self.vertex[:,2]
        r_vec=self.vertex[:,0:2]
        minh=np.average(z_vec)
        fresult=np.zeros((by+1,2*bx+1),dtype=complex)
        for i in range(0,by+1):
            for j in range(0,2*bx+1):
                n=i
                m=-bx + j
                q_vec=np.array([m*2*np.pi/self.Lx,n*2*np.pi/self.Ly])
                ax=np.dot(r_vec,q_vec)
                fresult[i,j]=np.sum(np.multiply(np.exp(-1j*np.dot(r_vec,q_vec)),(z_vec-minh)))

        return fresult/(np.sqrt(self.Lx*self.Ly))
    
    def ft_rho(self,bx,by):
        """Fourier transfrom of array A, where A[i,:]=r_i and r_i=(x,y,z). 
        The fourier transfrom is applied on the height, defined as h(x,y)=z. Lx and
        Ly define the size of the simulation box and bx and by are the boundaries of 
        q, defined as Lx/bx and Ly/by. """
        z_vec=self.vertex[:,2]
        r_vec=self.vertex[:,0:2]
        fresult=np.zeros((by+1,2*bx+1),dtype=complex)
        rho0=self.ninclusion/self.nvertex
        vertex_occupancy= self.vertex_occupancy > -0.5

        for i in range(0,by+1):
            for j in range(0,2*bx+1):
                n=i
                m=-bx + j
                q_vec=np.array([m*2*np.pi/self.Lx,n*2*np.pi/self.Ly])
                fresult[i,j]=np.sum(np.multiply(np.exp(-1j*np.dot(r_vec,q_vec)), (vertex_occupancy-rho0)))

        #return (self.nvertex/(self.Lx*self.Ly))*fresult/(np.sqrt(self.Lx*self.Ly))
        return fresult

    def qvec2(self,bx,by):

        fresult=np.zeros((by+1,2*bx+1),dtype=float)

        for i in range(0,by+1):
            for j in range(0,2*bx+1):
                n=i
                m=-bx + j
                fresult[i,j]=(m*2*np.pi/(self.Lx))**2 + (n*2*np.pi/(self.Ly))**2
        
        return fresult
    
    def ft_matrix(self,bx,by):

        ft_h=self.ft_height(bx,by)
        ft_r=self.ft_rho(bx,by)

        spectrum=np.zeros((by+1,2*bx+1,2,2),dtype=complex)

        spectrum[:,:,0,0]=ft_h*np.conjugate(ft_h)
        spectrum[:,:,1,1]=ft_r*np.conjugate(ft_r)
        spectrum[:,:,0,1]=ft_h*np.conjugate(ft_r)
        spectrum[:,:,1,0]=ft_r*np.conjugate(ft_h)


        q2vec=self.qvec2(bx,by)

        return spectrum,q2vec
    
    def ft_matrix_no_inc(self,bx,by):

        ft_h=self.ft_height(bx,by)

        spectrum=np.zeros((by+1,2*bx+1),dtype=complex)

        spectrum=ft_h*np.conjugate(ft_h)

        q2vec=self.qvec2(bx,by)

        return spectrum,q2vec

        




        

        
                   


"""
directory_path1=r'/home/adriarwin/analysis/output3102.tsi'

import time

start_time = time.time()

a=frame(directory_path1)
a.inclusion_connectivity()
b=a.inclusion_connectivity_matrix
a.inclusion_cluster()
a.area_calculation()
b=a.ft_matrix(6,6)

end_time = time.time()


elapsed_time = end_time - start_time

print("Elapsed time:",elapsed_time, "seconds")


print("Time for 5000 frames:",elapsed_time*5000/60, "minutes")"""




