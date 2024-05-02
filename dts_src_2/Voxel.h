#if !defined(TEM_Voxel_H_INCLUDED_)
#define TEM_Voxel_H_INCLUDED_

#include "SimDef.h"
#include "Voxelization.h"
/*
    Author: Weria Pezeshkian (weria.pezeshkian@gmail.com)
    Updated: 2024

    Description:
    The Voxel class represents a single voxel within a 3D grid used for voxelizing a box-shaped region. It inherits from the Voxelization class.

    Voxel Characteristics:
    - Each voxel contains a list of objects (vertices/points) that fall within its bounds.
    - The voxel grid is dynamic, adjusting its size and contents during simulation time.

    Dependencies:
    - SimDef.h: Contains definitions and constants used in the simulation.
    - Voxelization.h: Declaration of the Voxelization class, which handles the voxelization process.

    Usage Example:
    // Include necessary headers
    #include "Voxel.h"
    #include "SomeObjectType.h"

    // Create an instance of Voxel with appropriate parameters
    Voxel<SomeObjectType> voxel();

    // Access voxel properties and functionalities
    voxel.AddtoContentList();

    Class Template:
    - Voxel is templated on Type, allowing it to work with different types of objects.

    Important Notes:
    - The Voxel class assumes ownership of the objects added to its content list. Ensure proper memory management when removing objects.
*/
template<typename Type>
class Voxel : public Voxelization<Type> {
public:
    Voxel(int id, int n, int m, int k) {
        m_X_index = n;                                          //   -------
        m_Y_index = m;                                          //     | |
        m_Z_index = k;                                          //   -------
        m_ID=id;
        m_NeighbouringVoxel[1][1][1] = this;
    }
    ~Voxel(){
        for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
        for (int k = 0; k < 3; ++k) {
            delete m_NeighbouringVoxel[i][j][k];
        }}}
        m_List.clear();
    }// ~Voxel()

    inline const int          GetID()       const        { return m_ID; }      // voxel id
    inline int                GetXIndex()                { return m_X_index; } // x index,  m_X_index*m_Vox_Lx 
    inline int                GetYIndex()                { return m_Y_index; } //
    inline int                GetZIndex()                { return m_Z_index; } //
    inline std::vector<Type*> GetContentObjects()        { return m_List; }    //all the objects, e.g., vertices in the voxel

public:
    void AddtoContentList(Type * p_obj){ // adding a new vertex/point, when the vertex enters the domain
        m_List.push_back(p_obj);
        return;
    }
    void ClearList(){ // adding a new vertex/point, when the vertex enters the domain
        m_List.clear();
        return;
    }
    bool SetANeighbourCell(int i, int j, int k, Voxel *vox){// setting the nb voxels 26+1 voxels
        if(i > 1 || i < -1 || j > 1 || j < -1 || k > 1 || k < -1) {
            std::cout << " ---> error: Invalid indices. Indices must be within range [-1, 1]." << std::endl;
            return false;
        }

        m_NeighbouringVoxel[i + 1][j + 1][k + 1] = vox;
        return true;
    }
    void RemoveObjectFromContentList(Type * p_obj){// removing a  vertex, when the vertex exist the domain
        m_List.erase(std::remove(m_List.begin(), m_List.end(), p_obj), m_List.end());
        return;
    }
    Voxel * GetANeighbourCell(int i, int j, int k){// accessing a specific nighbouring cell
        if(i==0 && j==0 && k==0){
            return this;
        }
        else if(i>1 || i<-1 || j>1 || j<-1 || k>1 || k<-1  )
            std::cout<<" ---> error: such indices are not permitted "<<std::endl;
        
        return m_NeighbouringVoxel[i+1][j+1][k+1];
    }
private:
  int m_ID;
  std::vector <Type *> m_List;              // List of the vertices/points that the cell contains/ Content
    int m_X_index;                          // indices of the cell withing the box;
    int m_Y_index;
    int m_Z_index;
    Voxel *m_NeighbouringVoxel[3][3][3];
    


};

#endif
