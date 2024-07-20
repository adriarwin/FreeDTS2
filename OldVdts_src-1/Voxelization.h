#ifndef VOXELIZATION_H_INCLUDED
#define VOXELIZATION_H_INCLUDED

#include "SimDef.h"
#include "Voxel.h"
#include "Vec3D.h"
/*
 Weria Pezeshkian@ 2024,
 This class divides the box into a grid and assigns each object, such as vertices or points, to a voxel based on its position. In the case of a dynamic box, the voxel size decreases while the number of voxels remains constant until the Voxelize function is executed again. Therefore, if the voxel size becomes smaller than a certain threshold, lmin, the Voxelize function should be rerun.
 */

template<typename Type>
class Voxelization
{
public:
Voxelization(Vec3D pBox) { // Use Type here
        // Initialize the box object
        m_pBox = pBox;
        m_Nx = 2;
        m_Ny = 2;
        m_Nz = 2;

}
~Voxelization() {// to clear the memory allocated for this
        if (m_AllVoxel != nullptr) {
            for (int i = 0; i < m_Nx; ++i) {
            for (int j = 0; j < m_Ny; ++j) {
            for (int k = 0; k < m_Nz; ++k) {
                delete m_AllVoxel[i][j][k];
            }
                delete[] m_AllVoxel[i][j];
            }
                delete[] m_AllVoxel[i];
            }
            delete[] m_AllVoxel;
            m_AllVoxel = nullptr;
        }
}

inline int GetXSideVoxel() const { return (*m_pBox)(0)/double(m_Nx); }
inline int GetYSideVoxel() const { return (*m_pBox)(1)/double(m_Ny); }
inline int GetZSideVoxel() const { return (*m_pBox)(2)/double(m_Nz); }
inline int GetXVoxelNumber() const { return m_Nx; }
inline int GetYVoxelNumber() const { return m_Ny; }
inline int GetZVoxelNumber() const { return m_Nz; }



//----> Important and detailed function
bool Voxelize(std::vector<Type *> all_pObjects) {
        
        //--->find the appropriate voxel size and number of voxels
        m_Nx=int((*m_pBox)(0)/m_Lx);
        m_Ny=int((*m_pBox)(1)/m_Ly);
        m_Nz=int((*m_pBox)(2)/m_Lz);
        m_Lx = (*m_pBox)(0)/double(m_Nx);
        m_Ly = (*m_pBox)(1)/double(m_Ny);
        m_Lz = (*m_pBox)(2)/double(m_Nz);
//---> make m_AllVoxel empty
        if (m_AllVoxel != nullptr) {
            for (int i = 0; i < m_Nx; ++i) {
            for (int j = 0; j < m_Ny; ++j) {
            for (int k = 0; k < m_Nz; ++k) {
                delete m_AllVoxel[i][j][k];
            }
                delete[] m_AllVoxel[i][j];
            }
                delete[] m_AllVoxel[i];
            }
            delete[] m_AllVoxel;
            m_AllVoxel = nullptr;
        }
        //----> making m_AllVoxel an m_AllVoxel[m_Nx][m_Ny][m_Nz] and create all the cells
        int voxel_id = 0;
        m_AllVoxel = new Voxel<Type>***[m_Nx];
        for (int i = 0; i < m_Nx; ++i) {
            m_AllVoxel[i] = new Voxel<Type>**[m_Ny];
            for (int j = 0; j < m_Ny; ++j) {
                m_AllVoxel[i][j] = new Voxel<Type>*[m_Nz];
                for (int k = 0; k < m_Nz; ++k) {
                    m_AllVoxel[i][j][k] = new Voxel<Type>(voxel_id,i,j,k);
                    voxel_id++;
                }
            }
        }
        //--> Update voxel neighbours
        //----> making m_AllVoxel an m_AllVoxel[m_Nx][m_Ny][m_Nz] and create all the cells
        for (int i = 0; i < m_Nx; ++i) {
            for (int j = 0; j < m_Ny; ++j) {
                for (int k = 0; k < m_Nz; ++k) {
                    for (int n = -1; n <=1; ++n)
                    for (int m = -1; m <= 1; ++m)
                    for (int p = -1; p <= 1; ++p)
                    (m_AllVoxel[i][j][k])->SetANeighbourCell(n,m,p, m_AllVoxel[(i+n)%m_Nx][(j+m)%m_Ny][(k+p)%m_Nz]);
                }
            }
        }
        //--->allocate the objects to the voxels
        for (typename std::vector<Type *>::iterator it = all_pObjects.begin(); it != all_pObjects.end(); ++it) {
            double x = (*it)->GetXPos(); // Retrieve X position of each object
            double y = (*it)->GetYPos(); // Retrieve Y position of each object
            double z = (*it)->GetZPos(); // Retrieve Z position of each object
            // Process or allocate objects to voxels based on their Z position
            int nx = (*it)->GetXPos()/m_Lx;
            int ny = (*it)->GetYPos()/m_Ly;
            int nz = (*it)->GetZPos()/m_Lz;
            Voxel<Type> * pVoxel = m_AllVoxel[nx][ny][nz];
            (m_AllVoxel[nx][ny][nz])->AddtoContentList(*it);
            (*it)->UpdateVoxel(pVoxel);
        }

        return true; // Placeholder
}

private:
    int m_Nx; // Number of the voxels in the x direction
    int m_Ny; // Number of the voxels in the y direction
    int m_Nz; // Number of the voxels in the z direction
    double m_Lx; // voxel length in the x direction
    double m_Ly; // voxel length in the y direction
    double m_Lz; // voxel length in the z direction

    Vec3D *m_pBox;
    Voxel<Type> *m_AllVoxel[2][2][2]; 

public:
    // Add other member functions here
};

#endif // VOXELIZATION_H_INCLUDED
