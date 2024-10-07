#if !defined(AFX_FluctuationSpectrum_H_334B21B8_C13C_2248_BF23_124095086233__INCLUDED_)
#define AFX_FluctuationSpectrum_H_334B21B8_C13C_2248_BF23_124095086233__INCLUDED_

#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
#include "Inclusion_Interaction_Map.h"

/*
 * @brief Energy calculation based on FreeDTS1.1 force field.
 *
 * This class is responsible for calculating the energy of a system using the FreeDTS1.1 force field.
 * It includes methods to compute various energy contributions such as single vertex energy, energy
 * due to interactions between two inclusions.

 *
 * @note This class inherits from the AbstractEnergy interface, providing a common interface for energy
 * calculation modules in the simulation framework.
 *
 * @author Weria Pezeshkian (weria.pezeshkian@gmail.com)
 * @copyright Weria Pezeshkian
 */
class State;
class FluctuationSpectrum {
public:
    FluctuationSpectrum(State* pState, int Nx, int Ny);
	 ~FluctuationSpectrum();

public:
    void CalculateSpectrum();
private:
    void GenerateVectorOrder();
    void GenerateZeroAndNonZeroVectorOrder();
    double FourierTransform(std::vector<std::vector<double>> qvector);
    void AverageHeight();
    


private:
    State* m_pState;
    int m_Nx,m_Ny;
    std::vector<std::vector<int>> m_VectorOrder;
    int m_SpectrumSize;
    std::vector<std::vector<std::vector<int>>> m_MatrixOrder;
    double m_AverageHeight;




};


#endif
