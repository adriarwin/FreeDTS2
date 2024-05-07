#if !defined(AFX_CouplingtoFixedGlobalCurvature_H_334B21B8_C13C_2248_BF23_124095086255__INCLUDED_)
#define AFX_CouplingtoFixedGlobalCurvature_H_334B21B8_C13C_2248_BF23_124095086255__INCLUDED_
/*
 * CouplingtoFixedGlobalCurvature:
 * Author: Weria Pezeshkian (weria.pezeshkian@gmail.com)
 * Copyright (c) Weria Pezeshkian
 *
 * This object affects the simulation results in every step. It is a global coupling that penalizes the deviations
 * of total mean curvature from a target value. Additionally, it can model the energy associated with the difference
 * in the inner and outer monolayer area. The energy coupling models two effects:
 *
 * 1. Energy related to mean curvature deviation:
 *    E = k / (2A) * (∑(2h - gC0)av)^2
 *    Where:
 *      - E is the energy
 *      - k is the energy coupling constant (will be divided by 2 in the constructor)
 *      - A is the area
 *      - h is the total mean curvature
 *      - gC0 is the target mean curvature
 *
 * 2. Energy related to area difference:
 *    E = k / (2Ah^2) * (ΔA - ΔA0)^2
 *    Where:
 *      - ΔA is the difference in area
 *      - ΔA0 is the target area difference
 *
 * This class inherits from AbstractGlobalCurvature.
 */
#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
#include "AbstractGlobalCurvature.h"

class CouplingtoFixedGlobalCurvature : public  AbstractGlobalCurvature {
    
public:
    CouplingtoFixedGlobalCurvature(VAHGlobalMeshProperties *VHA, State *pstate, double Gkappa, double GlobalC0);
    ~CouplingtoFixedGlobalCurvature();

       inline double GetEnergy()                        {return m_Energy;} // Only part of energy asscoiated with this class
    
public:
    void Initialize();
    void UpdateEnergyChange(double delta_area, double delta_curvature);
    double CalculateEnergyChange(double delta_area, double delta_curvature);
    
    inline  std::string GetDerivedDefaultReadName()  {return "FixedGlobalCurvature";}
    inline static std::string GetDefaultReadName() {return "FixedGlobalCurvature";}

private:
    double m_K;        // energy coupling constant (in the constructor it will be devided by 2)
    double m_gC0;    // Global curvature
};


#endif
