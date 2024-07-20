#if !defined(AFX_CouplingTotalAreaToHarmonicPotential_H_334B21B8_C13C_2248_BF23_124095086233__INCLUDED_)
#define AFX_CouplingTotalAreaToHarmonicPotential_H_334B21B8_C13C_2248_BF23_124095086233__INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 this class try to force the membrane area to a targeted size
*/
#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
#include "TotalAreaCoupling.h"
class CouplingTotalAreaToHarmonicPotential : public  TotalAreaCoupling {
public:
    CouplingTotalAreaToHarmonicPotential(int eqsteps, double gamma,  double P0);
    ~CouplingTotalAreaToHarmonicPotential();



       inline double GetTotalArea()                  {return m_TotalArea;}
       inline bool GetState()                        {return m_State;}
    


public:
    void Initialize(std::vector<triangle *> &pTriangle);   ///
    double CalculateEnergyChange(int step, double oldarea,  double newarea);
    void UpdateArea(double oldarea, double newarea);

    
private:
    double AreaofTrianglesAroundVertex(vertex * pVeretx);   ///
    double m_TotalArea;
    int m_NoEQStep;
    double m_K0;
    double m_A0;
    double m_NT;  // number of the triangles in the system
    bool m_State;
    double m_Gamma;
    double m_6SQPI;   /// 1/6pi^1/2

};


#endif
