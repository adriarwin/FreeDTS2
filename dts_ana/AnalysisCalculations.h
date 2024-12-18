#ifndef ANALYSISCALCULATIONS_H
#define ANALYSISCALCULATIONS_H

#include <vector>
#include <string>
#include "SimDef.h"
#include "Nfunction.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"

class State;
class AnalysisCalculations {
public:

    AnalysisCalculations(State *pState);

    
    double GetArea() {return m_pArea;}

    double GetVolume() {return m_pVolume;}

    
    double GetEnergy() {return m_pEnergy;}

    
    double GetProjectedArea() {return m_pProjectedArea;}

    
    double GetMeanCurvature() {return m_pMeanCurvature;}

    
    double GetThickness() {return m_pThickness;}

    
    double GetGaussianCurvature() {return m_pGaussianCurvature;}

    void Calculate();

private:

    /*void CalculateArea();
    void CalculateEnergy();
    void CalculateProjectedArea();
    void CalculateMeanCurvature();
    void CalculateThickness();
    void CalculateGaussianCurvature();*/
    void InitializeMemberVariables();
    double CalculateSingleTriangleVolume(triangle *pTriangle);


private:

    double m_pGaussianCurvature;
    double m_pProjectedArea;
    double m_pEnergy;
    double m_pArea;
    double m_pThickness;
    double m_pMeanCurvature;
    double m_pVolume;
    State *m_pState;
    
};

#endif // READTRAJTSIFILES_H