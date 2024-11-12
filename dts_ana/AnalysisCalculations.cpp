#include <stdio.h>
#include "AnalysisCalculations.h"
#include "Nfunction.h"
#include "vertex.h"
#include "State.h"
#include "triangle.h"
#include "links.h"



AnalysisCalculations::AnalysisCalculations(State *pState){
    m_pState=pState;
}

void AnalysisCalculations::Calculate(){


    double area_v;
    double area_t=0.0;
    InitializeMemberVariables();


    if (m_pState->GetAnalysisVariables()->GetProjectedAreaCalculationActive()){
        m_pProjectedArea=(*(m_pState->GetMesh()->GetBox()))(1)*(*(m_pState->GetMesh()->GetBox()))(0);
    }

    if (m_pState->GetAnalysisVariables()->GetEnergyCalculationActive()){
        double totalE = m_pState->GetEnergyCalculator()->CalculateAllLocalEnergy();
        m_pState->GetEnergyCalculator()->UpdateTotalEnergy(totalE);
        m_pEnergy=m_pState->GetEnergyCalculator()->GetEnergy();
    }


    if (m_pState->GetAnalysisVariables()->GetAreaCalculationActive() || m_pState->GetAnalysisVariables()->GetVolumeCalculationActive()) {
        std::vector<triangle *>& all_tri = m_pState->GetMesh()->GetActiveT();
        for (std::vector<triangle *>::iterator it = all_tri.begin() ; it != all_tri.end(); ++it) {
             if (m_pState->GetAnalysisVariables()->GetAreaCalculationActive()){
                m_pArea += (*it)->GetArea();
             }
             if (m_pState->GetAnalysisVariables()->GetVolumeCalculationActive()){
                m_pVolume += CalculateSingleTriangleVolume((*it));}
            }
    }



    

    std::vector<vertex *>& all_vertex = m_pState->GetMesh()->GetActiveV();
    double minHeight=(*all_vertex.begin())->GetVZPos();
    double maxHeight=(*all_vertex.begin())->GetVZPos();

    
    for (std::vector<vertex *>::iterator it = all_vertex.begin() ; it != all_vertex.end(); ++it) {

        area_v=(*it)->GetArea();
        area_t+=area_v;
        
        if (m_pState->GetAnalysisVariables()->GetMeanCurvatureCalculationActive()) {
            m_pMeanCurvature+= area_v*0.5*((*it)->GetP1Curvature() + (*it)->GetP2Curvature());

        }
        if (m_pState->GetAnalysisVariables()->GetGaussianCurvatureCalculationActive()) {
            m_pGaussianCurvature+= area_v*((*it)->GetP1Curvature()*(*it)->GetP2Curvature());
        }

        if (m_pState->GetAnalysisVariables()->GetThicknessCalculationActive()) {
            double zPos = (*it)->GetVZPos(); // Get the Z-coordinate of the vertex
            if (zPos < minHeight) {
                minHeight = zPos; // Update min height
            }
            if (zPos > maxHeight) {
                maxHeight = zPos; // Update max height
            }

        }
    }
    if (m_pState->GetAnalysisVariables()->GetThicknessCalculationActive()) {
        m_pThickness=maxHeight - minHeight;}

    /*if (m_pState->GetAnalysisVariables()->GetMeanCurvatureCalculationActive()) {
            m_pMeanCurvature/=area_t;

        }
        if (m_pState->GetAnalysisVariables()->GetGaussianCurvatureCalculationActive()) {
            m_pGaussianCurvature/= area_t;
        }*/

}

void AnalysisCalculations::InitializeMemberVariables(){

    if (m_pState->GetAnalysisVariables()->GetEnergyCalculationActive()) {
            m_pEnergy=0;
        }
        if (m_pState->GetAnalysisVariables()->GetAreaCalculationActive()) {
            m_pArea=0;
        }
        if (m_pState->GetAnalysisVariables()->GetVolumeCalculationActive()) {
            m_pVolume=0;
        }
        if (m_pState->GetAnalysisVariables()->GetProjectedAreaCalculationActive()) {
            m_pProjectedArea=0;

        }
        if (m_pState->GetAnalysisVariables()->GetMeanCurvatureCalculationActive()) {
            m_pMeanCurvature=0;

        }
        if (m_pState->GetAnalysisVariables()->GetGaussianCurvatureCalculationActive()) {
            m_pGaussianCurvature=0;
        }

        if (m_pState->GetAnalysisVariables()->GetThicknessCalculationActive()) {
            m_pThickness=0;
        }

}


double AnalysisCalculations::CalculateSingleTriangleVolume(triangle *pTriangle){

    /*if(m_pState->GetMesh()->GetHasCrossedPBC()){
        *(m_pState->GetTimeSeriesLog()) << "---> the system has crossed the PBC while volume is being calculated.";
        *(m_pState->GetTimeSeriesLog()) << " SOLUTION: Restart the simulation and center the system. Also, activate the command for centering the box.";

         exit(0);
    }*/
    
    double T_area = pTriangle->GetArea();
    Vec3D Normal_v = pTriangle->GetNormalVector();
    Vec3D Pos = pTriangle->GetV1()->GetPos();

    // Compute triangle volume
    return T_area * (Vec3D::dot(Normal_v, Pos)) / 3.0;
}

