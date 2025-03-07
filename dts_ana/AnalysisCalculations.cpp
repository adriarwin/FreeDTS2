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

    //Calculating projected area
    if (m_pState->GetAnalysisVariables()->GetProjectedAreaCalculationActive()){
        m_pProjectedArea=(*(m_pState->GetMesh()->GetBox()))(1)*(*(m_pState->GetMesh()->GetBox()))(0);
    }
    //Calculating total energy
    if (m_pState->GetAnalysisVariables()->GetEnergyCalculationActive()){
        double totalE = m_pState->GetEnergyCalculator()->CalculateAllLocalEnergy();
        m_pState->GetEnergyCalculator()->UpdateTotalEnergy(totalE);
        m_pEnergy=m_pState->GetEnergyCalculator()->GetEnergy();
    }

    //Calculating area and volume
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

    //Calculating curvature and thickness
    for (std::vector<vertex *>::iterator it = all_vertex.begin() ; it != all_vertex.end(); ++it) {

        area_v=(*it)->GetArea();
        area_t+=area_v;
        
        //Mean Curvature loop
        if (m_pState->GetAnalysisVariables()->GetMeanCurvatureCalculationActive()) {
            m_pMeanCurvature+= area_v*0.5*((*it)->GetP1Curvature() + (*it)->GetP2Curvature());

        }
        //Gaussian Curvature loop
        if (m_pState->GetAnalysisVariables()->GetGaussianCurvatureCalculationActive()) {
            m_pGaussianCurvature+= area_v*((*it)->GetP1Curvature()*(*it)->GetP2Curvature());
        }
        //Thickness loop
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

    if (m_pState->GetAnalysisVariables()->GetInclusionCalculationActive()) {
        CalculateInclusion();
    }


    /*if (m_pState->GetAnalysisVariables()->GetMeanCurvatureCalculationActive()) {
            m_pMeanCurvature/=area_t;

        }
        if (m_pState->GetAnalysisVariables()->GetGaussianCurvatureCalculationActive()) {
            m_pGaussianCurvature/= area_t;
        }*/

}

void AnalysisCalculations::CalculateInclusion(){

    std::vector<inclusion *>& all_inclusions = m_pState->GetMesh()->GetInclusion();

    double area_v;

    for (std::vector<inclusion *>::iterator it = all_inclusions.begin() ; it != all_inclusions.end(); ++it){

        //Calculating inclusion neighbours
        std::vector<vertex *> neighbour_list= (*it)->Getvertex()->GetVNeighbourVertex();
        for (std::vector<vertex *>::iterator it2 = neighbour_list.begin() ; it2 != neighbour_list.end(); ++it2){
            if ((*it2)->VertexOwnInclusion()){
                    m_pInclusionNeighbour++;
            }
        }

        area_v=(*it)->Getvertex()->GetArea();
        m_pInclusionMeanCurvature+= area_v*0.5*((*it)->Getvertex()->GetP1Curvature() + (*it)->Getvertex()->GetP2Curvature());
        m_pInclusionGaussianCurvature+= area_v*((*it)->Getvertex()->GetP1Curvature()*(*it)->Getvertex()->GetP2Curvature());
        m_pInclusionEnergy+=(*it)->Getvertex()->GetEnergy();
    }
    
    m_pInclusionNeighbour/=static_cast<double>(all_inclusions.size());

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

        if (m_pState->GetAnalysisVariables()->GetInclusionCalculationActive()) {
            m_pInclusionNeighbour=0;
            m_pInclusionEnergy=0;
            m_pInclusionMeanCurvature=0;
            m_pInclusionGaussianCurvature=0;
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

