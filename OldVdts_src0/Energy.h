#if !defined(AFX_Energy_H_334B21B8_C13C_2248_BF23_124095086233__INCLUDED_)
#define AFX_Energy_H_334B21B8_C13C_2248_BF23_124095086233__INCLUDED_

#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
#include "Inclusion_Interaction_Map.h"

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 This object calculate energy of the system
 */
class Energy
{
public:
    Energy();
	Energy(Inclusion_Interaction_Map * );
	 ~Energy();

     inline Inclusion_Interaction_Map  * GetIncIntMap()                const  {return m_pInt;}
private:
    double m_Kappa;
    double m_KappaG;
    double m_mem_c0;
    
    //== vertex area
    double m_Kva;
    double m_av0;
    double m_KvaEdge;
    double m_av0Edge;

    std::vector<double> m_Membrane_model_parameters;
    int m_NO_Membrane_model_parameters;
    
    //=== edge
    double m_Lambda;
    double m_KnEdge;
    double m_KgEdge;
    
    // field info
     Vec3D m_FieldDirection;
     double m_FieldStrength;

public:
    double TotalEnergy(std::vector<vertex *> pVeretx, std::vector<links *> plink);   ///
    double Energy_OneVertexMove(vertex * pVeretx);   ///
    double Energy_OneLinkFlip(links * pLinks);
    double SingleVertexEnergy(vertex *p);
    double TwoInclusionsInteractionEnergy(links *);
    double InteractionFunction(double N2, double A, double B, double theta);
    double SingleEdgeVertexEnergy(vertex *p);
private:

    Inclusion_Interaction_Map * m_pInt;
    double m_Angle3D;
    double m_Angle2D;
    

private:
    double Geo_Theta(vertex *v1, vertex *v2);
    double F10(vertex *v1, vertex *v2,std::vector<double>);
    double F2(vertex *v1, vertex *v2,std::vector<double>);
    double F11(vertex *v1, vertex *v2,std::vector<double>);




    





};


#endif
